"""Pure-Python circular binary segmentation derived from DNAcopy.

DNAcopy's CBS implementation is GPL-licensed and implemented mostly in
Fortran. This module reconstructs the weighted segmentation algorithm in
Python/NumPy so the rebuild project does not need R, DNAcopy, or a compiled
extension at runtime.

Reference implementation:
https://git.bioconductor.org/packages/DNAcopy

The statistical shape follows DNAcopy's ``changepoints.R``,
``changepoints-wtd.f`` and ``cbsWtstats.f``:

* scan circular arcs for the maximum weighted two-sample statistic;
* recursively accept binary or ternary changepoints;
* test candidate splits against a permutation reference distribution;
* re-check internal arc edges with a binary permutation test;
* immediately accept strong, sufficiently wide signals.

DNAcopy's Fortran block-pruning and hybrid tail approximation are performance
optimizations. Here the arc scan is vectorized with NumPy and the permutation
loop stops early once a candidate cannot pass the requested alpha threshold.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import special, stats

if TYPE_CHECKING:
    from numpy import ndarray

    from ..cnary import CopyNumArray

DEFAULT_NPERM = 10000
MIN_WIDTH = 2
KMAX = 25
NMIN = 200
ETA = 0.05
NGRID = 100
TOL = 1e-6
STRONG_T = 7.0
STRONG_MIN_WIDTH = 10
RNG_SEED = 0xA5EED


@dataclass(frozen=True)
class ArcStatistic:
    """Maximum circular-arc statistic and its half-open bin interval."""

    tstat: float
    start: int
    end: int


def _weighted_mean(values: ndarray, weights: ndarray) -> float:
    return float(np.average(values, weights=weights))


def _nu(x: float, tol: float = TOL) -> float:
    """Translate ``nu`` from DNAcopy ``src/tailprobs.f``."""
    if x <= 0.01:
        return float(np.exp(-0.583 * x))
    log_nu = np.log(2.0) - 2.0 * np.log(x)
    previous = log_nu
    block = 2
    k = 0
    while True:
        for _ in range(block):
            k += 1
            log_nu -= 2.0 * stats.norm.cdf(-x * np.sqrt(k) / 2.0) / k
        if abs((log_nu - previous) / log_nu) <= tol:
            return float(np.exp(log_nu))
        previous = log_nu
        block *= 2


def _it1tsq(x: float, width: float) -> float:
    """Translate ``it1tsq`` from DNAcopy ``src/tailprobs.f``."""

    def antiderivative(value: float) -> float:
        shifted = value - 0.5
        return 8.0 * shifted / (1.0 - 4.0 * shifted**2) + 2.0 * np.log(
            (1.0 + 2.0 * shifted) / (1.0 - 2.0 * shifted)
        )

    return float(antiderivative(x + width) - antiderivative(x))


def _tail_probability(
    statistic: float,
    delta: float,
    n: int,
    ngrid: int = NGRID,
    tol: float = TOL,
) -> float:
    """Translate DNAcopy's Siegmund/Yao hybrid tail approximation."""
    increment = (0.5 - delta) / ngrid
    total = 0.0
    left = 0.5 - increment
    midpoint = 0.5 - 0.5 * increment
    scaled = statistic / np.sqrt(n)
    for _ in range(ngrid):
        left += increment
        midpoint += increment
        x = scaled / np.sqrt(midpoint * (1.0 - midpoint))
        total += _nu(float(x), tol) ** 2 * _it1tsq(left, increment)
    return float(2.0 * 9.973557e-2 * statistic**3 * np.exp(-(statistic**2) / 2.0) * total)


def _log_choose(n: float, k: float) -> float:
    return float(special.gammaln(n + 1.0) - special.gammaln(k + 1.0) - special.gammaln(n - k + 1.0))


def _eta_boundary(nperm: int, eta: float, successes: int) -> list[int]:
    """Translate ``etabdry`` from DNAcopy ``src/getbdry.f``."""
    boundaries: list[int] = []
    accepted = 0
    for draw in range(1, nperm + 1):
        probability = stats.hypergeom.cdf(accepted, nperm, successes, draw)
        if probability <= eta:
            accepted += 1
            boundaries.append(draw)
            if accepted == successes:
                break
    return boundaries


def _probability_exceeding_boundary(nperm: int, boundaries: list[int]) -> float:
    """Translate ``pexceed`` from DNAcopy ``src/getbdry.f``."""
    successes = len(boundaries)
    log_total = _log_choose(float(nperm), float(successes))
    probability = np.exp(_log_choose(float(nperm - boundaries[0]), float(successes)) - log_total)
    if successes >= 2:
        probability += np.exp(
            np.log(boundaries[0])
            + _log_choose(float(nperm - boundaries[1]), float(successes - 1))
            - log_total
        )
    if successes >= 3:
        remaining = float(nperm - boundaries[2])
        choose = _log_choose(remaining, float(successes - 2)) - log_total
        probability += np.exp(np.log(boundaries[0]) + np.log(boundaries[0] - 1.0) - np.log(2.0) + choose)
        probability += np.exp(np.log(boundaries[0]) + np.log(boundaries[1] - boundaries[0]) + choose)
    for index in range(3, successes):
        first = float(boundaries[index - 3])
        second = float(boundaries[index - 2])
        third = float(boundaries[index - 1])
        choose = _log_choose(float(nperm - boundaries[index]), float(successes - index)) - log_total
        probability += np.exp(_log_choose(first, float(index)) + choose)
        probability += np.exp(_log_choose(first, float(index - 1)) + np.log(third - first) + choose)
        probability += np.exp(
            _log_choose(first, float(index - 2))
            + np.log(second - first)
            + np.log(third - second)
            + choose
        )
        probability += np.exp(
            _log_choose(first, float(index - 2))
            + np.log(second - first)
            - np.log(2.0)
            + np.log(second - first - 1.0)
            + choose
        )
    return float(probability)


@lru_cache(maxsize=None)
def _sequential_boundary(nperm: int, accepted_rejections: int, eta: float = ETA) -> tuple[int, ...]:
    """Return the DNAcopy sequential-stopping row for one rejection limit."""
    successes = accepted_rejections + 1
    if successes == 1:
        return (nperm - int(nperm * eta),)
    eta0 = eta
    candidate: list[int] = []
    for row_size in range(2, successes + 1):
        low_eta = eta0 * 0.25
        high_eta = eta0 * 1.1
        low = _eta_boundary(nperm, low_eta, row_size)
        high = _eta_boundary(nperm, high_eta, row_size)
        low_probability = _probability_exceeding_boundary(nperm, low)
        high_probability = _probability_exceeding_boundary(nperm, high)
        candidate = low
        while (high_eta - low_eta) / low_eta > 1e-2:
            eta0 = low_eta + (high_eta - low_eta) * (eta - low_probability) / (
                high_probability - low_probability
            )
            candidate = _eta_boundary(nperm, eta0, row_size)
            probability = _probability_exceeding_boundary(nperm, candidate)
            if probability > eta:
                high_eta = eta0
                high_probability = probability
            else:
                low_eta = eta0
                low_probability = probability
    return tuple(candidate)


def _max_tstat(values: ndarray, weights: ndarray, center: bool = True) -> ArcStatistic:
    """Return DNAcopy's maximum weighted circular-arc t statistic.

    A circular split is represented by a contiguous arc and its complement.
    Enumerating half-open arcs ``[start:end]`` with both sides at least
    ``MIN_WIDTH`` bins therefore covers the same candidate split space.
    """
    n = len(values)
    if n < 2 * MIN_WIDTH or np.ptp(values) == 0:
        return ArcStatistic(0.0, 0, n)

    if center:
        values = values - _weighted_mean(values, weights)
    weighted = values * weights
    prefix_sum = np.concatenate(([0.0], np.cumsum(weighted)))
    prefix_weight = np.concatenate(([0.0], np.cumsum(weights)))
    total_weight = prefix_weight[-1]
    tss = float(np.sum(weights * values**2) - np.sum(weighted) ** 2 / total_weight)
    best_bss = 0.0
    best_start = 0
    best_end = n

    for width in range(MIN_WIDTH, n - MIN_WIDTH + 1):
        arc_sums = prefix_sum[width:] - prefix_sum[:-width]
        arc_weights = prefix_weight[width:] - prefix_weight[:-width]
        denom = arc_weights * (total_weight - arc_weights)
        valid = denom > 0
        if not valid.any():
            continue
        bss = np.zeros_like(arc_sums)
        bss[valid] = arc_sums[valid] ** 2 * total_weight / denom[valid]
        idx = int(np.argmax(bss))
        if bss[idx] > best_bss:
            best_bss = float(bss[idx])
            best_start = idx
            best_end = idx + width

    if tss <= best_bss + 0.0001:
        tss = best_bss + 1.0
    return ArcStatistic(best_bss / ((tss - best_bss) / (n - 2.0)), best_start, best_end)


def _max_small_arc_tstat(values: ndarray, weights: ndarray, kmax: int) -> float:
    """Translate the weighted small-arc statistic ``hwtmaxp``."""
    n = len(values)
    total_weight = float(weights.sum())
    tss = float(np.sum(weights * values**2) - np.sum(values * weights) ** 2 / total_weight)
    circular_values = np.concatenate((values, values))
    circular_weights = np.concatenate((weights, weights))
    prefix_sum = np.concatenate(([0.0], np.cumsum(circular_values * circular_weights)))
    prefix_weight = np.concatenate(([0.0], np.cumsum(circular_weights)))
    best_bss = 0.0
    for width in range(MIN_WIDTH, min(kmax, n - MIN_WIDTH) + 1):
        arc_sums = prefix_sum[width : n + width] - prefix_sum[:n]
        arc_weights = prefix_weight[width : n + width] - prefix_weight[:n]
        denom = arc_weights * (total_weight - arc_weights)
        valid = denom > 0
        if valid.any():
            best_bss = max(best_bss, float(np.max(arc_sums[valid] ** 2 * total_weight / denom[valid])))
    if tss <= best_bss + 0.0001:
        tss = best_bss + 1.0
    return best_bss / ((tss - best_bss) / (n - 2.0))


def _minimum_arc_weight_fraction(weights: ndarray, width: int) -> float:
    """Translate the ``delta`` result from DNAcopy ``getmncwt``."""
    total = float(weights.sum())
    circular = np.concatenate((weights, weights))
    prefix = np.concatenate(([0.0], np.cumsum(circular)))
    direct = prefix[width : len(weights) + width] - prefix[: len(weights)]
    return float(np.min(np.minimum(direct, total - direct)) / total)


def _permute_weighted(values: ndarray, weights: ndarray, rng: np.random.Generator) -> ndarray:
    """Translate DNAcopy's weighted Fisher-Yates routine ``wxperm``."""
    roots = np.sqrt(weights)
    permuted = values * roots
    for index in range(len(permuted) - 1, -1, -1):
        selected = int(rng.integers(0, index + 1))
        held = permuted[index]
        permuted[index] = permuted[selected] / roots[index]
        permuted[selected] = held
    return permuted


def _candidate_passes(
    values: ndarray,
    weights: ndarray,
    observed: ArcStatistic,
    alpha: float,
    nperm: int,
    rng: np.random.Generator,
) -> bool:
    """Translate weighted ``wfindcpt`` acceptance and early stopping."""
    minor_width = min(observed.end - observed.start, len(values) - observed.end + observed.start)
    if observed.tstat**0.5 <= 0.1:
        return False
    if observed.tstat**0.5 >= STRONG_T and minor_width >= STRONG_MIN_WIDTH:
        return True

    hybrid = len(values) > NMIN
    if hybrid:
        delta = _minimum_arc_weight_fraction(weights, KMAX + 1)
        approximation = _tail_probability(observed.tstat**0.5, delta, len(values))
        if approximation > alpha:
            return False
        rejection_limit = int((alpha - approximation) * nperm)
    else:
        rejection_limit = int(alpha * nperm)
    boundary = _sequential_boundary(nperm, rejection_limit)
    rejected = 0
    observed_tstat = observed.tstat * 0.99999
    for draw in range(1, nperm + 1):
        permuted = _permute_weighted(values, weights, rng)
        if hybrid:
            permutation_tstat = _max_small_arc_tstat(permuted, weights, KMAX)
        else:
            permutation_tstat = _max_tstat(permuted, weights, center=False).tstat
        if permutation_tstat >= observed_tstat:
            rejected += 1
            if rejected > rejection_limit:
                return False
        if draw >= boundary[rejected]:
            return True
    return True


def _binary_permutation_pvalue(
    values: ndarray,
    weights: ndarray,
    split: int,
    nperm: int,
    rng: np.random.Generator,
) -> float:
    """Re-check one internal arc edge as DNAcopy's ``wtpermp`` does."""
    n = len(values)
    if split <= 1 or split >= n - 1:
        return 1.0
    left_weight = float(weights[:split].sum())
    right_weight = float(weights[split:].sum())
    total_weight = left_weight + right_weight
    mean = _weighted_mean(values, weights)
    left_delta = abs(_weighted_mean(values[:split], weights[:split]) - mean)
    right_delta = abs(_weighted_mean(values[split:], weights[split:]) - mean)
    if split <= n - split:
        selected = split
        selected_weight = left_weight
        observed = left_delta
    else:
        selected = n - split
        selected_weight = right_weight
        observed = right_delta

    tss = float(np.sum(weights * (values - mean) ** 2))
    explained = observed**2 * selected_weight * total_weight / (total_weight - selected_weight)
    tstat = explained / max((tss - explained) / (n - 2.0), np.finfo(float).eps)
    if tstat > 25 and selected >= STRONG_MIN_WIDTH:
        return 0.0

    roots = np.sqrt(weights)
    permuted = values.copy()
    permuted[:split] *= roots[:split]
    rejected = 0
    for _ in range(nperm):
        sample_sum = 0.0
        for index in range(n - 1, n - selected - 1, -1):
            chosen = int(rng.integers(0, index + 1))
            held = permuted[index]
            permuted[index] = permuted[chosen]
            permuted[chosen] = held
            sample_sum += permuted[index] * roots[index]
        sampled_delta = abs(sample_sum / selected_weight - mean)
        if sampled_delta >= observed * 0.99999:
            rejected += 1
    return rejected / nperm


def _find_changepoints(
    values: ndarray,
    weights: ndarray,
    alpha: float,
    nperm: int,
    rng: np.random.Generator,
) -> tuple[int, ...]:
    """Return zero, one, or two DNAcopy-style breakpoints for one segment."""
    if len(values) < 2 * MIN_WIDTH or np.ptp(values) == 0:
        return ()
    centered = values - _weighted_mean(values, weights)
    observed = _max_tstat(centered, weights)
    if not _candidate_passes(centered, weights, observed, alpha, nperm, rng):
        return ()
    if observed.end == len(values):
        return (observed.start,)
    if observed.start == 0:
        return (observed.end,)

    changepoints: list[int] = []
    if _binary_permutation_pvalue(
        centered[: observed.end],
        weights[: observed.end],
        observed.start,
        nperm,
        rng,
    ) <= alpha:
        changepoints.append(observed.start)
    if _binary_permutation_pvalue(
        centered[observed.start :],
        weights[observed.start :],
        observed.end - observed.start,
        nperm,
        rng,
    ) <= alpha:
        changepoints.append(observed.end)
    return tuple(changepoints)


def _recursive_breakpoints(
    values: ndarray,
    weights: ndarray,
    alpha: float,
    nperm: int,
    rng: np.random.Generator,
) -> list[int]:
    """Recursively segment one chromosome arm and return bin boundaries."""

    def visit(lo: int, hi: int) -> list[int]:
        local = _find_changepoints(values[lo:hi], weights[lo:hi], alpha, nperm, rng)
        if not local:
            return [hi]
        boundaries: list[int] = []
        start = lo
        for breakpoint in local:
            end = lo + breakpoint
            boundaries.extend(visit(start, end))
            start = end
        boundaries.extend(visit(start, hi))
        return boundaries

    return visit(0, len(values))


def segment_pycbs(cnarr: CopyNumArray, alpha: float) -> CopyNumArray:
    """Segment one CNVkit chromosome arm with the Python CBS reconstruction."""
    if not len(cnarr):
        return cnarr
    nperm = int(os.environ.get("CNVKIT_PYCBS_NPERM", DEFAULT_NPERM))
    if nperm < 1:
        raise ValueError("CNVKIT_PYCBS_NPERM must be a positive integer")
    values = cnarr["log2"].to_numpy(dtype=float)
    if "weight" in cnarr:
        weights = cnarr["weight"].to_numpy(dtype=float)
    else:
        weights = np.ones(len(cnarr), dtype=float)
    rng = np.random.default_rng(RNG_SEED)
    breakpoints = _recursive_breakpoints(values, weights, alpha, nperm, rng)
    chrom = cnarr.chromosome.iat[0]
    rows: list[tuple[object, ...]] = []
    lo = 0
    for hi in breakpoints:
        rows.append(
            (
                chrom,
                cnarr.start.iat[lo],
                cnarr.end.iat[hi - 1],
                _weighted_mean(values[lo:hi], weights[lo:hi]),
                "-",
                hi - lo,
            )
        )
        lo = hi
    logging.debug("Python CBS segmented %s into %d regions", chrom, len(rows))
    table = pd.DataFrame.from_records(
        rows,
        columns=["chromosome", "start", "end", "log2", "gene", "probes"],
    )
    segarr = cnarr.as_dataframe(table)
    segarr.sort_columns()
    return segarr
