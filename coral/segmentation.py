#!/usr/bin/env python3
"""Segment CNVkit .cnr files and write .cns output with Python CBS.

The preferred backend is a full pure-Python ``pycbs`` reconstruction of the
DNAcopy CBS algorithm. This script first uses a CNVkit installation that already
registers ``method="pycbs"``; if that is unavailable, it uses the bundled CoRAL
``coral.pycbs`` module in this branch. A simple weighted binary-segmentation
fallback remains available for emergency use only.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


LOGGER = logging.getLogger("segmentation")
RNG_SEED = 0xA5EED


@dataclass(frozen=True)
class Segment:
    chromosome: str
    start: int
    end: int
    gene: str
    log2: float
    depth: float
    probes: int
    weight: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Python CBS segmentation on a CNVkit .cnr file.",
    )
    parser.add_argument("cnr", type=Path, help="Input CNVkit .cnr file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output .cns path; default replaces .cnr suffix with .cns",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.0001,
        help="CBS significance threshold, default: 0.0001",
    )
    parser.add_argument(
        "--nperm",
        type=int,
        default=1000,
        help="Permutation count for fallback CBS, default: 1000",
    )
    parser.add_argument(
        "--min-probes",
        type=int,
        default=2,
        help="Minimum probes on each side of a fallback split, default: 2",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Worker processes for CNVkit pycbs backend, default: 1",
    )
    parser.add_argument(
        "--backend",
        choices=("auto", "cnvkit-pycbs", "fallback"),
        default="auto",
        help="Segmentation backend, default: auto",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress messages to stderr",
    )
    return parser.parse_args()


def run_cnvkit_pycbs(
    cnr_path: Path,
    output_path: Path,
    alpha: float,
    processes: int,
) -> bool:
    try:
        from skgenome import tabio
    except Exception as exc:
        LOGGER.info("CNVkit tab I/O is not importable: %s", exc)
        return False

    cnarr = tabio.read(str(cnr_path), "tab")
    try:
        cnarr.sample_id = cnr_path.stem
    except AttributeError:
        pass

    try:
        from cnvlib.segmentation import SEGMENT_METHODS, do_segmentation

        if "pycbs" in SEGMENT_METHODS:
            segarr = do_segmentation(
                cnarr,
                method="pycbs",
                threshold=alpha,
                processes=processes,
            )
            tabio.write(segarr, str(output_path), "tab", verbose=False)
            return True
        LOGGER.info("CNVkit installation does not expose method 'pycbs'")
    except Exception as exc:
        LOGGER.info("CNVkit registered pycbs backend is not usable: %s", exc)

    try:
        try:
            from coral import pycbs as bundled_pycbs
        except Exception:
            import pycbs as bundled_pycbs  # type: ignore[no-redef]

        segarr = bundled_pycbs.segment_pycbs(cnarr, alpha)
    except Exception as exc:
        LOGGER.info("Bundled CoRAL pycbs backend is not usable: %s", exc)
        return False

    tabio.write(segarr, str(output_path), "tab", verbose=False)
    return True


def clean_cnr(df: pd.DataFrame) -> pd.DataFrame:
    required = {"chromosome", "start", "end", "log2"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required .cnr columns: {sorted(missing)}")

    df = df.copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["log2"] = pd.to_numeric(df["log2"], errors="coerce")
    if "depth" in df:
        df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    else:
        df["depth"] = np.nan
    if "weight" in df:
        df["weight"] = pd.to_numeric(df["weight"], errors="coerce")
    else:
        df["weight"] = 1.0

    keep = (
        df["chromosome"].notna()
        & np.isfinite(df["start"])
        & np.isfinite(df["end"])
        & np.isfinite(df["log2"])
        & np.isfinite(df["weight"])
        & (df["weight"] > 0)
    )
    df = df.loc[keep].copy()
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    return df.sort_values(["chromosome", "start", "end"], kind="mergesort")


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(values, weights=weights))


def best_split(
    values: np.ndarray,
    weights: np.ndarray,
    min_probes: int,
) -> tuple[int | None, float]:
    n = len(values)
    if n < min_probes * 2 or np.ptp(values) == 0:
        return None, 0.0

    prefix_w = np.concatenate(([0.0], np.cumsum(weights)))
    prefix_x = np.concatenate(([0.0], np.cumsum(values * weights)))
    total_w = prefix_w[-1]
    total_x = prefix_x[-1]
    mean = total_x / total_w
    total_ss = float(np.sum(weights * (values - mean) ** 2))
    best_idx: int | None = None
    best_stat = 0.0

    for idx in range(min_probes, n - min_probes + 1):
        left_w = prefix_w[idx]
        right_w = total_w - left_w
        if left_w <= 0 or right_w <= 0:
            continue
        left_mean = prefix_x[idx] / left_w
        right_mean = (total_x - prefix_x[idx]) / right_w
        between = left_w * right_w / total_w * (left_mean - right_mean) ** 2
        within = max(total_ss - between, np.finfo(float).eps)
        stat = between / (within / max(n - 2, 1))
        if stat > best_stat:
            best_stat = float(stat)
            best_idx = idx

    return best_idx, best_stat


def permutation_pvalue(
    values: np.ndarray,
    weights: np.ndarray,
    observed: float,
    min_probes: int,
    nperm: int,
    alpha: float,
    rng: np.random.Generator,
) -> float:
    if nperm < 1:
        return 1.0
    exceed = 0
    for draw in range(1, nperm + 1):
        permuted = rng.permutation(values)
        _, stat = best_split(permuted, weights, min_probes)
        if stat >= observed * 0.99999:
            exceed += 1
            if exceed / draw > alpha and draw >= min(100, nperm):
                break
    return (exceed + 1.0) / (draw + 1.0)


def recursive_boundaries(
    values: np.ndarray,
    weights: np.ndarray,
    alpha: float,
    nperm: int,
    min_probes: int,
    rng: np.random.Generator,
) -> list[int]:
    def visit(lo: int, hi: int) -> list[int]:
        local_values = values[lo:hi]
        local_weights = weights[lo:hi]
        split, stat = best_split(local_values, local_weights, min_probes)
        if split is None:
            return [hi]
        pvalue = permutation_pvalue(
            local_values,
            local_weights,
            stat,
            min_probes,
            nperm,
            alpha,
            rng,
        )
        if pvalue > alpha:
            return [hi]
        mid = lo + split
        return visit(lo, mid) + visit(mid, hi)

    return visit(0, len(values))


def summarize_gene(series: pd.Series) -> str:
    genes: list[str] = []
    for value in series.dropna().astype(str):
        for gene in value.split(","):
            gene = gene.strip()
            if gene and gene != "Background" and gene not in genes:
                genes.append(gene)
            if len(genes) >= 50:
                return ",".join(genes)
    return ",".join(genes) if genes else "-"


def build_segments(
    df: pd.DataFrame,
    alpha: float,
    nperm: int,
    min_probes: int,
) -> list[Segment]:
    rng = np.random.default_rng(RNG_SEED)
    segments: list[Segment] = []

    for chrom, chrom_df in df.groupby("chromosome", sort=False):
        values = chrom_df["log2"].to_numpy(dtype=float)
        weights = chrom_df["weight"].to_numpy(dtype=float)
        boundaries = recursive_boundaries(
            values,
            weights,
            alpha,
            nperm,
            min_probes,
            rng,
        )
        lo = 0
        for hi in boundaries:
            part = chrom_df.iloc[lo:hi]
            part_weights = part["weight"].to_numpy(dtype=float)
            part_depth = part["depth"].to_numpy(dtype=float)
            finite_depth = np.isfinite(part_depth)
            if finite_depth.any():
                depth = weighted_mean(part_depth[finite_depth], part_weights[finite_depth])
            else:
                depth = math.nan
            segments.append(
                Segment(
                    chromosome=str(chrom),
                    start=int(part["start"].iloc[0]),
                    end=int(part["end"].iloc[-1]),
                    gene=summarize_gene(part["gene"]) if "gene" in part else "-",
                    log2=weighted_mean(part["log2"].to_numpy(dtype=float), part_weights),
                    depth=depth,
                    probes=int(len(part)),
                    weight=float(part_weights.sum()),
                )
            )
            lo = hi

    return segments


def run_fallback(
    cnr_path: Path,
    output_path: Path,
    alpha: float,
    nperm: int,
    min_probes: int,
) -> None:
    raw = pd.read_csv(cnr_path, sep="\t")
    df = clean_cnr(raw)
    if df.empty:
        pd.DataFrame(
            columns=["chromosome", "start", "end", "gene", "log2", "depth", "probes", "weight"],
        ).to_csv(output_path, sep="\t", index=False)
        return
    segments = build_segments(df, alpha, nperm, min_probes)
    pd.DataFrame([seg.__dict__ for seg in segments]).to_csv(
        output_path,
        sep="\t",
        index=False,
        float_format="%.6g",
    )


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    cnr_path = args.cnr.resolve()
    output_path = args.output or cnr_path.with_suffix(".cns")
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.backend in {"auto", "cnvkit-pycbs"}:
        ok = run_cnvkit_pycbs(cnr_path, output_path, args.alpha, args.processes)
        if ok:
            LOGGER.info("Wrote %s using CNVkit pycbs backend", output_path)
            return 0
        if args.backend == "cnvkit-pycbs":
            raise SystemExit("Requested backend cnvkit-pycbs is unavailable")

    run_fallback(cnr_path, output_path, args.alpha, args.nperm, args.min_probes)
    LOGGER.info("Wrote %s using fallback CBS backend", output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
