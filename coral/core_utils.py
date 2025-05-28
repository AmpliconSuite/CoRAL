from __future__ import annotations

import functools
import logging
import pathlib
import re
import tempfile
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, NamedTuple, ParamSpec, TypeVar
from collections import defaultdict

import colorama
import memray
import numpy as np
import pandas as pd
import pyranges

from coral import datatypes, global_state

if TYPE_CHECKING:
    from coral.datatypes import EdgeId, Node, Walk


P = ParamSpec("P")
T = TypeVar("T")

logger = logging.getLogger(__name__)


def path_to_edge_count(
    path: Walk,
) -> dict[EdgeId, int]:
    edge_counts: dict[EdgeId, int] = defaultdict(int)
    for edge_idx in range(len(path)):
        # Only track edge counts, not nodes
        if edge_idx % 2 == 0:
            edge_counts[path[edge_idx]] += 1  # type: ignore[index]
    return edge_counts


def path_to_str(path: Walk, edge_counts: dict[EdgeId, int]) -> str:
    """Convert path to string for outputting to file (e.g.,
    *_graph.txt, *_cycle.txt).

        ex:
            path = [
                EdgeId(SEQUENCE, 11), chr2-171,805,866[+],
                EdgeId(DISCORDANT, 5), chr2-180,359,131[-],
                EdgeId(SEQUENCE, 15), chr2-180,409,353[+],
                EdgeId(DISCORDANT, 6), chr2-186,985,081[-],
                EdgeId(SEQUENCE, 18)
            ]
        "12+,6-,16+,7-,19+"
    """
    # Reverse path if start node's pos is > than end node's
    if path[0][1] > path[-1][1]:
        path = path[::-1]

    path_str = ""
    for i in range(len(path)):
        if i % 2 == 0:
            edge_id: EdgeId = path[i]  # type: ignore[assignment]
            edge_count = edge_counts[edge_id]
            path_str += f"{edge_id.type.value}{edge_id.idx+1}"
            if i < len(path) - 1:
                next_node: Node = path[i + 1]  # type: ignore[assignment]
                path_str += f"{next_node.strand.value}:{edge_count},"
            else:
                prev_node: Node = path[i - 1]  # type: ignore[assignment]
                path_str += f"{prev_node.strand.inverse.value}:{edge_count}"
    return path_str


def path_to_str__old(path: Walk) -> str:
    """(DEPRECATED) Convert path to string for outputting to file
    (e.g., *_graph.txt, *_cycle.txt).

        ex:
            path = [
                EdgeId(SEQUENCE, 11), chr2-171,805,866[+],
                EdgeId(DISCORDANT, 5), chr2-180,359,131[-],
                EdgeId(SEQUENCE, 15), chr2-180,409,353[+],
                EdgeId(DISCORDANT, 6), chr2-186,985,081[-],
                EdgeId(SEQUENCE, 18)
            ]
            "12+,16+,19+"
    """
    # Reverse path if start node's pos is > than end node's
    if path[0][1] > path[-1][1]:
        path = path[::-1]

    path_str = ""
    for i in range(len(path)):
        if i % 4 == 0:
            edge_id: EdgeId = path[i]  # type: ignore[assignment]
            path_str += f"{edge_id.idx+1}"
            if i < len(path) - 1:
                next_node: Node = path[i + 1]  # type: ignore[assignment]
                path_str += f"{next_node.strand.value},"
            else:
                prev_node: Node = path[i - 1]  # type: ignore[assignment]
                path_str += f"{prev_node.strand.inverse.value}\t"
    return path_str


def profile_fn_with_call_counter(fn: Callable[P, T]) -> Callable[P, T]:
    call_ctr = 0

    @functools.wraps(fn)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        if not global_state.STATE_PROVIDER.should_profile:
            return fn(*args, **kwargs)

        nonlocal call_ctr
        call_ctr += 1
        fn_call = datatypes.FnCall(fn.__name__, call_ctr)

        start = time.perf_counter()
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                profile_path = f"{temp_dir}/{fn.__name__}_profile.bin"
                with memray.Tracker(profile_path):
                    result = fn(*args, **kwargs)
                mem_stats = memray._memray.compute_statistics(profile_path)
                runtime_s = time.perf_counter() - start
                global_state.PROFILED_FN_CALLS[fn_call] = (
                    datatypes.ProfileResult(
                        peak_ram_gb=mem_stats.peak_memory_allocated / 1e9,
                        runtime_s=runtime_s,
                    )
                )
                print(
                    f"{colorama.Style.DIM}{colorama.Fore.LIGHTMAGENTA_EX}"
                    f"Profile result for {fn_call}: "
                    f"{global_state.PROFILED_FN_CALLS[fn_call]}"
                    f"{colorama.Style.RESET_ALL}"
                )
                return result
        except AttributeError:  # Avoid crashing if private API changes
            logger.error("Failed to profile function due to memray error")

        return fn(*args, **kwargs)

    return wrapper


def profile_fn(fn: Callable[P, T]) -> Callable[P, T]:
    @functools.wraps(fn)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        if not global_state.STATE_PROVIDER.should_profile:
            return fn(*args, **kwargs)

        fn_call = datatypes.FnCall(fn.__name__, None)

        start = time.perf_counter()
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                profile_path = f"{temp_dir}/{fn.__name__}_profile.bin"
                with memray.Tracker(profile_path):
                    result = fn(*args, **kwargs)
                mem_stats = memray._memray.compute_statistics(profile_path)
                global_state.PROFILED_FN_CALLS[fn_call] = (
                    datatypes.ProfileResult(
                        peak_ram_gb=mem_stats.peak_memory_allocated / 1e9,
                        runtime_s=time.perf_counter() - start,
                    )
                )
                return result
        except AttributeError:  # Avoid crashing if private API changes
            logger.error("Failed to profile function due to memray error")

        return fn(*args, **kwargs)

    return wrapper


def get_amplicon_id_from_filename(filename: str) -> int:
    return int(filename.split("_")[-2].split("amplicon")[1])


class ReconstructionPaths(NamedTuple):
    graph_path: pathlib.Path
    cycle_path: pathlib.Path | None = None


def get_reconstruction_paths_from_shared_dir(
    reconstruction_dir: pathlib.Path,
) -> list[ReconstructionPaths]:
    reconstruction_paths = []
    for bp_filepath in reconstruction_dir.glob("*_graph.txt"):
        cycle_filename = re.sub(
            r"amplicon(\d+)_graph.txt",
            r"amplicon\1_cycles.txt",
            bp_filepath.name,
        )
        cycle_filepath = reconstruction_dir / cycle_filename
        reconstruction_paths.append(
            ReconstructionPaths(
                graph_path=bp_filepath,
                cycle_path=cycle_filepath if cycle_filepath.exists() else None,
            )
        )
    return sorted(reconstruction_paths, key=lambda x: x.graph_path.name)


def get_reconstruction_paths_from_separate_dirs(
    cycle_dir: pathlib.Path,
    graph_dir: pathlib.Path,
) -> list[ReconstructionPaths]:
    reconstruction_paths = []
    for bp_filepath in graph_dir.glob("*_graph.txt"):
        cycle_filename = re.sub(
            r"amplicon(\d+)_graph.txt",
            r"amplicon\1_cycles.txt",
            bp_filepath.name,
        )
        cycle_filepath = cycle_dir / cycle_filename
        reconstruction_paths.append(
            ReconstructionPaths(
                graph_path=bp_filepath,
                cycle_path=cycle_filepath if cycle_filepath.exists() else None,
            )
        )
    return sorted(reconstruction_paths, key=lambda x: x.graph_path.name)


def bin_genome(
    true_segments: pd.DataFrame,
    reconstructed_segments: pd.DataFrame,
    margin_size: int = 10000,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Bin the intervals by the breakpoints union.
    Warning: Setting a margin size can effect the output of different distances
    """

    df_bins = pd.DataFrame(
        np.concatenate(
            (
                true_segments[["Chromosome", "Start"]].values,
                true_segments[["Chromosome", "End"]].values,
                reconstructed_segments[["Chromosome", "Start"]].values,
                reconstructed_segments[["Chromosome", "End"]].values,
            ),
            axis=0,
        )
    )

    df_bins.columns = ["Chromosome", "Start"]
    # df_bins['Start'] = df_bins['Start']
    df_bins = df_bins.drop_duplicates().sort_values(by=["Chromosome", "Start"])

    # rotate with 1 up the start column
    df_bins_suffix = df_bins.tail(-1)
    df_bins_suffix = pd.concat([df_bins_suffix, df_bins.head(1)])

    df_bins.reset_index(drop=True, inplace=True)
    df_bins_suffix.reset_index(drop=True, inplace=True)
    df_bins = pd.concat([df_bins, df_bins_suffix], axis=1, ignore_index=True)
    df_bins.columns = ["Chromosome", "Start", "Chromosome2", "End"]

    # keep only rows with same chr and non-negative distance
    df_bins = df_bins[
        (df_bins["Chromosome"] == df_bins["Chromosome2"])
        & (abs(df_bins["Start"] - df_bins["End"]) >= 0)
    ]

    df_bins = df_bins[(df_bins["End"] - df_bins["Start"]) > 0]
    df_bins["length"] = abs(df_bins["End"] - df_bins["Start"])
    # chrlist = df_bins["Chromosome"].drop_duplicates().tolist()

    return df_bins[["Chromosome", "Start", "End", "length"]]  # , chrlist
