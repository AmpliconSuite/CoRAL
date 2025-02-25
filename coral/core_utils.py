from __future__ import annotations

from dataclasses import dataclass
import functools
import pathlib
import time
from typing import TYPE_CHECKING, Callable, Generic, NamedTuple, ParamSpec, TypeVar

import memray

from coral.datatypes import SolverOptions

if TYPE_CHECKING:
    from coral.datatypes import EdgeId, Node, Walk
    from coral.breakpoint.breakpoint_graph import BreakpointGraph
P = ParamSpec("P")
T = TypeVar("T")


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

class ProfileResult(NamedTuple,Generic[T]):
    peak_ram_gb: float | None
    runtime_s: float | None
    result: T

def profile_fn(fn: Callable[P, T]) -> Callable[P, T]:
    @functools.wraps(fn)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        if not kwargs.get("should_profile"):
            return fn(*args, **kwargs)
        
        if not kwargs.get("bp_graph") or not kwargs.get("solver_options"):
            raise ValueError("bp_graph and solver_options must be provided")

        solver_options: SolverOptions = kwargs.get("solver_options") # type: ignore[assignment]
        bp_graph: BreakpointGraph = kwargs.get("bp_graph") # type: ignore[assignment]

        profile_path = f"{solver_options.output_dir}" \
            f"/amplicon_{bp_graph.amplicon_idx}_mem_profile.bin"
        # Tracker errors if file already exists
        if pathlib.Path(profile_path).exists():
            pathlib.Path(profile_path).unlink()

        start = time.perf_counter()
        try:
            with memray.Tracker(profile_path):
                result = fn(*args, **kwargs)
            mem_stats = memray._memray.compute_statistics(profile_path)
            bp_graph.peak_ram_gb = mem_stats.peak_memory_allocated / 1e9
            end = time.perf_counter()
            bp_graph.runtime_s = end - start
            return result
        except:
            pass # Avoid crashing if private API changes

        return fn(*args, **kwargs)
    return wrapper
