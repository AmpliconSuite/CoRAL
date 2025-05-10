from __future__ import annotations

import logging

from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.output.traversal import eulerian_cycle_t, eulerian_path_t

logger = logging.getLogger(__name__)


def get_single_cycle_str(
    bp_graph: BreakpointGraph,
    cycle_idx: int,
    weight_sorted_cycle_idx: int,
    *,
    output_path_constraints: bool = True,
) -> str:
    """Generate the output string for a given cycle produced by optimization."""

    cycle_weights = bp_graph.walk_weights.cycles
    logger.debug(f"Traversing next cycle, CN = {cycle_weights[cycle_idx]}")

    satisfied_pc = bp_graph.path_constraints_satisfied.cycles[cycle_idx]

    path_constraints_satisfied_cycle = []
    path_constraints_support_cycle = []
    for pathi in satisfied_pc:
        pathi_ = bp_graph.longest_path_constraints[pathi].pc_idx
        path_constraints_satisfied_cycle.append(
            bp_graph.path_constraints[pathi_].path
        )
        path_constraints_support_cycle.append(
            bp_graph.longest_path_constraints[pathi].support,
        )
    cycle_seg_list = eulerian_cycle_t(
        bp_graph,
        bp_graph.walks.cycles[cycle_idx],
        path_constraints_satisfied_cycle,
        path_constraints_support_cycle,
    )
    assert cycle_seg_list[0] == cycle_seg_list[-1]

    output_str = f"Cycle={weight_sorted_cycle_idx};"
    output_str += f"Copy_count={cycle_weights[cycle_idx]};"
    output_str += "Segments="
    for segi in range(len(cycle_seg_list) - 2):
        output_str += f"{cycle_seg_list[segi][0]}{cycle_seg_list[segi][1]},"
    output_str += f"{cycle_seg_list[-2][0]}{cycle_seg_list[-2][1]}"

    if not output_path_constraints:
        output_str += "\n"
        return output_str

    output_str += ";Path_constraints_satisfied="
    for pathi in range(len(satisfied_pc) - 1):
        output_str += f"{satisfied_pc[pathi]+1},"

    # Avoid trailing comma
    if len(satisfied_pc) > 0:
        output_str += f"{satisfied_pc[-1]+1}\n"
    else:
        output_str += "\n"
    return output_str


def get_single_path_str(
    bp_graph: BreakpointGraph,
    path_idx: int,
    weight_sorted_path_idx: int,
    *,
    output_path_constraints: bool = True,
) -> str:
    """Generate the output string for a given path produced by optimization."""
    if not bp_graph.walks.paths[path_idx]:
        logger.error(f"Path {path_idx} has no edges.")
        return ""

    path_weights = bp_graph.walk_weights.paths
    logger.debug(f"Traversing next path, CN = {path_weights[path_idx]}")

    satisfied_pc = bp_graph.path_constraints_satisfied.paths[path_idx]

    path_constraints_satisfied_path = []
    path_constraints_support_path = []
    for pathi in satisfied_pc:
        pathi_ = bp_graph.longest_path_constraints[pathi].pc_idx
        path_constraints_satisfied_path.append(
            bp_graph.path_constraints[pathi_].path,
        )
        path_constraints_support_path.append(
            bp_graph.longest_path_constraints[pathi].support,
        )
    path_seg_list = eulerian_path_t(
        bp_graph,
        bp_graph.walks.paths[path_idx],
        path_constraints_satisfied_path,
        path_constraints_support_path,
    )
    output_str = f"Path={weight_sorted_path_idx};"
    output_str += f"Copy_count={path_weights[path_idx]};"
    output_str += "Segments=0+,"
    for segi in range(len(path_seg_list) - 1):
        output_str += f"{path_seg_list[segi][0]}{path_seg_list[segi][1]},"
    output_str += f"{path_seg_list[-1][0]}{path_seg_list[-1][1]},0-"

    if not output_path_constraints:
        output_str += "\n"
        return output_str

    output_str += ";Path_constraints_satisfied="
    for pathi in range(len(satisfied_pc) - 1):
        output_str += f"{satisfied_pc[pathi]+1},"

    # Avoid trailing comma
    if len(satisfied_pc) > 0:
        output_str += f"{satisfied_pc[-1]+1}\n"
    else:
        output_str += "\n"

    return output_str
