import io
import logging

from coral.breakpoint import infer_breakpoint_graph
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.models.output import eulerian_cycle_t, eulerian_path_t

logger = logging.getLogger(__name__)


def get_single_cycle_output(
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


def get_single_path_output(
    bp_graph: BreakpointGraph,
    path_idx: int,
    weight_sorted_path_idx: int,
    *,
    output_all_paths: bool = True,
) -> str:
    """Generate the output string for a given path produced by optimization."""

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

    if not output_all_paths:
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


def output_amplicon_solution(
    bp_graph: BreakpointGraph, output_file: io.TextIOWrapper
) -> None:
    output_file.write("List of cycle segments:\n")

    walk_indices = sorted(
        [(0, i) for i in range(len(bp_graph.walk_weights.cycles))]
        + [(1, i) for i in range(len(bp_graph.walk_weights.paths))],
        key=lambda item: bp_graph.walk_weights[item[0]][item[1]],
        reverse=True,
    )
    heaviest_walk = walk_indices[0]
    if heaviest_walk[0] == 1:
        output_file.write("Largest graph walk solved was a path.\n")
        path_output = get_single_path_output(
            bp_graph, heaviest_walk[1], 1, output_all_paths=False
        )
        output_file.write(path_output)
        return

    output_file.write("Largest graph walk solved was a cycle.\n")

    for cycle in bp_graph.walks.cycles:
        output_file.write(f"{cycle}\n")


def output_amplicon_info(
    bp_graph: BreakpointGraph, output_file: io.TextIOWrapper, was_solved: bool
) -> None:
    output_file.write(f"AmpliconID = {bp_graph.amplicon_idx+1}\n")
    output_file.write(f"#Intervals = {len(bp_graph.amplicon_intervals)}\n")
    output_file.write("AmpliconIntervals:\n")
    for interval in bp_graph.amplicon_intervals:
        output_file.write(f"{interval}\n")
    output_file.write(f"Total Amplicon Size: {bp_graph.total_interval_size}\n")

    output_file.write(f"# Chromosomes: {bp_graph.num_chromosomes}\n")
    output_file.write(f"# Sequence Edges: {bp_graph.num_seq_edges}\n")
    output_file.write(f"# Concordant Edges: {bp_graph.num_conc_edges}\n")
    output_file.write(f"# Discordant Edges: {bp_graph.num_disc_edges}\n")
    output_file.write(f"# Non-Source Edges: {bp_graph.num_nonsrc_edges}\n")
    output_file.write(f"# Source Edges: {bp_graph.num_src_edges}\n")

    if was_solved:
        output_amplicon_solution(bp_graph, output_file)
    else:
        output_file.write("Amplicon was unsolved.\n")
