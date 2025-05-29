from __future__ import annotations

import io
import logging
import pathlib

from coral import core_utils
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.constants import CHR_TAG_TO_IDX
from coral.datatypes import FinalizedPathConstraint, OutputPCOptions, Walk
from coral.output.utils import (
    get_single_cycle_str,
    get_single_path_str,
)

logger = logging.getLogger(__name__)


def output_amplicon_walks(
    bp_graph: BreakpointGraph,
    output_prefix: str,
    *,
    pc_output_option: OutputPCOptions = OutputPCOptions.LONGEST,
) -> None:
    """Write the result from cycle decomposition into *.cycles files"""

    amplicon_idx = bp_graph.amplicon_idx
    logger.info(f"Output cycles for amplicon {amplicon_idx+1}.")
    cycle_path = pathlib.Path(output_prefix + f"_amplicon{amplicon_idx + 1}_cycles.txt")
    fp = cycle_path.open("w")

    interval_num = 1
    ai_amplicon = sorted(
        bp_graph.amplicon_intervals,
        key=lambda ai: (CHR_TAG_TO_IDX[ai.chr], ai.start),
    )
    for ai in ai_amplicon:
        fp.write(f"Interval\t{interval_num}\t{ai.chr}\t{ai.start}\t{ai.end}\n")
        interval_num += 1

    fp.write("List of cycle segments\n")
    for seqi in range(len(bp_graph.sequence_edges)):
        sseg = bp_graph.sequence_edges[seqi]
        fp.write(f"Segment\t{seqi + 1}\t{sseg.chr}\t{sseg.start}\t{sseg.end}\n")

    satisfied_path_constraints = bp_graph.path_constraints_satisfied
    longest_path_constraints = bp_graph.longest_path_constraints

    if pc_output_option == OutputPCOptions.ALL and len(bp_graph.path_constraints) > 0:
        fp.write("List of all subpath constraints\n")
        for path_idx, basic_pc in enumerate(bp_graph.path_constraints):
            edge_counts_idx = core_utils.path_to_edge_count(basic_pc.path)
            write_path_constraint_to_file(
                path_idx,
                FinalizedPathConstraint(
                    edge_counts=edge_counts_idx,
                    pc_idx=path_idx,
                    support=basic_pc.support,
                ),
                basic_pc.path,
                fp,
            )
            fp.write("\n")
    elif pc_output_option == OutputPCOptions.LONGEST and len(bp_graph.longest_path_constraints) > 0:
        fp.write("List of longest subpath constraints\n")
        satisfied_pc_idxs = {
            path_idx
            for path in (
                satisfied_path_constraints.cycles
                + satisfied_path_constraints.paths
            )
            for path_idx in path
        }
        for constraint_i, finalized_pc in enumerate(longest_path_constraints):
            basic_pc = bp_graph.path_constraints[finalized_pc.pc_idx]
            basic_pc.support = finalized_pc.support
            write_path_constraint_to_file(constraint_i, finalized_pc, basic_pc.path, fp)
            if constraint_i in satisfied_pc_idxs:
                fp.write("\tSatisfied\n")
            else:
                fp.write("\tUnsatisfied\n")

    walk_weights = bp_graph.walk_weights

    # sort walks according to weights
    walk_indices = sorted(
        [(0, i) for i in range(len(walk_weights.cycles))]
        + [(1, i) for i in range(len(walk_weights.paths))],
        key=lambda item: walk_weights[item[0]][item[1]],
        reverse=True,
    )

    fp.write("List of extracted cycles/paths\n")
    for walk_type, walk_idx in walk_indices:
        if walk_type == 0:  # cycles
            if pc_output_option == OutputPCOptions.LONGEST and len(bp_graph.longest_path_constraints) > 0:
                output_str = get_single_cycle_str(
                    bp_graph, walk_idx, walk_indices.index((0, walk_idx)) + 1, 
                    output_path_constraints = True
                )
            else:
                output_str = get_single_cycle_str(
                    bp_graph, walk_idx, walk_indices.index((0, walk_idx)) + 1, 
                    output_path_constraints = False
                )
            fp.write(output_str)
        else:  # paths
            if pc_output_option == OutputPCOptions.LONGEST and len(bp_graph.longest_path_constraints) > 0:
                output_str = get_single_path_str(
                    bp_graph, walk_idx, walk_indices.index((1, walk_idx)) + 1,
                    output_path_constraints = True
                )
            else:
                output_str = get_single_path_str(
                    bp_graph, walk_idx, walk_indices.index((1, walk_idx)) + 1,
                    output_path_constraints = False
                )
            fp.write(output_str)
    fp.close()


def write_path_constraint_to_file(
    idx: int, finalized_pc: FinalizedPathConstraint, path: Walk, fp: io.TextIOWrapper
) -> None:
    fp.write(f"Path constraint {idx + 1}\t")
    fp.write(core_utils.path_to_str(path, finalized_pc.edge_counts))
    fp.write(f"\tSupport={finalized_pc.support}")
