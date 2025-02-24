from __future__ import annotations

import io
import logging
import random
from pathlib import Path
from typing import Dict

from coral import constants, core_utils
from coral.breakpoint import infer_breakpoint_graph
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.constants import CHR_TAG_TO_IDX
from coral.datatypes import FinalizedPathConstraint, Walk
from coral.output.amplicon_summary import (
    get_single_cycle_output,
    get_single_path_output,
    output_amplicon_info,
)

logger = logging.getLogger(__name__)


def output_amplicon_walks(
    bp_graph: BreakpointGraph,
    output_dir: str,
    *,
    output_all_path_constraints: bool = True,
) -> None:
    """Write the result from cycle decomposition into *.cycles files"""

    amplicon_idx = bp_graph.amplicon_idx
    logger.info(f"Output cycles for amplicon {amplicon_idx+1}.")
    fp = open(f"{output_dir}/amplicon{amplicon_idx + 1}_cycles.txt", "w")
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

    if output_all_path_constraints:
        fp.write("List of all subpath constraints\n")
        for path_idx, basic_pc in enumerate(bp_graph.path_constraints):
            write_path_constraint_to_file(
                FinalizedPathConstraint(
                    edge_counts={},
                    pc_idx=path_idx,
                    support=0,
                ),
                basic_pc.path,
                fp,
            )
    else:
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
            write_path_constraint_to_file(finalized_pc, basic_pc.path, fp)
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
            output_str = get_single_cycle_output(
                bp_graph, walk_idx, walk_indices.index((0, walk_idx)) + 1
            )
            fp.write(output_str)
        else:  # paths
            output_str = get_single_path_output(
                bp_graph, walk_idx, walk_indices.index((1, walk_idx)) + 1
            )
            fp.write(output_str)
    fp.close()


def output_summary_amplicon_stats(
    was_amplicon_solved: Dict[int, bool],
    bp_graphs: list[BreakpointGraph],
    output_dir: str,
) -> None:
    logger.info("Outputting solution info for all amplicons.")

    with Path(f"{output_dir}/amplicon_summary.txt").open("w") as fp:
        fp.write(
            f"{sum(was_amplicon_solved.values())}/{len(bp_graphs)} amplicons "
            "solved.\n"
        )
        for amplicon_idx, bp_graph in enumerate(bp_graphs):
            fp.write(
                "------------------------------------------------------------\n"
            )
            output_amplicon_info(
                bp_graph, fp, was_amplicon_solved[amplicon_idx]
            )
    fp.close()


def write_path_constraint_to_file(
    finalized_pc: FinalizedPathConstraint, path: Walk, fp: io.TextIOWrapper
) -> None:
    fp.write(f"Path constraint\t{finalized_pc.pc_idx+1}\t")
    fp.write(core_utils.path_to_str(path, finalized_pc.edge_counts))
    fp.write(f"\tSupport={finalized_pc.support}")
