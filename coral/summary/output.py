from __future__ import annotations

import dataclasses
import importlib.metadata
import io
import logging
import os

from coral import datatypes, global_state, text_utils
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.output.utils import get_single_cycle_str, get_single_path_str

logger = logging.getLogger(__name__)


def output_amplicon_solution(
    bp_graph: BreakpointGraph, output_file: io.TextIOWrapper
) -> None:
    walk_indices = sorted(
        [(0, i) for i in range(len(bp_graph.walk_weights.cycles))]
        + [(1, i) for i in range(len(bp_graph.walk_weights.paths))],
        key=lambda item: bp_graph.walk_weights[item[0]][item[1]],
        reverse=True,
    )
    heaviest_walk = walk_indices[0]
    if heaviest_walk[0] == 1:
        output_file.write("\tHeaviest graph walk solved was a path.\n")
        path_output = get_single_path_str(
            bp_graph, heaviest_walk[1], 1, output_path_constraints=False
        )
        output_file.write(f"\t{path_output}")
    else:
        output_file.write("\tHeaviest graph walk solved was a cycle.\n")
        cycle_output = get_single_cycle_str(
            bp_graph, heaviest_walk[1], 1, output_path_constraints=False
        )
        output_file.write(f"\t{cycle_output}")


def output_amplicon_info(
    bp_graph: BreakpointGraph, output_file: io.TextIOWrapper, was_solved: bool
) -> None:
    output_file.write(text_utils.AMPLICON_SEPARATOR + "\n")
    output_file.write(f"AmpliconID = {bp_graph.amplicon_idx+1}\n")
    output_file.write(f"#Intervals = {len(bp_graph.amplicon_intervals)}\n")
    output_file.write("Amplicon Intervals:\n")
    for interval in bp_graph.amplicon_intervals:
        output_file.write(f"\t{interval}\n")
    output_file.write(
        f"Total Amplicon Size: {bp_graph.total_interval_size:,d}\n"
    )

    output_file.write(f"# Chromosomes: {bp_graph.num_chromosomes}\n")
    output_file.write(f"# Sequence Edges: {bp_graph.num_seq_edges}\n")
    output_file.write(f"# Concordant Edges: {bp_graph.num_conc_edges}\n")
    output_file.write(f"# Discordant Edges: {bp_graph.num_disc_edges}\n")
    output_file.write(f"# Non-Source Edges: {bp_graph.num_nonsrc_edges}\n")
    output_file.write(f"# Source Edges: {len(bp_graph.amplicon_intervals) * 4}\n")

    if was_solved:
        output_file.write(
            f"{text_utils.CYCLE_DECOMP_STATUS_TEMPLATE.format(status='SUCCESS')}\n"
        )
        if model_metadata := bp_graph.model_metadata:
            output_file.write(f"\t{model_metadata.to_output_str()}\n")
        if (
            bp_graph.relative_mip_gap is not None
            and bp_graph.relative_mip_gap > 1  # 1% gap
        ):
            output_file.write(
                f"\t{text_utils.SUBOPTIMAL_WARNING} (Relative MIP Gap: "
                f"{bp_graph.relative_mip_gap:.5f}%)\n"
            )
        if bp_graph.upper_bound is not None:
            output_file.write(f"\tUpper Bound: {bp_graph.upper_bound:.5f}\n")
        output_amplicon_solution(bp_graph, output_file)
    else:
        output_file.write(
            f"{text_utils.CYCLE_DECOMP_STATUS_TEMPLATE.format(status='FAILURE')}\n"
        )
        if model_metadata := bp_graph.model_metadata:
            output_file.write(f"\t{model_metadata.to_output_str()}\n")


def add_resource_usage_summary(solver_options: datatypes.SolverOptions) -> None:
    with global_state.STATE_PROVIDER.summary_filepath.open("a") as fp:
        fp.write(text_utils.AMPLICON_SEPARATOR + "\n")
        fp.write("Solver Settings: \n")
        fp.write(f"Solver: {solver_options.solver.name}\n")
        threads_used = (
            solver_options.num_threads
            if solver_options.num_threads != -1
            else os.cpu_count()
        )
        fp.write(f"Threads: {threads_used}\n")
        fp.write(f"Time Limit: {solver_options.time_limit_s} s\n")
        fp.write(text_utils.AMPLICON_SEPARATOR + "\n")
        fp.write("Resource Usage Summary:\n")
        for fn_call, profile in sorted(global_state.PROFILED_FN_CALLS.items()):
            fn_tag = (
                f"{fn_call.fn_name}/{fn_call.call_ctr}"
                if fn_call.call_ctr is not None
                else fn_call.fn_name
            )
            fp.write(
                f"{fn_tag} | Peak RAM: {profile.peak_ram_gb:.5f} GB | "
                f"Runtime: {profile.runtime_s:.5f} s\n"
            )


def get_summary_header(was_amplicon_solved: dict[int, bool]) -> str:
    header_str = f"{text_utils.VERSION_TEMPLATE.format(
        version=importlib.metadata.version("coral")
    )}\n"
    header_str += f"{sum(was_amplicon_solved.values())}/{len(was_amplicon_solved)} amplicons solved.\n"
    header_str += (
        f"Runtime Limit: {global_state.STATE_PROVIDER.time_limit_s} s\n"
    )
    header_str += f"{text_utils.PROFILE_ENABLED_TEMPLATE.format(
        enabled=global_state.STATE_PROVIDER.should_profile)}\n"

    return header_str


def output_summary_amplicon_stats(
    was_amplicon_solved: dict[int, bool],
    bp_graphs: list[BreakpointGraph],
) -> None:
    logger.info("Outputting solution info for all amplicons.")

    with global_state.STATE_PROVIDER.summary_filepath.open("w") as fp:
        fp.write(get_summary_header(was_amplicon_solved))
        for bp_graph in bp_graphs:
            output_amplicon_info(
                bp_graph, fp, was_amplicon_solved[bp_graph.amplicon_idx]
            )
