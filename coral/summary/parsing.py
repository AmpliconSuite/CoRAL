from __future__ import annotations

import dataclasses
import functools
import io
import pathlib
import re
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pandera as pa
import pandera.typing as pat

from coral import cycle2bed, datatypes, text_utils
from coral.text_utils import COMMA_NUMBER_REGEX

THREAD_PATTERN = re.compile(r"Threads: (\d+)")
RESOURCE_PATTERN_MULTICALL_FN = re.compile(
    r"(\S+)/(\d+) \| Peak RAM: ([\d.]+) GB \| Runtime: ([\d.]+) s"
)
RESOURCE_PATTERN_SINGLECALL_FN = re.compile(
    r"(\S+)(?<!/\d) \| Peak RAM: ([\d.]+) GB \| Runtime: ([\d.]+) s"
)


@dataclass
class AmpliconSummary:
    total_amplicon_size: int | None = None
    num_breakpoints: int | None = None
    profiles_by_fn: dict[str, datatypes.ProfileResult] = field(
        default_factory=dict
    )
    are_cycles_solved: bool = False
    was_suboptimal_solution: bool = False
    model_metadata: datatypes.ModelMetadata | None = None


@dataclass
class FullProfileSummary:
    version: str | None = None
    profiling_enabled: bool | None = None
    threads_used: int | None = None
    solver_used: datatypes.Solver | None = None
    solver_time_limit: int | None = None

    profiles_by_fn: dict[str, datatypes.ProfileResult] = field(
        default_factory=dict
    )
    amplicon_summaries: dict[int, AmpliconSummary] = field(
        default_factory=lambda: defaultdict(AmpliconSummary)
    )

    @functools.cached_property
    def max_overall_ram_usage(self) -> float:
        return max(
            *(
                profile_result.peak_ram_gb
                for profile_result in self.profiles_by_fn.values()
                if profile_result.peak_ram_gb is not None
            )
        )


class FullSummaryModel(pa.DataFrameModel):
    # Match Full Profile Summary

    version: pat.Series[str]
    profiling_enabled: pat.Series[bool]
    threads_used: pat.Series[int]
    solver_used: pat.Series[str]
    solver_time_limit: pat.Series[int]
    profiles_by_fn: pat.Series[dict[str, datatypes.ProfileResult]]
    amplicon_summaries: pat.Series[dict[int, AmpliconSummary]]

    n_graphs: pat.Series[int]
    n_solved: pat.Series[int]
    runtime_s: pat.Series[float]
    peak_ram_gb: pat.Series[float]


class AmpliconSummaryModel(pa.DataFrameModel):
    # Match Amplicon Summary
    total_amplicon_size: pat.Series[int]
    num_breakpoints: pat.Series[int]
    profiles_by_fn: pat.Series[dict[str, datatypes.ProfileResult]]
    are_cycles_solved: pat.Series[bool]
    was_suboptimal_solution: pat.Series[bool]
    model_metadata: pat.Series[datatypes.ModelMetadata] = pa.Field(
        nullable=True
    )


def get_fn_call_resource_info(
    pattern: re.Pattern[str], line: str
) -> tuple[datatypes.FnCall, datatypes.ProfileResult] | None:
    if not (matched_groups := pattern.match(line)):
        return None
    if pattern == RESOURCE_PATTERN_MULTICALL_FN:
        fn_call = datatypes.FnCall(
            fn_name=matched_groups.group(1),
            call_ctr=int(matched_groups.group(2)),
        )
        ram_gb, runtime_s = matched_groups.groups()[2:]
    else:
        fn_call = datatypes.FnCall(
            fn_name=matched_groups.group(1),
            call_ctr=None,
        )
        ram_gb, runtime_s = matched_groups.groups()[1:]
    profile_result = datatypes.ProfileResult(
        peak_ram_gb=float(ram_gb),
        runtime_s=float(runtime_s),
    )
    return fn_call, profile_result


def parse_header(header_chunk: str) -> FullProfileSummary:
    if not (version := text_utils.VERSION_PATTERN.search(header_chunk)):
        raise ValueError("Could not parse version from header")
    if not (enabled := text_utils.PROFILE_ENABLED_PATTERN.search(header_chunk)):
        raise ValueError("Could not parse enabled from header")
    return FullProfileSummary(
        version=version.group(1),
        profiling_enabled=enabled.group(1) == "True",
    )


def parse_solver_summary(
    summary_str: str, full_profile_summary: FullProfileSummary
) -> None:
    """Parse solver info from summary string.

    ex: '\nSolver Settings: \nSolver: GUROBI\nThreads: 2\nTime Limit: 28800 s\n'
    """
    if not (match := text_utils.SOLVER_PATTERN.search(summary_str)):
        # print(summary_str)
        raise ValueError(f"Could not parse solver from {summary_str}")
    full_profile_summary.solver_used = datatypes.Solver[match.group(1)]

    if not (match := text_utils.THREADS_PATTERN.search(summary_str)):
        raise ValueError(f"Could not parse threads from {summary_str}")
    full_profile_summary.threads_used = int(match.group(1))

    if not (match := text_utils.TIME_LIMIT_PATTERN.search(summary_str)):
        raise ValueError(f"Could not parse time limit from {summary_str}")
    full_profile_summary.solver_time_limit = int(match.group(1))


def parse_resource_summary(
    summary_str: str, full_profile_summary: FullProfileSummary
) -> None:
    """Parse resource info from summary string.

    ex: '
        Resource Usage Summary:
        reconstruct_graphs | Peak RAM: 1.370868797 GB | Runtime: 1274.0610159100033 s
        solve_single_graph/1 | Peak RAM: 0.016169349 GB | Runtime: 28830.582218587006 s\n
        solve_single_graph/2 | Peak RAM: 0.004381079 GB | Runtime: 2.2074499869922874 s\n'
    """

    for line in summary_str.splitlines():
        # Single-call FN, i.e., pertaining to shared graph generation logic
        if fn_profile := get_fn_call_resource_info(
            RESOURCE_PATTERN_SINGLECALL_FN, line
        ):
            fn_call, profile_result = fn_profile
            full_profile_summary.profiles_by_fn[fn_call.fn_name] = (
                profile_result
            )
        # Multi-call FN, i.e., pertaining to amplicon-specific solving logic
        if fn_profile := get_fn_call_resource_info(
            RESOURCE_PATTERN_MULTICALL_FN, line
        ):
            fn_call, profile_result = fn_profile
            full_profile_summary.amplicon_summaries[
                fn_call.call_ctr  # type: ignore[index]
            ].profiles_by_fn[fn_call.fn_name] = profile_result
            continue


def parse_amplicon_summary(summary_str: str) -> AmpliconSummary:
    amplicon_summary = AmpliconSummary()
    for line in summary_str.splitlines():
        if match := text_utils.AMPLICON_SIZE_PATTERN.search(line):
            amplicon_size = int(match.group(1).replace(",", ""))
            amplicon_summary.total_amplicon_size = amplicon_size
            continue
        if match := text_utils.AMPLICON_BREAKPOINTS_PATTERN.search(line):
            amplicon_summary.num_breakpoints = int(
                match.group(1).replace(",", "")
            )
            continue
        if match := text_utils.CYCLE_DECOMP_STATUS_PATTERN.search(line):
            amplicon_summary.are_cycles_solved = match.group(1) == "SUCCESS"
            continue
        if match := text_utils.MODEL_METADATA_PATTERN.search(line):
            alpha, weights, resolution, num_path_constraints = match.groups()[
                2:
            ]
            amplicon_summary.model_metadata = datatypes.ModelMetadata(
                model_type=datatypes.ModelType[match.group(1)],
                k=int(match.group(2)),
                alpha=float(alpha) if alpha != "None" else None,
                total_weights=float(weights) if weights != "None" else None,
                resolution=float(resolution) if resolution != "None" else None,
                num_path_constraints=int(num_path_constraints),
            )
            continue
        if match := text_utils.MODEL_METADATA_PATTERN_V2_1.search(line):
            alpha, weights, resolution = match.groups()[2:]
            amplicon_summary.model_metadata = datatypes.ModelMetadata(
                model_type=datatypes.ModelType[match.group(1)],
                k=int(match.group(2)),
                alpha=float(alpha) if alpha != "None" else None,
                total_weights=float(weights) if weights != "None" else None,
                resolution=float(resolution) if resolution != "None" else None,
            )
            continue
        if text_utils.SUBOPTIMAL_WARNING in line:
            if match := text_utils.MIP_GAP_PATTERN.search(line):
                prev_metadata_dict = (
                    amplicon_summary.model_metadata._asdict()  # type: ignore[union-attr]
                )
                prev_metadata_dict["mip_gap"] = float(match.group(1))
                amplicon_summary.model_metadata = datatypes.ModelMetadata(
                    **prev_metadata_dict
                )
            amplicon_summary.was_suboptimal_solution = True
            continue
    return amplicon_summary


def parse_full_summary(path: pathlib.Path) -> FullProfileSummary:
    summary_chunks = path.read_text().split(text_utils.AMPLICON_SEPARATOR)
    full_profile_summary = parse_header(summary_chunks.pop(0))

    resource_summary = None
    if full_profile_summary.profiling_enabled:
        resource_summary = summary_chunks.pop()

    # Update solver info in-place
    parse_solver_summary(summary_chunks.pop(), full_profile_summary)

    amplicon_summaries = summary_chunks
    for i, amplicon_summary_str in enumerate(amplicon_summaries):
        amplicon_summary = parse_amplicon_summary(amplicon_summary_str)
        full_profile_summary.amplicon_summaries[i + 1] = amplicon_summary

    # Defer resource parsing until after amplicon summaries are parsed, since
    # we track per-amplicon usages.
    if resource_summary:
        parse_resource_summary(resource_summary, full_profile_summary)

    return full_profile_summary


def get_full_summaries_from_directory(
    directory: pathlib.Path,
) -> dict[str, FullProfileSummary]:
    summary_stats_by_dataset = {}
    for dataset_path in directory.iterdir():
        summary_path = dataset_path / "amplicon_summary.txt"
        if summary_path.exists():
            try:
                summary_stats = parse_full_summary(summary_path)
                summary_stats_by_dataset[dataset_path.name] = summary_stats
            except Exception as e:
                print(f"Unable to parse summary at {summary_path}")
    return summary_stats_by_dataset


def get_summary_stats_df(
    directory: pathlib.Path,
) -> pat.DataFrame[FullSummaryModel]:
    """Adds n_graphs, n_solved, runtime and peak RAM columns to a summary
     stats DataFrame.

    Args:
        stats_df: DataFrame containing amplicon summary statistics

    Returns:
        DataFrame with generated summary statistics columns
    """
    summary_stats_by_dataset = get_full_summaries_from_directory(directory)
    raw_df = pd.DataFrame(
        [
            dataclasses.asdict(stats)
            for stats in summary_stats_by_dataset.values()
        ],
        index=summary_stats_by_dataset.keys(),
    )
    stats_df = FullSummaryModel.validate(raw_df)

    stats_df["n_graphs"] = stats_df.amplicon_summaries.apply(lambda x: len(x))
    stats_df["n_solved"] = stats_df.amplicon_summaries.apply(
        lambda x: sum(i.are_cycles_solved for i in x.values())
    )

    return stats_df


def get_amplicon_summary_df(
    full_summary_by_name: dict[str, FullProfileSummary],
) -> pat.DataFrame[AmpliconSummaryModel]:
    default_per_amplicon = pd.DataFrame(
        [
            {**dataclasses.asdict(summary), "dataset": dataset, "amplicon": idx}
            for dataset, full_stats in full_summary_by_name.items()
            for idx, summary in full_stats.amplicon_summaries.items()
        ]
    )
    default_per_amplicon.index = default_per_amplicon.apply(
        lambda x: (x["dataset"], x["amplicon"]), axis=1
    )
    default_per_amplicon["runtime_h"] = default_per_amplicon[
        "profiles_by_fn"
    ].apply(
        lambda x: x["solve_single_graph"].runtime_s / 3600
        if "solve_single_graph" in x
        else None
    )
    default_per_amplicon["peak_ram_gb"] = default_per_amplicon[
        "profiles_by_fn"
    ].apply(
        lambda x: x["solve_single_graph"].peak_ram_gb
        if "solve_single_graph" in x
        else None
    )
    default_per_amplicon["model_used"] = default_per_amplicon[
        "model_metadata"
    ].apply(
        lambda x: x.model_type
        if x is not None
        else datatypes.ModelType.GREEDY  # Remove this default once latest batch of summaries produced
    )
    return AmpliconSummaryModel.validate(default_per_amplicon)


# Add per-amplicon runtime and RAM usage
def get_solve_graph_stats(x: pd.DataFrame, stat_name: str) -> float | None:
    return (
        x.profiles_by_fn["solve_single_graph"][stat_name]
        if "solve_single_graph" in x.profiles_by_fn
        else None
    )


def parse_cycle_file(
    cycle_file: io.TextIOWrapper,
    output_prefix: str,
    num_cycles: int | None,
    print_command: bool = False,
) -> dict[int, datatypes.ReconstructedCycle]:
    """
    Parse a '*_cycles.{txt,bed}' file  into a dict of
    cycle ID to list of cycle segments:

    Example .bed file:
        #chr	start	end	orientation	cycle_id	iscyclic	weight
        chr6	39545281	41545280	+	1	True	6.338641016903567
    """

    cycle_path = pathlib.Path(cycle_file.name).resolve()
    # check if it ends with .bed, if not convert it
    bed_path = None
    if cycle_path.suffix == ".txt":
        # convert it to a bed
        num_str = str(num_cycles) if num_cycles else "all"
        prefix = (
            output_prefix
            if output_prefix.endswith("/")
            else f"{output_prefix}_"
        )
        bed_filename = f"{prefix}converted_{num_str}_cycles.bed"
        cycle2bed.convert_cycles_to_bed(
            cycle_file, bed_filename, num_cycles=num_cycles, print_command = print_command
        )
        bed_path = pathlib.Path(bed_filename)
    elif cycle_path.suffix == ".bed":
        bed_path = cycle_path
    else:
        raise ValueError(
            f"Unsupported cycle file extension: {cycle_path.suffix}"
        )

    cycles: dict[int, datatypes.ReconstructedCycle] = {}
    with bed_path.open("r") as bed_file:
        for line in bed_file:
            chunks = line.strip().split("\t")
            if chunks[0][0] == "#":
                continue
            cycle_id = int(chunks[4])
            intv = datatypes.CycleInterval(
                chunks[0],
                int(chunks[1]),
                int(chunks[2]),
                datatypes.Strand(chunks[3]),
                cn=float(chunks[6]),
            )
            if cycle_id not in cycles:
                cycles[cycle_id] = datatypes.ReconstructedCycle(
                    cycle_id=cycle_id,
                    is_cyclic=chunks[5] == "True",
                )
            cycles[cycle_id].segments.append(intv)

    return cycles
