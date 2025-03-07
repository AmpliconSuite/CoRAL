from __future__ import annotations

import functools
import pathlib
import re
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
import pandera as pa

from coral import datatypes
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
    profiles_by_fn: dict[str, datatypes.ProfileResult] = field(
        default_factory=dict
    )


@dataclass
class FullProfileSummary:
    threads_used: int | None = None
    profiles_by_fn: dict[str, datatypes.ProfileResult] = field(
        default_factory=dict
    )
    amplicon_summary_info: dict[int, AmpliconSummary] = field(
        default_factory=dict
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


def get_fn_call_resource_info(
    line: str,
) -> tuple[datatypes.FnCall, datatypes.ProfileResult] | None:
    if matched_groups := RESOURCE_PATTERN_MULTICALL_FN.match(line):
        fn_call = datatypes.FnCall(
            fn_name=matched_groups.group(1),
            call_ctr=int(matched_groups.group(2)),
        )
        profile_result = datatypes.ProfileResult(
            peak_ram_gb=float(matched_groups.group(3)),
            runtime_s=float(matched_groups.group(4)),
        )
        return fn_call, profile_result
    if matched_groups := RESOURCE_PATTERN_SINGLECALL_FN.match(line):
        fn_call = datatypes.FnCall(
            fn_name=matched_groups.group(1),
            call_ctr=None,
        )
        profile_result = datatypes.ProfileResult(
            peak_ram_gb=float(matched_groups.group(2)),
            runtime_s=float(matched_groups.group(3)),
        )
        return fn_call, profile_result
    return None


def parse_amplicon_summary(path: pathlib.Path) -> tuple[int, AmpliconSummary]:
    full_text = path.read_text().split()
    with path.open("r") as f:
        amplicon_chunks = []
        for line in f:
            if total_amplicon_size := COMMA_NUMBER_REGEX.search(line):
                amplicon_summary = AmpliconSummary()
                amplicon_summary.total_amplicon_size = int(
                    total_amplicon_size.group()
                )


def get_resource_usage(
    reconstruction_dir: pathlib.Path,
) -> dict[str, FullProfileSummary]:
    full_profiles_by_dataset: dict[str, FullProfileSummary] = {}
    for dataset_path in reconstruction_dir.iterdir():
        if not (
            summary_path := next(iter(dataset_path.glob("*_summary.txt")), None)
        ):
            continue

        full_profile_summary = FullProfileSummary()
        with summary_path.open("r") as f:
            for line in f:
                if total_amplicon_size := COMMA_NUMBER_REGEX.search(line):
                    amplicon_summary = (
                        full_profile_summary.amplicon_summary_info[
                            int(total_amplicon_size.group())
                        ]
                    )
                    amplicon_summary.total_amplicon_size = int(
                        total_amplicon_size.group()
                    )
                    continue
                if threads_used := THREAD_PATTERN.search(line):
                    full_profile_summary.threads_used = int(
                        threads_used.group(1)
                    )
                    continue
                if not (fn_profile := get_fn_call_resource_info(line)):
                    continue
                fn_call, profile_result = fn_profile
                if fn_call.call_ctr is None:
                    full_profile_summary.profiles_by_fn[fn_call.fn_name] = (
                        profile_result
                    )
                else:
                    full_profile_summary.amplicon_summary_info[
                        fn_call.call_ctr
                    ].profiles_by_fn[fn_call.fn_name] = profile_result
        full_profiles_by_dataset[dataset_path.name] = full_profile_summary
    return full_profiles_by_dataset


def plot_resource_usage(
    reconstruction_dir: pathlib.Path, output_dir: pathlib.Path
) -> None:
    full_profiles_by_dataset = get_resource_usage(reconstruction_dir)
    ram_usage = np.array(
        [
            full_profile_summary.max_overall_ram_usage
            for full_profile_summary in full_profiles_by_dataset.values()
        ]
    )
    print(ram_usage)
    for dataset_name, full_profile_summary in full_profiles_by_dataset.items():
        if full_profile_summary.threads_used is None:
            continue

        for (
            fn_name,
            profile_result,
        ) in full_profile_summary.profiles_by_fn.items():
            print(
                f"{dataset_name} {fn_name} {profile_result.peak_ram_gb} {profile_result.runtime_s}"
            )
