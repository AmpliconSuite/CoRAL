"""Global state variables"""

from __future__ import annotations

# Updated upon initialization during cli.py and referenced by
# `core_utils.profile_fn` decorator
import pathlib
import time
from dataclasses import dataclass, field

from coral import datatypes
from coral.datatypes import FnCall

# Time required to parse solution and generate summary file when we are unable
# to run to completion (find optimal/best feasible solution for all provided
# breakpoint graphs).
CLEAN_EXIT_BUFFER_S = 10


@dataclass
class GlobalStateProvider:
    should_profile: bool = False
    output_prefix: str = str(pathlib.Path.cwd())
    #output_dir: pathlib.Path = field(default_factory=pathlib.Path.cwd)
    time_limit_s: int = 21600  # 6 hrs in seconds
    start_time: float = field(default_factory=time.time)

    @property
    def summary_filepath(self) -> pathlib.Path:
        if self.output_prefix == str(pathlib.Path.cwd()):
            return pathlib.Path(self.output_prefix) / "amplicon_summary.txt"
        return pathlib.Path(self.output_prefix + "_amplicon_summary.txt")

    @property
    def remaining_time_s(self) -> float:
        return max(
            (
                self.time_limit_s
                - (time.time() - self.start_time)
                - CLEAN_EXIT_BUFFER_S
            ),
            0,
        )


STATE_PROVIDER = GlobalStateProvider()


PROFILED_FN_CALLS: dict[FnCall, datatypes.ProfileResult] = {}
