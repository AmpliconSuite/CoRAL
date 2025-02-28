"""Global state variables"""

from __future__ import annotations

# Updated upon initialization during cli.py and referenced by
# `core_utils.profile_fn` decorator
import pathlib
from dataclasses import dataclass, field

from coral import datatypes
from coral.datatypes import FnCall


@dataclass
class GlobalStateProvider:
    should_profile: bool = False
    output_dir: pathlib.Path = field(default_factory=pathlib.Path.cwd)

    @property
    def summary_filepath(self) -> pathlib.Path:
        return self.output_dir / "amplicon_summary.txt"


STATE_PROVIDER = GlobalStateProvider()


PROFILED_FN_CALLS: dict[FnCall, datatypes.ProfileResult] = {}
