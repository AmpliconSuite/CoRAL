"""Tests for the run summary file (``*_summary.txt``).

The summary doubles as a marker that the reconstruct stage ran (AmpliconClassifier
keys off its presence), so it must be written even when no amplicons were found.

The optional ``Solver Settings`` and ``Resource Usage Summary`` sections are
appended only under ``--profile``, and only once a run completes; the parser must
therefore treat both as optional rather than assuming a fixed section layout.
"""

from __future__ import annotations

import pathlib
from collections.abc import Iterator

import pytest

from coral import datatypes, global_state, text_utils
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.summary import output as summary_output
from coral.summary import parsing as summary_parsing


@pytest.fixture
def summary_prefix(tmp_path: pathlib.Path) -> Iterator[pathlib.Path]:
    """Point global state at a temp prefix, restoring it afterwards."""
    provider = global_state.STATE_PROVIDER
    saved = (provider.output_prefix, provider.should_profile)
    saved_profiles = dict(global_state.PROFILED_FN_CALLS)
    provider.output_prefix = str(tmp_path / "sample")
    provider.should_profile = False
    yield tmp_path / "sample"
    provider.output_prefix, provider.should_profile = saved
    global_state.PROFILED_FN_CALLS.clear()
    global_state.PROFILED_FN_CALLS.update(saved_profiles)


def make_graphs(num_amplicons: int) -> list[BreakpointGraph]:
    return [BreakpointGraph(amplicon_idx=i) for i in range(num_amplicons)]


def write_summary(bp_graphs: list[BreakpointGraph], solved: bool | None = False) -> pathlib.Path:
    summary_output.output_summary_amplicon_stats(
        {bp_graph.amplicon_idx: solved for bp_graph in bp_graphs}, bp_graphs
    )
    return global_state.STATE_PROVIDER.summary_filepath


def test_unperformed_amplicons(
    summary_prefix: pathlib.Path,
) -> None:
    text = write_summary(make_graphs(2), solved = None).read_text()

    assert text.count("Cycle Decomposition Status: UNPERFORMED") == 2
    assert "FAILURE" not in text
    assert "0/2 amplicons solved." in text


def add_resource_usage(output_dir: pathlib.Path) -> None:
    global_state.STATE_PROVIDER.should_profile = True
    global_state.PROFILED_FN_CALLS[
        datatypes.FnCall("reconstruct_graphs", None)
    ] = datatypes.ProfileResult(peak_ram_gb=1.5, runtime_s=42.0)
    summary_output.add_resource_usage_summary(
        datatypes.SolverOptions(
            num_threads=4,
            time_limit_s=7200,
            output_dir=output_dir,
            output_prefix="sample",
            model_prefix="pyomo",
            solver=datatypes.Solver.GUROBI,
        )
    )


# ---------------------------------------------------------------------------
# The file must exist even with no amplicons (e.g. an empty seeds file).
# ---------------------------------------------------------------------------


def test_summary_written_with_no_amplicons(summary_prefix: pathlib.Path) -> None:
    summary_path = write_summary([])

    assert summary_path.exists()
    assert "0/0 amplicons solved." in summary_path.read_text()


def test_summary_header_counts_all_amplicons(
    summary_prefix: pathlib.Path,
) -> None:
    """Unsolved amplicons still count towards the denominator."""
    summary_path = write_summary(make_graphs(3))

    assert "0/3 amplicons solved." in summary_path.read_text()


# ---------------------------------------------------------------------------
# Optional trailing sections: present only under --profile, and only once the
# run finishes.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("num_amplicons", [0, 1, 2])
def test_parse_unprofiled_summary(
    summary_prefix: pathlib.Path, num_amplicons: int
) -> None:
    summary_path = write_summary(make_graphs(num_amplicons))

    parsed = summary_parsing.parse_full_summary(summary_path)

    assert len(parsed.amplicon_summaries) == num_amplicons
    assert parsed.solver_used is None
    assert parsed.profiles_by_fn == {}


@pytest.mark.parametrize("num_amplicons", [0, 1, 2])
def test_parse_profiled_summary(
    summary_prefix: pathlib.Path, tmp_path: pathlib.Path, num_amplicons: int
) -> None:
    summary_path = write_summary(make_graphs(num_amplicons))
    add_resource_usage(tmp_path)

    parsed = summary_parsing.parse_full_summary(summary_path)

    assert len(parsed.amplicon_summaries) == num_amplicons
    assert parsed.solver_used == datatypes.Solver.GUROBI
    assert parsed.threads_used == 4
    assert parsed.solver_time_limit == 7200
    assert "reconstruct_graphs" in parsed.profiles_by_fn


def test_parse_interrupted_profiled_summary(
    summary_prefix: pathlib.Path,
) -> None:
    """Profiling on in the header, but the run died before the trailing sections.

    The amplicon sections must survive rather than being consumed as resource
    usage.
    """
    global_state.STATE_PROVIDER.should_profile = True
    summary_path = write_summary(make_graphs(2))
    assert text_utils.RESOURCE_USAGE_HEADER not in summary_path.read_text()

    parsed = summary_parsing.parse_full_summary(summary_path)

    assert parsed.profiling_enabled is True
    assert len(parsed.amplicon_summaries) == 2
    assert parsed.solver_used is None
