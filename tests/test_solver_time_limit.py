"""Tests for solver results returned without a solution after a time limit.

When the global time limit (``--global-time-limit``) is exhausted, later solves
are launched with a time limit of ~0s. SCIP then stops before finding any
solution and omits the "Dual Bound"/"Gap" lines from the end of its log, so
Pyomo never sets ``results.solver.gap``. Reading that attribute unguarded raised
``AttributeError`` and crashed cycle decomposition, discarding every cycle the
greedy decomposition had already found.
"""

from __future__ import annotations

import pathlib
import shutil
from collections.abc import Iterator
from types import SimpleNamespace
from typing import Any

import pyomo.environ as pyo
import pyomo.opt
import pytest

from coral import cycle_decomposition, datatypes, global_state
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.models import cycle_utils

TEST_GRAPH = pathlib.Path(__file__).parent / "data" / "amplicon1_graph.txt"

requires_scip = pytest.mark.skipif(
    shutil.which("scip") is None, reason="SCIP executable not on PATH"
)


def _no_solution_scip_results() -> pyomo.opt.SolverResults:
    """Mirror what Pyomo returns when SCIP is interrupted with 0 solutions."""
    results = pyomo.opt.SolverResults()
    results.solver.status = pyo.SolverStatus.ok
    results.solver.termination_condition = pyo.TerminationCondition.maxTimeLimit
    results.solver.message = "time limit reached"
    return results


def _raise_no_value() -> float:
    raise ValueError("No value for uninitialized VarData object x[0,0]")


def test_parse_time_limit_result_without_gap_is_aborted() -> None:
    results = _no_solution_scip_results()
    with pytest.raises(AttributeError):
        _ = results.solver.gap  # The precondition for the original crash.

    # Pyomo still registers one (empty) solution, and no variable values are
    # loaded, so evaluating the objective raises.
    model = SimpleNamespace(solutions=[object()], Objective=_raise_no_value)

    parsed = cycle_utils.parse_solver_output(
        results,
        model,  # type: ignore[arg-type]
        SimpleNamespace(amplicon_idx=0),  # type: ignore[arg-type]
        k=1,
        pc_list=[],
        total_weights=1.0,
        model_type=datatypes.CycleDecompOptions.MAX_WEIGHT,
    )

    assert parsed.solver_status == pyo.SolverStatus.aborted
    assert parsed.termination_condition == pyo.TerminationCondition.maxTimeLimit
    assert parsed.num_cycles == 0
    assert parsed.num_paths == 0


@pytest.fixture
def time_limit_state() -> Iterator[global_state.GlobalStateProvider]:
    provider = global_state.STATE_PROVIDER
    saved = (provider.time_limit_s, provider.start_time)
    yield provider
    provider.time_limit_s, provider.start_time = saved


def _load_test_graph() -> BreakpointGraph:
    with TEST_GRAPH.open() as graph_file:
        bp_graph = parse_breakpoint_graph(graph_file)
    bp_graph.amplicon_idx = 0
    bp_graph.infer_amplicon_intervals(0)
    return bp_graph


def _greedy(
    bp_graph: BreakpointGraph, tmp_path: pathlib.Path
) -> datatypes.CycleSolution:
    (tmp_path / "models").mkdir(parents=True, exist_ok=True)
    return cycle_decomposition.maximize_weights_greedy(
        bp_graph=bp_graph,
        total_weights=sum(len(s) * s.cn for s in bp_graph.sequence_edges),
        node_order={
            node: i for i, node in enumerate(bp_graph.node_adjacencies)
        },
        solver_options=datatypes.SolverOptions(
            output_dir=tmp_path,
            num_threads=1,
            time_limit_s=60,
            solver=datatypes.Solver.SCIP,
        ),
        cycle_id=0,
    )


@requires_scip
@pytest.mark.parametrize(
    "mode",
    [
        datatypes.CycleDecompOptions.MAX_WEIGHT,
        datatypes.CycleDecompOptions.MIN_CYCLES,
    ],
)
def test_exhausted_global_time_limit_does_not_crash(
    mode: datatypes.CycleDecompOptions,
    time_limit_state: global_state.GlobalStateProvider,
    tmp_path: pathlib.Path,
) -> None:
    # Every solve now gets a time limit of 0s.
    time_limit_state.time_limit_s = 0
    (tmp_path / "models").mkdir()
    bp_graph = _load_test_graph()
    cycle_decomposition.cycle_decomposition_single_graph(
        bp_graph,
        datatypes.SolverOptions(
            output_dir=tmp_path,
            output_prefix="test",
            num_threads=1,
            time_limit_s=60,
            solver=datatypes.Solver.SCIP,
        ),
        should_force_greedy=mode == datatypes.CycleDecompOptions.MAX_WEIGHT,
    )
    assert not bp_graph.walks.cycles
    assert not bp_graph.walks.paths


@requires_scip
def test_greedy_keeps_walks_found_before_time_limit(
    time_limit_state: global_state.GlobalStateProvider,
    tmp_path: pathlib.Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    time_limit_state.time_limit_s = 3600
    full = _greedy(_load_test_graph(), tmp_path / "full")
    assert (
        full.num_cycles + full.num_paths >= 2
    ), "test graph should need more than one greedy iteration"

    # Exhaust the global time limit right after the first solve.
    original_solve = cycle_utils.PyomoSolverWrapper.solve

    def solve_then_exhaust(
        self: cycle_utils.PyomoSolverWrapper, *args: Any, **kwargs: Any
    ) -> pyomo.opt.SolverResults:
        results = original_solve(self, *args, **kwargs)
        time_limit_state.time_limit_s = 0
        return results

    time_limit_state.time_limit_s = 3600
    monkeypatch.setattr(
        cycle_utils.PyomoSolverWrapper, "solve", solve_then_exhaust
    )
    partial = _greedy(_load_test_graph(), tmp_path / "partial")

    assert partial.num_cycles + partial.num_paths == 1
    assert partial.solver_status == pyo.SolverStatus.ok
    assert full.early_stop is None
    assert partial.early_stop is not None
    assert partial.early_stop.reason.startswith("global time limit")
    assert partial.early_stop.num_iterations == 1
    assert 0 < partial.total_weights_included < full.total_weights_included


@requires_scip
@pytest.mark.parametrize("postprocess", [False, True], ids=["greedy", "post"])
def test_partial_summary_end_to_end(
    postprocess: bool,
    time_limit_state: global_state.GlobalStateProvider,
    tmp_path: pathlib.Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the `cycle_all` driver with the clock exhausted after one solve.

    With --postprocess-greedy-sol, postprocessing then gets no time; the greedy
    walk must still be kept.
    """
    monkeypatch.setattr(
        global_state.STATE_PROVIDER, "output_prefix", str(tmp_path / "test")
    )
    original_solve = cycle_utils.PyomoSolverWrapper.solve

    def solve_then_exhaust(
        self: cycle_utils.PyomoSolverWrapper, *args: Any, **kwargs: Any
    ) -> pyomo.opt.SolverResults:
        results = original_solve(self, *args, **kwargs)
        time_limit_state.time_limit_s = 0
        return results

    time_limit_state.time_limit_s = 3600
    monkeypatch.setattr(
        cycle_utils.PyomoSolverWrapper, "solve", solve_then_exhaust
    )
    (tmp_path / "models").mkdir()
    cycle_decomposition.cycle_decomposition_all_graphs(
        [_load_test_graph()],
        datatypes.SolverOptions(
            output_dir=tmp_path,
            output_prefix="test",
            num_threads=1,
            time_limit_s=60,
            solver=datatypes.Solver.SCIP,
        ),
        should_force_greedy=True,
        should_postprocess=postprocess,
    )
    text = (tmp_path / "test_summary.txt").read_text()

    assert "0/1 amplicons solved." in text
    assert "Cycle Decomposition Status: PARTIAL" in text
    assert (
        "Stopped early: global time limit (0 s) reached after 1 greedy "
        "iterations." in text
    )
    assert "(greedy normally stops at 90%)" in text
    # The walk found before the time limit is still written out.
    assert (tmp_path / "test_amplicon1_cycles.txt").exists()
    # Leave a complete CoRAL-style output (summary + graph + cycles) so the
    # directory can be fed to AmpliconClassifier when run with --basetemp.
    shutil.copy(TEST_GRAPH, tmp_path / "test_amplicon1_graph.txt")
