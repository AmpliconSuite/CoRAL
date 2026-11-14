"""Tests for the completeness lines and PARTIAL status in the run summary.

When greedy cycle decomposition is stopped by a time limit before reaching its
own targets, the walks found so far are written out, but the amplicon is
reported as PARTIAL (not counted as solved) along with why it stopped and how
much of the amplicon its walks explain.
"""

from __future__ import annotations

import io

from coral import datatypes
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import WalkData
from coral.summary import output as summary_output

EARLY_STOP = datatypes.GreedyEarlyStop(
    reason="global time limit (28800 s)",
    num_iterations=6,
    target_weight_fraction=0.9,
    target_pc_fraction=0.9,
)


def _graph(
    idx: int,
    *,
    explained: float | None,
    early_stop: datatypes.GreedyEarlyStop | None = None,
    num_pc: int = 0,
    satisfied: WalkData[list[int]] | None = None,
) -> BreakpointGraph:
    bp_graph = BreakpointGraph(amplicon_idx=idx)
    bp_graph.explained_weight_fraction = explained
    bp_graph.early_stop = early_stop
    bp_graph.longest_path_constraints = [None] * num_pc  # type: ignore[list-item]
    if satisfied is not None:
        bp_graph.path_constraints_satisfied = satisfied
    return bp_graph


def _completeness(bp_graph: BreakpointGraph) -> list[str]:
    buf = io.StringIO()
    summary_output.output_completeness(bp_graph, buf)  # type: ignore[arg-type]
    return buf.getvalue().splitlines()


def test_complete_decomposition_reports_explained_fraction_only() -> None:
    assert _completeness(_graph(0, explained=0.9421, num_pc=4)) == [
        "\tCN-weighted length explained: 94.2%"
    ]


def test_partial_decomposition_reports_stop_reason_and_targets() -> None:
    bp_graph = _graph(
        0,
        explained=0.413,
        early_stop=EARLY_STOP,
        num_pc=40,
        satisfied=WalkData([[0, 1], [2]], [[2, 3]]),
    )
    assert _completeness(bp_graph) == [
        "\tStopped early: global time limit (28800 s) reached after 6 greedy "
        "iterations.",
        "\tCN-weighted length explained: 41.3% (greedy normally stops at 90%)",
        "\tPath constraints satisfied: 4/40 (greedy normally stops at 90%)",
    ]


def test_partial_without_path_constraints_omits_that_line() -> None:
    lines = _completeness(_graph(0, explained=0.5, early_stop=EARLY_STOP))
    assert not any("Path constraints" in line for line in lines)


def test_header_does_not_count_partial_amplicons_as_solved() -> None:
    graphs = [
        _graph(0, explained=0.95),
        _graph(1, explained=0.41, early_stop=EARLY_STOP),
        _graph(2, explained=None),
    ]
    header = summary_output.get_summary_header(
        {0: True, 1: True, 2: False}, graphs
    )
    assert "1/3 amplicons solved." in header


def test_amplicon_status_line() -> None:
    def status(bp_graph: BreakpointGraph, was_solved: bool | None) -> str:
        buf = io.StringIO()
        summary_output.output_amplicon_info(bp_graph, buf, was_solved)  # type: ignore[arg-type]
        return next(
            line
            for line in buf.getvalue().splitlines()
            if line.startswith("Cycle Decomposition Status")
        )

    # Walks are needed to print the heaviest walk; skip that part here.
    original = summary_output.output_amplicon_solution
    summary_output.output_amplicon_solution = lambda *_a, **_k: None  # type: ignore[assignment]
    try:
        assert status(_graph(0, explained=0.95), True).endswith("SUCCESS")
        assert status(
            _graph(0, explained=0.4, early_stop=EARLY_STOP), True
        ).endswith("PARTIAL")
        assert status(_graph(0, explained=None), False).endswith("FAILURE")
        assert status(_graph(0, explained=None), None).endswith("UNPERFORMED")
    finally:
        summary_output.output_amplicon_solution = original
