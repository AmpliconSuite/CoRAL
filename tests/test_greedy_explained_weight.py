"""Tests for the CN-weighted length explained by greedy cycle decomposition.

The greedy loop subtracted each walk's CN-weighted length from its running
remainder but never added it to the solution it returns, so
``total_weights_included`` was always 0. ``--postprocess-greedy-sol`` uses that
value as the minimum coverage its re-solve must keep, so with 0 the re-solve
was free to drop every walk, and did.
"""

from __future__ import annotations

import pathlib
import shutil
from collections.abc import Iterator
from types import SimpleNamespace
from typing import Any

import pyomo.environ as pyo
import pytest

from coral import cycle_decomposition, datatypes, global_state
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.models import cycle_utils

TEST_GRAPH = pathlib.Path(__file__).parent / "data" / "amplicon1_graph.txt"

requires_scip = pytest.mark.skipif(
    shutil.which("scip") is None, reason="SCIP executable not on PATH"
)


@pytest.fixture
def ample_time() -> Iterator[None]:
    provider = global_state.STATE_PROVIDER
    saved = (provider.time_limit_s, provider.start_time)
    provider.time_limit_s = 3600
    yield
    provider.time_limit_s, provider.start_time = saved


def _load_test_graph() -> BreakpointGraph:
    with TEST_GRAPH.open() as graph_file:
        bp_graph = parse_breakpoint_graph(graph_file)
    bp_graph.amplicon_idx = 0
    bp_graph.infer_amplicon_intervals(0)
    return bp_graph


def _total_weights(bp_graph: BreakpointGraph) -> float:
    return sum(len(s) * s.cn for s in bp_graph.sequence_edges)


def _node_order(bp_graph: BreakpointGraph) -> dict[datatypes.Node, int]:
    return {node: i for i, node in enumerate(bp_graph.node_adjacencies)}


def _solver_options(tmp_path: pathlib.Path) -> datatypes.SolverOptions:
    (tmp_path / "models").mkdir(parents=True, exist_ok=True)
    return datatypes.SolverOptions(
        output_dir=tmp_path,
        num_threads=1,
        time_limit_s=60,
        solver=datatypes.Solver.SCIP,
    )


def _one_walk_solution(
    weight_included: float, *, is_cycle: bool
) -> datatypes.CycleSolution:
    sol = datatypes.CycleSolution(
        pyo.SolverStatus.ok,
        pyo.TerminationCondition.optimal,
        total_weights_included=weight_included,
    )
    walks = sol.walks.cycles if is_cycle else sol.walks.paths
    weights = sol.walk_weights.cycles if is_cycle else sol.walk_weights.paths
    pcs = sol.satisfied_pc.cycles if is_cycle else sol.satisfied_pc.paths
    walks.append({})
    weights.append(2.0)
    pcs.append([])
    return sol


def test_greedy_accumulates_explained_weight_without_solver(
    ample_time: None, monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    """Each greedy iteration's covered weight is summed into the result."""
    per_iteration = [600.0, 250.0, 100.0]  # of a 1000.0 total
    canned = iter(
        [
            _one_walk_solution(per_iteration[0], is_cycle=True),
            _one_walk_solution(per_iteration[1], is_cycle=False),
            _one_walk_solution(per_iteration[2], is_cycle=True),
        ]
    )

    def fake_parse(*_args: Any, **_kwargs: Any) -> datatypes.CycleSolution:
        return next(canned)

    monkeypatch.setattr(
        cycle_decomposition.models.concrete,
        "generate_model",
        lambda *_a, **_k: None,
    )
    monkeypatch.setattr(
        cycle_utils,
        "get_solver",
        lambda *_a, **_k: SimpleNamespace(solve=lambda *_a, **_k: None),
    )
    monkeypatch.setattr(cycle_utils, "parse_solver_output", fake_parse)

    bp_graph = SimpleNamespace(
        amplicon_idx=0,
        longest_path_constraints=[],
        sequence_edges=[],
        concordant_edges=[],
        discordant_edges=[],
        source_edges=[],
    )
    monkeypatch.setattr(
        datatypes.EdgeToCN, "from_graph", staticmethod(lambda _g: None)
    )

    sol = cycle_decomposition.maximize_weights_greedy(
        bp_graph=bp_graph,  # type: ignore[arg-type]
        total_weights=1000.0,
        node_order={},
        cycle_id=0,
        solver_options=_solver_options(tmp_path),
    )

    # Stops once >= 90% is explained (600 + 250 + 100 = 950).
    assert sol.num_cycles == 2
    assert sol.num_paths == 1
    assert sol.total_weights_included == pytest.approx(sum(per_iteration))


@requires_scip
def test_greedy_reports_explained_weight(
    ample_time: None, tmp_path: pathlib.Path
) -> None:
    bp_graph = _load_test_graph()
    total = _total_weights(bp_graph)
    sol = cycle_decomposition.maximize_weights_greedy(
        bp_graph=bp_graph,
        total_weights=total,
        node_order=_node_order(bp_graph),
        cycle_id=0,
        solver_options=_solver_options(tmp_path),
    )
    assert sol.num_cycles + sol.num_paths >= 1
    assert 0.9 <= sol.total_weights_included / total <= 1.0 + 1e-6


@requires_scip
def test_postprocess_keeps_greedy_walks(
    ample_time: None, tmp_path: pathlib.Path
) -> None:
    bp_graph = _load_test_graph()
    total = _total_weights(bp_graph)
    kwargs = dict(
        bp_graph=bp_graph,
        total_weights=total,
        node_order=_node_order(bp_graph),
        cycle_id=0,
        solver_options=_solver_options(tmp_path),
    )
    greedy = cycle_decomposition.greedy_solve(**kwargs)
    post = cycle_decomposition.greedy_solve(**kwargs, should_postprocess=True)

    assert greedy.num_cycles + greedy.num_paths >= 1
    assert (
        post.num_cycles + post.num_paths >= 1
    ), "postprocessing dropped every greedy walk"
    assert post.total_weights_included >= 0.9999 * min(
        greedy.total_weights_included, 0.9 * total
    )
