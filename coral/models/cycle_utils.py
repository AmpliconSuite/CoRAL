from __future__ import annotations

import logging
import os
from typing import Any, List, Optional

import pyomo
import pyomo.contrib.appsi
import pyomo.core
import pyomo.environ as pyo
import pyomo.opt
import pyomo.opt.base
import pyomo.opt.results.solver
import pyomo.solvers
import pyomo.solvers.plugins
import pyomo.solvers.plugins.solvers
import pyomo.solvers.plugins.solvers.GUROBI
import pyomo.util.infeasible

from coral import datatypes, global_state
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    CycleSolution,
    EdgeId,
    EdgeToCN,
    EdgeType,
    OptimizationWalk,
)
from coral.pyomo_wrappers import PyomoResults

logger = logging.getLogger(__name__)


def process_walk_edge(
    walk: OptimizationWalk,
    model: pyo.Model,
    edge_idx: int,
    edge_count: int,
    bp_graph: BreakpointGraph,
    is_cycle: bool,
    remaining_cn: EdgeToCN | None = None,
    resolution: float = 0.0,
) -> None:
    """Update `cycle` parameter with the appropriate entry for a given edge
    count after solving the LP."""
    src_node_offset = bp_graph.num_nonsrc_edges + 2 * bp_graph.num_src_edges

    if remaining_cn:
        assert resolution, "Resolution must be provided when processing a greedy solution using remaining CN."

    # Is sequence edge
    if edge_idx < bp_graph.num_seq_edges:
        walk[EdgeId(EdgeType.SEQUENCE, edge_idx)] = edge_count
        if remaining_cn:
            remaining_cn.sequence[edge_idx] -= edge_count * model.w[0].value
            if remaining_cn.sequence[edge_idx] < resolution:
                remaining_cn.sequence[edge_idx] = 0.0
    # Is concordant edge
    elif edge_idx < bp_graph.num_seq_edges + bp_graph.num_conc_edges:
        conc_edge_idx = edge_idx - bp_graph.num_seq_edges
        walk[EdgeId(EdgeType.CONCORDANT, conc_edge_idx)] = edge_count
        if remaining_cn:
            remaining_cn.concordant[conc_edge_idx] -= (
                edge_count * model.w[0].value
            )
            if remaining_cn.concordant[conc_edge_idx] < resolution:
                remaining_cn.concordant[conc_edge_idx] = 0.0
    # Is discordant edge
    elif edge_idx < bp_graph.num_nonsrc_edges:
        disc_edge_idx = (
            edge_idx - bp_graph.num_seq_edges - bp_graph.num_conc_edges
        )
        walk[EdgeId(EdgeType.DISCORDANT, disc_edge_idx)] = edge_count
        if remaining_cn:
            remaining_cn.discordant[disc_edge_idx] -= (
                edge_count * model.w[0].value
            )
            if remaining_cn.discordant[disc_edge_idx] < resolution:
                remaining_cn.discordant[disc_edge_idx] = 0.0
    elif is_cycle:
        logger.error("Cyclic path cannot connect to source nodes.")
        logger.warning("Aborted.")
        os.abort()
    # Is source edge
    elif edge_idx < src_node_offset:
        src_edge_idx = edge_idx - bp_graph.num_nonsrc_edges
        if src_edge_idx % 2 == 0:
            s_edge_idx = src_edge_idx // 2
            # source edge connected to s
            walk[EdgeId(EdgeType.SOURCE, s_edge_idx)] = 1
            if remaining_cn:
                remaining_cn.source[s_edge_idx] -= edge_count * model.w[0].value
                if remaining_cn.source[s_edge_idx] < resolution:
                    remaining_cn.source[s_edge_idx] = 0.0
        else:
            t_edge_idx = (src_edge_idx - 1) // 2
            # source edge connected to t
            walk[EdgeId(EdgeType.SINK, t_edge_idx)] = 1
            if remaining_cn:
                remaining_cn.source[t_edge_idx] -= edge_count * model.w[0].value
                if remaining_cn.source[t_edge_idx] < resolution:
                    remaining_cn.source[t_edge_idx] = 0.0
    # Is synthetic end node
    elif (edge_idx - src_node_offset) % 2 == 0:
        nsi = (edge_idx - src_node_offset) // 2
        # source edge connected to s
        walk[EdgeId(EdgeType.SYNTHETIC_SOURCE, nsi)] = 1
    else:
        nti = (edge_idx - src_node_offset - 1) // 2
        # source edge connected to t
        walk[EdgeId(EdgeType.SYNTHETIC_SINK, nti)] = 1


def parse_solver_output(
    results: PyomoResults,
    model: pyo.Model,
    bp_graph: BreakpointGraph,
    k: int,
    pc_list: List,
    total_weights: float,
    model_type: datatypes.ModelType,
    remaining_cn: EdgeToCN | None = None,
    resolution: float = 0.0,
    is_pc_unsatisfied: List[bool] | None = None,
    alpha: float | None = None,
) -> CycleSolution:
    termination_condition = results.solver.termination_condition
    parsed_sol = CycleSolution(
        results.solver.status,
        termination_condition,
        model_metadata=datatypes.ModelMetadata(
            model_type=model_type,
            k=k,
            alpha=alpha,
            total_weights=total_weights,
            resolution=resolution,
            num_path_constraints=len(pc_list),
        ),
    )

    if termination_condition == pyo.TerminationCondition.infeasible:
        # pyomo.util.infeasible.log_infeasible_constraints(
        #     model, log_expression=True, log_variables=True
        # )
        logger.warning(
            f"{bp_graph.amplicon_idx}: Unable to parse infeasible solution."
        )
        return parsed_sol
    # SCIP can return an infeasible pseudo-solution when the time limit is used
    # to interrupt the solver, so we need to check that the solution is actually
    # feasible.
    if (
        termination_condition == pyo.TerminationCondition.maxTimeLimit
        and results.solver.status == pyo.SolverStatus.ok
        and len(model.solutions) != 0
        and results.solver.gap == float("inf")
    ):
        logger.warning(
            f"{bp_graph.amplicon_idx}: Solver returned infeasible solution."
        )
        parsed_sol.solver_status = pyo.SolverStatus.aborted
        return parsed_sol
    if (
        results.solver.status == pyo.SolverStatus.aborted
        and termination_condition == pyo.TerminationCondition.maxTimeLimit
    ):
        if len(model.solutions) == 0:
            logger.warning(
                f"{bp_graph.amplicon_idx}: Solver reached time limit without "
                "a feasible solution."
            )
            return parsed_sol

        parsed_sol.solver_status = pyo.SolverStatus.ok
        logger.warning(
            f"{bp_graph.amplicon_idx}: Solver reached time limit with a feasible "
            "but suboptimal solution."
        )

    try:
        parsed_sol.upper_bound = model.Objective()
    except ValueError:
        # Gurobi can return an infeasible solution when the time limit is used
        # to terminate the solver.
        parsed_sol.solver_status = pyo.SolverStatus.aborted
        return parsed_sol

    try:
        # SCIP stores MIP gap in results, doesn't populate model gap
        parsed_sol.mip_gap = results.solver.gap
    except AttributeError:
        # Gurobi_Direct (gurobipy) stores MIP gap only in model
        parsed_sol.mip_gap = model.solutions[0].gap

    lseg = len(bp_graph.sequence_edges)
    lc = len(bp_graph.concordant_edges)
    ld = len(bp_graph.discordant_edges)
    lsrc = len(bp_graph.source_edges)
    nnodes = len(bp_graph.node_adjacencies)  # Does not include s and t
    nedges = lseg + lc + ld + 2 * lsrc + 2 * len(bp_graph.endnode_adjacencies)

    for i in range(k):
        if model.z[i].value < 0.9:
            logger.debug(f"Walk {i} does not exist; CN = {model.w[i].value}.")
            continue

        logger.debug(f"Walk {i} exists; CN = {model.w[i].value}.")
        if resolution and (walk_weight := model.w[i].value) < resolution:
            parsed_sol.walk_weights.cycles.append(walk_weight)
            logger.debug(
                "\tCN < resolution, iteration terminated successfully."
            )
            break
        found_cycle = False
        for node_idx in range(nnodes):
            if model.c[node_idx, i].value >= 0.9:
                found_cycle = True
                break
        walk_edge: OptimizationWalk = {}
        path_constraints_s = []
        for pi in range(len(pc_list)):
            if model.r[pi, i].value >= 0.9:
                path_constraints_s.append(pi)
                # Only used for greedy, flip to False when PC satisfied
                if is_pc_unsatisfied:
                    is_pc_unsatisfied[pi] = False
        for edge_idx in range(nedges):
            if (edge_count := model.x[edge_idx, i].value) >= 0.9:
                edge_count = round(edge_count)
                # Update cycle in-place via helper
                process_walk_edge(
                    walk_edge,
                    model,
                    edge_idx,
                    edge_count,
                    bp_graph,
                    is_cycle=found_cycle,
                    remaining_cn=remaining_cn,
                    resolution=resolution,
                )
        if (walk_weight := model.w[i].value) > 0.0:
            if not walk_edge:
                logger.error(f"Walk {i} has no edges.")
                continue
            if found_cycle:
                parsed_sol.walks.cycles.append(walk_edge)
                parsed_sol.walk_weights.cycles.append(walk_weight)
                parsed_sol.satisfied_pc.cycles.append(path_constraints_s)
                parsed_sol.satisfied_pc_set |= set(path_constraints_s)
            else:
                parsed_sol.walks.paths.append(walk_edge)
                parsed_sol.walk_weights.paths.append(walk_weight)
                parsed_sol.satisfied_pc.paths.append(path_constraints_s)
                parsed_sol.satisfied_pc_set |= set(path_constraints_s)
        for seqi in range(lseg):
            parsed_sol.total_weights_included += (
                model.x[seqi, i].value
                * model.w[i].value
                * bp_graph.sequence_edges[seqi].gap
            )

    logger.debug(
        f"Total length weighted CN from cycles/paths = {parsed_sol.total_weights_included}/{total_weights}."
    )
    logger.debug(
        f"Total num subpath constraints satisfied = {len(parsed_sol.satisfied_pc_set)}/{len(pc_list)}."
    )
    return parsed_sol


class PyomoSolverWrapper:
    def __init__(self, solver: pyomo.opt.base.OptSolver) -> None:
        self.solver: pyomo.opt.base.OptSolver = solver

    def solve(
        self, model: pyo.Model, model_filepath: str, **kwargs: dict[str, Any]
    ) -> PyomoResults:
        try:
            # Gurobi_Direct plugin doesn't respect the `logfile` option, so we
            # need to set it manually, and require `tee`.
            # SCIP / Gurobi Shell just need the `logfile` option.
            self.solver.options["LogFile"] = f"{model_filepath}.log"
            results = self.solver.solve(
                model,
                **kwargs,
                logfile=f"{model_filepath}.log",
                tee=True,
                # keepfiles=True,
                # load_solutions=True,
            )
            # We potentially need to load because the Pyomo flag will delete the
            # solver results object, which we need to extract MIP Gap
            # information for certain solvers where the Pyomo log-parsing fails.
            # model.solutions.load_from(results)
        except ValueError as e:
            logger.error(f"Solver {self.solver.name} raised an exception: {e}")
            failed_solver_info = pyomo.opt.results.solver.SolverInformation()
            failed_solver_info.status = pyomo.opt.SolverStatus.aborted
            failed_solver_info.termination_condition = (
                pyomo.opt.TerminationCondition.maxTimeLimit
            )

            failed_solver_results = pyomo.opt.SolverResults()
            failed_solver_results.solver = failed_solver_info

            return failed_solver_results
        return results


def get_solver(
    solver_options: datatypes.SolverOptions,
) -> PyomoSolverWrapper:
    solver_type = solver_options.solver
    raw_solver = pyo.SolverFactory(solver_type.value)
    wrapped_solver = PyomoSolverWrapper(raw_solver)

    solver_time_limit = min(
        global_state.STATE_PROVIDER.remaining_time_s,
        solver_options.time_limit_s,
    )

    if solver_type == datatypes.Solver.GUROBI:
        if solver_options.num_threads > 0:
            raw_solver.options["threads"] = solver_options.num_threads
        raw_solver.options["NonConvex"] = 2
        raw_solver.options["timelimit"] = solver_time_limit
    elif solver_type == datatypes.Solver.SCIP:
        raw_solver.options["limits/time"] = solver_time_limit
        raw_solver.options["display/freq"] = 1
        raw_solver.options["display/lpinfo"] = True

        # Threads are only respected when using "Fiber-SCIP" executable
        # raw_solver.options["lp/threads"] = solver_options.num_threads
        # solver.options["parallel/maxnthreads"] = solver_options.num_threads
        # solver.options["propagating/nlobbt/nlptimelimit"] = solver_time_limit
    return wrapped_solver
