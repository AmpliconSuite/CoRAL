from __future__ import annotations

import logging
import os
from typing import List, Optional

import pyomo
import pyomo.contrib.appsi
import pyomo.core
import pyomo.environ as pyo
import pyomo.opt
import pyomo.solvers
import pyomo.solvers.plugins
import pyomo.solvers.plugins.solvers
import pyomo.solvers.plugins.solvers.GUROBI
import pyomo.util.infeasible

from coral import datatypes
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    CycleSolution,
    EdgeId,
    EdgeToCN,
    EdgeType,
    OptimizationWalk,
)

logger = logging.getLogger(__name__)


def process_walk_edge(
    walk: OptimizationWalk,
    model: pyo.Model,
    edge_idx: int,
    edge_count: int,
    bp_graph: BreakpointGraph,
    is_cycle: bool,
    remaining_cn: Optional[EdgeToCN] = None,
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
    solver_status: str,
    solver_termination_condition: pyo.TerminationCondition,
    model: pyo.Model,
    bp_graph: BreakpointGraph,
    k: int,
    pc_list: List,
    total_weights: float,
    remaining_cn: EdgeToCN | None = None,
    resolution: float = 0.0,
    is_pc_unsatisfied: List[bool] | None = None,
) -> CycleSolution:
    parsed_sol = CycleSolution(solver_status, solver_termination_condition)

    if solver_termination_condition == pyo.TerminationCondition.infeasible:
        pyomo.util.infeasible.log_infeasible_constraints(
            model, log_expression=True, log_variables=True
        )
        logger.debug("Unable to parse infeasible solution.")
        return parsed_sol

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
            parsed_sol.walk_weights[0].append(walk_weight)
            logger.debug(
                "\tCN < resolution, iteration terminated successfully."
            )
            break
        found_cycle = False
        for node_idx in range(nnodes):
            if model.c[node_idx, i].value >= 0.9:
                found_cycle = True
                break
        cycle: dict = {}
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
                    cycle,
                    model,
                    edge_idx,
                    edge_count,
                    bp_graph,
                    is_cycle=found_cycle,
                    remaining_cn=remaining_cn,
                    resolution=resolution,
                )
        if (walk_weight := model.w[i].value) > 0.0:
            if not found_cycle:
                parsed_sol.walks[1].append(cycle)
                parsed_sol.walk_weights[1].append(walk_weight)
                parsed_sol.satisfied_pc[1].append(path_constraints_s)
                parsed_sol.satisfied_pc_set |= set(path_constraints_s)
            else:
                parsed_sol.walks[0].append(cycle)
                parsed_sol.walk_weights[0].append(walk_weight)
                parsed_sol.satisfied_pc[0].append(path_constraints_s)
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


def get_solver(
    solver_options: datatypes.SolverOptions,
) -> pyomo.solvers.plugins.solvers:
    solver_type = solver_options.solver
    solver = pyo.SolverFactory(solver_type.value)
    if solver_type == datatypes.Solver.GUROBI:
        if solver_options.num_threads > 0:
            solver.options["threads"] = solver_options.num_threads
        solver.options["NonConvex"] = 2
        solver.options["timelimit"] = solver_options.time_limit_s
    elif solver_type == datatypes.Solver.SCIP:
        solver.options["lp/threads"] = solver_options.num_threads
        solver.options["parallel/maxnthreads"] = solver_options.num_threads
        solver.options["limits/time"] = solver_options.time_limit_s
        solver.options["display/freq"] = 1
        solver.options["display/lpinfo"] = True
        solver.options["propagating/nlobbt/nlptimelimit"] = (
            solver_options.time_limit_s
        )
    return solver
