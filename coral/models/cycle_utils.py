from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional

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

from coral import datatypes, types
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import CycleSolution, EdgeToCN

logger = logging.getLogger(__name__)


def process_cycle_edge(
    cycle: types.AmpliconWalk,
    model: pyo.Model,
    edge_idx: int,
    edge_count: int,
    bp_graph: BreakpointGraph,
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
        cycle[("e", edge_idx)] = edge_count
        if remaining_cn:
            remaining_cn.sequence[edge_idx] -= edge_count * model.w[0].value
            if remaining_cn.sequence[edge_idx] < resolution:
                remaining_cn.sequence[edge_idx] = 0.0
    # Is concordant edge
    elif edge_idx < bp_graph.num_seq_edges + bp_graph.num_conc_edges:
        conc_edge_idx = edge_idx - bp_graph.num_seq_edges
        cycle[("c", conc_edge_idx)] = edge_count
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
        cycle[("d", disc_edge_idx)] = edge_count
        if remaining_cn:
            remaining_cn.discordant[disc_edge_idx] -= (
                edge_count * model.w[0].value
            )
            if remaining_cn.discordant[disc_edge_idx] < resolution:
                remaining_cn.discordant[disc_edge_idx] = 0.0
    # Is source edge
    elif edge_idx < src_node_offset:
        assert edge_count == 1
        src_edge_idx = edge_idx - bp_graph.num_nonsrc_edges
        if src_edge_idx % 2 == 0:
            s_edge_idx = src_edge_idx // 2
            cycle[("s", s_edge_idx)] = 1  # source edge connected to s
            if remaining_cn:
                remaining_cn.source[s_edge_idx] -= edge_count * model.w[0].value
                if remaining_cn.source[s_edge_idx] < resolution:
                    remaining_cn.source[s_edge_idx] = 0.0
        else:
            t_edge_idx = (src_edge_idx - 1) // 2
            cycle[("t", t_edge_idx)] = 1  # source edge connected to t
            if remaining_cn:
                remaining_cn.source[t_edge_idx] -= edge_count * model.w[0].value
                if remaining_cn.source[t_edge_idx] < resolution:
                    remaining_cn.source[t_edge_idx] = 0.0
    else:
        # Is synthetic end node
        assert edge_count == 1
        if (edge_idx - src_node_offset) % 2 == 0:
            nsi = (edge_idx - src_node_offset) // 2
            cycle[("ns", nsi)] = 1  # source edge connected to s
        else:
            nti = (edge_idx - src_node_offset - 1) // 2
            cycle[("nt", nti)] = 1  # source edge connected to t


def parse_solver_output(
    solver_status: str,
    solver_termination_condition: pyo.TerminationCondition,
    model: pyo.Model,
    bp_graph: BreakpointGraph,
    k: int,
    pc_list: List,
    total_weights: float,
    remaining_cn: Optional[EdgeToCN] = None,
    resolution: float = 0.0,
    is_pc_unsatisfied: List[bool] | None = None,
) -> CycleSolution:
    parsed_sol = CycleSolution(solver_status, solver_termination_condition)

    lseg = len(bp_graph.sequence_edges)
    lc = len(bp_graph.concordant_edges)
    ld = len(bp_graph.discordant_edges)
    lsrc = len(bp_graph.source_edges)
    nnodes = len(bp_graph.nodes)  # Does not include s and t
    nedges = lseg + lc + ld + 2 * lsrc + 2 * len(bp_graph.endnodes)

    for i in range(k):
        logger.debug(f"Walk {i} checking ; CN = {model.w[i].value}.")
        if model.z[i].value >= 0.9:
            logger.debug(f"Cycle/Path {i} exists; CN = {model.w[i].value}.")
            if resolution and (walk_weight := model.w[i].value) < resolution:
                parsed_sol.walk_weights[0].append(walk_weight)
                logger.debug(
                    "\tCN less than resolution, iteration terminated successfully."
                )
                break
            found_cycle = False
            for node_idx in range(nnodes):
                if model.c[node_idx, i].value >= 0.9:
                    found_cycle = True
                    break
            if not found_cycle:
                cycle: dict = {}
                for edge_idx in range(nedges):
                    if (edge_count := model.x[edge_idx, i].value) >= 0.9:
                        edge_count = round(edge_count)
                        # Update cycle in-place via helper
                        process_cycle_edge(
                            cycle,
                            model,
                            edge_idx,
                            edge_count,
                            bp_graph,
                            remaining_cn,
                            resolution,
                        )
                path_constraints_s = [
                    pi
                    for pi in range(len(pc_list))
                    if model.r[pi, i].value >= 0.9
                ]
                if (walk_weight := model.w[i].value) > 0.0:
                    parsed_sol.walks[1].append(cycle)
                    parsed_sol.walk_weights[1].append(walk_weight)
                    parsed_sol.satisfied_pc[1].append(path_constraints_s)
                    parsed_sol.satisfied_pc_set |= set(path_constraints_s)
            else:
                cycle = {}
                path_constraints_s = []
                for edge_idx in range(nedges):
                    if (edge_count := model.x[edge_idx, i].value) >= 0.9:
                        edge_count = round(edge_count)
                        if edge_idx < lseg:
                            cycle[("e", edge_idx)] = edge_count
                        elif edge_idx < lseg + lc:
                            cycle[("c", edge_idx - lseg)] = edge_count
                        elif edge_idx < lseg + lc + ld:
                            cycle[("d", edge_idx - lseg - lc)] = edge_count
                        else:
                            logger.debug(
                                "Error: Cyclic path cannot connect to source nodes."
                            )
                            logger.debug("Aborted.")
                            os.abort()
                for pi in range(len(pc_list)):
                    if model.r[pi, i].value >= 0.9:
                        path_constraints_s.append(pi)

                        # Only used for greedy, flip to False when PC satisfied
                        if is_pc_unsatisfied:
                            is_pc_unsatisfied[pi] = False
                if (walk_weight := model.w[i].value) > 0.0:
                    parsed_sol.walks[0].append(cycle)
                    parsed_sol.walk_weights[0].append(walk_weight)
                    parsed_sol.satisfied_pc[0].append(path_constraints_s)
                    parsed_sol.satisfied_pc_set |= set(path_constraints_s)
            for seqi in range(lseg):
                parsed_sol.total_weights_included += (
                    model.x[seqi, i].value
                    * model.w[i].value
                    * bp_graph.sequence_edges[seqi][-2]
                )

    logger.debug(
        f"Total length weighted CN from cycles/paths = {parsed_sol.total_weights_included}/{total_weights}."
    )
    logger.debug(
        f"Total num subpath constraints satisfied = {len(parsed_sol.satisfied_pc_set)}/{len(pc_list)}."
    )
    return parsed_sol


def get_solver(
    solver_type: datatypes.Solver, num_threads: int, time_limit_s: int
) -> pyomo.solvers.plugins.solvers:
    solver = pyo.SolverFactory(solver_type.value)
    if solver_type == datatypes.Solver.GUROBI:
        solver = pyo.SolverFactory(solver_type.value)
        if num_threads > 0:
            solver.options["threads"] = num_threads
        solver.options["NonConvex"] = 2
        solver.options["timelimit"] = time_limit_s
    elif solver_type == datatypes.Solver.SCIP:
        solver.options["lp/threads"] = num_threads
        solver.options["parallel/maxnthreads"] = num_threads
        solver.options["alg/concurrent"] = 1
        solver.options["display/freq"] = 1
        solver.options["display/lpinfo"] = True
        solver.options["propagating/nlobbt/nlptimelimit"] = time_limit_s
    return solver
