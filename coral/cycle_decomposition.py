"""Functions used for cycle decomposition"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

import numpy as np
import pyomo
import pyomo.contrib.appsi
import pyomo.core
import pyomo.environ
import pyomo.environ as pyo
import pyomo.opt
import pyomo.solvers
import pyomo.solvers.plugins
import pyomo.solvers.plugins.solvers
import pyomo.solvers.plugins.solvers.GUROBI
import pyomo.util.infeasible

from coral import datatypes, models
from coral.breakpoint import infer_breakpoint_graph
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.models import cycle_utils, output
from coral.models.concrete import CycleLPModel
from coral.models.path_constraints import longest_path_dict

logger = logging.getLogger(__name__)


def minimize_cycles(
    amplicon_id: int,
    bp_graph: BreakpointGraph,
    k: int,
    total_weights: float,
    node_order,
    pc_list,
    p_total_weight=0.9,
    p_bp_cn=0.9,
    num_threads=-1,
    time_limit=7200,
    model_prefix="pyomo",
    solver_to_use: datatypes.Solver = datatypes.Solver.GUROBI,
):
    """Cycle decomposition by minimizing the number of cycles/paths.

    Standard model for cycle detection used when generated BreakpointGraph is tractable.

    Args:
        amplicon_id: integer, amplicon ID
        g: breakpoint graph (object)
        k: integer, maximum mumber of cycles/paths allowed in cycle decomposition
        total_weights: float, total length-weighted CN in breakpoint graph g
        node_order: dict maps each node in the input breakpoint graph to a distinct integer, indicating a total order of the nodes in g
        pc_list: list of subpath constraints to be satisfied, each as a dict that maps an edge to its multiplicity
                *** note that all subpath constraints in this list are required to be satisfied ***
                *** otherwise will return infeasible ***
        p_total_weight: float between (0, 1), minimum proportion of length-weighted CN to be covered by the resulting cycles or paths,
                        default value is 0.9
        p_bp_cn: float float between (0, 1), minimum proportion of CN for each discordant edge to be covered by the resulting cycles or paths,
                        default value is 0.9
        num_threads: integer, number of working threads for gurobipy, by default it tries to use up all available cores
        time_limit: integer, maximum allowed running time, in seconds, default is 7200 (2 hour)
        model_prefix: output prefix for gurobi *.lp model

    Returns:
        A solution container parsed from the MIQCP solver output with the following attributes:
            (1) Pyomo status of optimization model
            (2) Total length weighted CN in resulting cycles/paths
            (3) Total num subpath constraints satisfied by resulting cycles/paths
            (4) List of cycles, each as a dict which maps an edge to its multiplicity in the cycle
            (5) List of the corresponding CN of the above cycles
            (6) Subpath constraints (indices) satisfied by each cycle
            (7) List of paths, each as a dict which maps an edge to its multiplicity in the path
            (8) List of the corresponding CN of the above paths
            (9) Subpath constraints (indices) satisfied by each path
    """
    logger.debug(
        f"Regular cycle decomposition with at most {k} cycles/paths allowed."
    )

    model = models.concrete.get_model(
        bp_graph,
        k,
        total_weights,
        node_order,
        pc_list,
        model_name=f"{model_prefix}/amplicon_{amplicon_id}_cycle_decomposition_{k=}",
    )

    model_name = f"{model_prefix}/amplicon_{amplicon_id}_model"
    model.write(f"{model_name}.lp", io_options={"symbolic_solver_labels": True})
    logger.debug(f"Completed model setup, wrote to {model_name}.lp.")

    solver = cycle_utils.get_solver(
        solver_type=solver_to_use,
        num_threads=num_threads,
        time_limit_s=max(
            time_limit, bp_graph.num_disc_edges * 300
        ),  # each breakpoint edge is assigned 5 minutes
    )
    results: pyomo.opt.SolverResults = solver.solve(model, tee=True)

    logger.debug(
        f"Completed optimization with status {results.solver.status}, condition {results.solver.termination_condition}."
    )
    if (
        results.solver.termination_condition
        == pyo.TerminationCondition.infeasible
    ):
        pyomo.util.infeasible.log_infeasible_constraints(
            model, log_expression=True, log_variables=True
        )

    return cycle_utils.parse_solver_output(
        results.solver.status,
        results.solver.termination_condition,
        model,
        bp_graph,
        k,
        pc_list,
        total_weights,
    )


def minimize_cycles_post(
    amplicon_id: int,
    bp_graph: BreakpointGraph,
    total_weights: float,
    node_order: Dict[tuple[Any, Any, Any], int],
    pc_list: list,
    init_sol: datatypes.InitialSolution,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    num_threads: int = -1,
    time_limit: int = 7200,
    model_prefix: str = "pyomo",
    solver_to_use: datatypes.Solver = datatypes.Solver.GUROBI,
) -> datatypes.CycleSolution:
    """Cycle decomposition by postprocessing the greedy MIQCP solution.

    Refines the initial solution generated by the greedy cycle model by
    minimizing the total # cycles out of the proposed cycles/paths.

    Args:
        g: breakpoint graph (object)
        total_weights: float, total length-weighted CN in breakpoint graph g
        node_order: dict maps each node in the input breakpoint graphg to a
            distinct integer, indicating a total order of the nodes in g
        pc_list: list of subpath constraints to be satisfied, each as a dict that
            maps an edge to its multiplicity
        init_sol: initial solution returned by maximize_weights_greedy
        p_total_weight: float between (0, 1), minimum proportion of length-weighted
            CN to be covered by the resulting cycles or paths, default value is 0.9
        resolution: float, minimum CN for each cycle or path, default value is 0.1
        num_threads: integer, number of working threads for gurobipy, by default it tries to use up all available cores
        time_limit: integer, maximum allowed running time, in seconds, default is 7200 (2 hour)
        model_prefix: output prefix for gurobi *.lp model

    Returns:
        A solution container parsed from the MIQCP solver output with the following attributes:
            (1) Pyomo status of optimization model
            (2) Total length weighted CN in resulting cycles/paths
            (3) Total num subpath constraints satisfied by resulting cycles/paths
            (4) List of cycles, each as a dict which maps an edge to its multiplicity in the cycle
            (5) List of the corresponding CN of the above cycles
            (6) Subpath constraints (indices) satisfied by each cycle
            (7) List of paths, each as a dict which maps an edge to its multiplicity in the path
            (8) List of the corresponding CN of the above paths
            (9) Subpath constraints (indices) satisfied by each path
    """
    logger.debug(
        "Cycle decomposition with initial solution from the greedy strategy."
    )

    k = len(init_sol.cycles[0]) + len(init_sol.cycles[1])
    logger.debug(f"Reset k (num cycles) to {k}.")
    p_path_constraints = 0.0
    path_constraint_indices_ = []
    for paths in (
        init_sol.satisfied_path_constraints[0]
        + init_sol.satisfied_path_constraints[1]
    ):
        for pathi in paths:
            if pathi not in path_constraint_indices_:
                path_constraint_indices_.append(pathi)
    if len(pc_list) > 0:
        p_path_constraints = (
            len(path_constraint_indices_) * 0.9999 / len(pc_list)
        )
        logger.debug(
            f"Required proportion of subpath constraints to be satisfied: {p_path_constraints}."
        )
    else:
        logger.debug("Proceed without subpath constraints.")

    model_name = f"{model_prefix}/amplicon_{amplicon_id}_cycle_decomposition_postprocessing_{k=}"
    model = models.concrete.get_model(
        bp_graph,
        k,
        total_weights,
        node_order,
        pc_list,
        model_name=model_name,
        is_post=True,
    )

    initialize_post_processing_solver(model, init_sol)
    model.write(f"{model_name}.lp", io_options={"symbolic_solver_labels": True})

    logger.debug(f"Completed model setup, wrote to {model_name}.lp.")
    solver = cycle_utils.get_solver(
        solver_type=solver_to_use,
        num_threads=num_threads,
        time_limit_s=max(
            time_limit, bp_graph.num_disc_edges * 300
        ),  # each breakpoint edge is assigned 5 minutes
    )
    results: pyomo.opt.SolverResults = solver.solve(model, tee=True)
    return cycle_utils.parse_solver_output(
        results.solver.status,
        results.solver.termination_condition,
        model,
        bp_graph,
        k,
        pc_list,
        total_weights,
    )


def initialize_post_processing_solver(
    model: CycleLPModel, init_sol: datatypes.InitialSolution
) -> None:
    """
    bp_graph: breakpoint graph (object)
    init_sol: initial solution returned by maximize_weights_greedy
    """
    for i in range(len(init_sol.cycles[0])):
        model.z[i] = 1
        model.w[i] = init_sol.cycle_weights[0][i]
        for var_name, var_idx in init_sol.cycles[0][i]:
            if var_name == "x":
                model.x[var_idx, i] = init_sol.cycles[0][i][(var_name, var_idx)]
            elif var_name == "c":
                model.c[var_idx, i] = init_sol.cycles[0][i][(var_name, var_idx)]
            elif var_name == "d":
                model.d[var_idx, i] = init_sol.cycles[0][i][(var_name, var_idx)]
            elif var_name == "y1":
                model.y1[var_idx, i] = init_sol.cycles[0][i][
                    (var_name, var_idx)
                ]
            elif var_name == "y2":
                model.y2[var_idx, i] = init_sol.cycles[0][i][
                    (var_name, var_idx)
                ]
    for i in range(len(init_sol.cycles[1])):
        i_ = i + len(init_sol.cycles[0])
        model.z[i_] = 1
        model.w[i_] = init_sol.cycle_weights[1][i]
        for v, vi in init_sol.cycles[1][i].keys():
            if v == "x":
                model.x[vi, i_] = init_sol.cycles[1][i][(v, vi)]
            elif v == "c":
                model.c[vi, i_] = init_sol.cycles[1][i][(v, vi)]
            elif v == "d":
                model.d[vi, i_] = init_sol.cycles[1][i][(v, vi)]
            elif v == "y1":
                model.y1[vi, i_] = init_sol.cycles[1][i][(v, vi)]
            elif v == "y2":
                model.y2[vi, i_] = init_sol.cycles[1][i][(v, vi)]
    for i in range(len(init_sol.satisfied_path_constraints[0])):
        for pi in init_sol.satisfied_path_constraints[0][i]:
            model.r[pi, i] = 1
            model.R[pi] = 1
    for i in range(len(init_sol.satisfied_path_constraints[1])):
        i_ = i + len(init_sol.satisfied_path_constraints[0])
        for pi in init_sol.satisfied_path_constraints[1][i]:
            model.r[pi, i_] = 1
            model.R[pi] = 1


def maximize_weights_greedy(
    amplicon_id: int,
    bp_graph: BreakpointGraph,
    total_weights: float,
    node_order: Dict[tuple[Any, Any, Any], int],
    pc_list: List,
    cycle_id: int,
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    cn_tol: float = 0.005,
    p_subpaths: float = 0.9,
    num_threads: int = -1,
    postprocess: int = 0,
    time_limit: int = 7200,
    model_prefix: str = "pyomo",
    solver_to_use: datatypes.Solver = datatypes.Solver.GUROBI,
):
    """Greedy cycle decomposition by maximizing the total length-weighted CN.

    The basis of this model can essentially be considered a case of the standard model
    used in `minimize_cycles` with fixed k = 1. The objective function of the
    solver maximizes total CN rather than minimizing cycles, and this is iteratively reduced

    Returns:
        A solution container parsed from the MIQCP solver output with the following attributes:
            (1) Pyomo status of optimization model
            (2) Total length weighted CN in resulting cycles/paths
            (3) Total num subpath constraints satisfied by resulting cycles/paths
            (4) List of cycles, each as a dict which maps an edge to its multiplicity in the cycle
            (5) List of the corresponding CN of the above cycles
            (6) Subpath constraints (indices) satisfied by each cycle
            (7) List of paths, each as a dict which maps an edge to its multiplicity in the path
            (8) List of the corresponding CN of the above paths
            (9) Subpath constraints (indices) satisfied by each path
    """

    remaining_weights = total_weights
    is_pc_unsatisfied = [True for i in range(len(pc_list))]
    num_unsatisfied_pc = sum(is_pc_unsatisfied)
    remaining_cn = datatypes.EdgeToCN.from_graph(bp_graph)
    next_w = resolution * 1.1

    full_solution = datatypes.CycleSolution(
        solver_status=pyo.SolverStatus.unknown,
        termination_condition=pyo.TerminationCondition.unknown,
    )

    logger.debug(
        f"Greedy cycle decomposition with length weighted CN = {remaining_weights} and num subpath constraints = {num_unsatisfied_pc}."
    )

    while next_w >= resolution and (
        remaining_weights > (1.0 - p_total_weight) * total_weights
        or num_unsatisfied_pc > np.floor((1.0 - p_subpaths) * len(pc_list))
    ):
        pp = 1.0
        if alpha > 0 and num_unsatisfied_pc > 0:
            pp = (
                alpha * remaining_weights / num_unsatisfied_pc
            )  # multi - objective optimization parameter
        logger.debug(
            f"Iteration {cycle_id + 1} with remaining CN = {remaining_weights} and num subpath constraints = {num_unsatisfied_pc}/{len(pc_list)}."
        )
        logger.debug(f"Multiplication factor for subpath constraints = {pp}.")
        model_name = f"{model_prefix}/amplicon_{amplicon_id}_cycle_decomposition_greedy_{cycle_id + 1}_{alpha=}"

        model = models.concrete.get_model(
            bp_graph,
            k=1,
            total_weights=total_weights,
            node_order=node_order,
            pc_list=pc_list,
            model_name=model_name,
            is_greedy=True,
            pp=pp,
            is_pc_unsatisfied=is_pc_unsatisfied,
            remaining_cn=remaining_cn,
        )
        model.write(
            f"{model_name}.lp", io_options={"symbolic_solver_labels": True}
        )
        logger.debug(f"Completed model setup, wrote to {model_name}.lp.")

        solver = cycle_utils.get_solver(
            solver_type=solver_to_use,
            num_threads=num_threads,
            time_limit_s=max(
                time_limit, bp_graph.num_disc_edges * 300
            ),  # each breakpoint edge is assigned 5 minutes
        )
        results: pyomo.opt.SolverResults = solver.solve(model, tee=True)
        curr_sol = cycle_utils.parse_solver_output(
            results.solver.status,
            results.solver.termination_condition,
            model,
            bp_graph,
            1,
            pc_list,
            total_weights,
            remaining_cn,
            resolution,
            is_pc_unsatisfied=is_pc_unsatisfied,
        )

        if (
            curr_sol.termination_condition
            == pyo.TerminationCondition.infeasible
        ):
            pyomo.util.infeasible.log_infeasible_constraints(
                model, log_expression=True, log_variables=True
            )
            logger.debug("Greedy cycle decomposition is infeasible, stopping.")
            break
        cycle_id += 1

        # Check if solver produced a cycle or path
        if len(curr_sol.walks.cycles) > 0:
            curr_walk_weight = curr_sol.walk_weights.cycles[0]
            full_solution.walks.cycles.append(curr_sol.walks.cycles[0])
            full_solution.walk_weights.cycles.append(curr_walk_weight)
            full_solution.satisfied_pc.cycles.append(
                curr_sol.satisfied_pc.cycles[0]
            )
        else:
            curr_walk_weight = curr_sol.walk_weights.paths[0]
            full_solution.walks.paths.append(curr_sol.walks.paths[0])
            full_solution.walk_weights.paths.append(curr_walk_weight)
            full_solution.satisfied_pc.paths.append(
                curr_sol.satisfied_pc.paths[0]
            )
        full_solution.satisfied_pc_set |= curr_sol.satisfied_pc_set

        # Update greedy stop conditions based on latest solution
        next_w = curr_walk_weight
        num_unsatisfied_pc = len(pc_list) - len(full_solution.satisfied_pc_set)
        remaining_weights -= curr_sol.total_weights_included


def cycle_decomposition(
    bb: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    num_threads: int = -1,
    postprocess: int = 0,
    time_limit: int = 7200,
    solver_to_use: datatypes.Solver = datatypes.Solver.GUROBI,
    *,
    output_all_path_constraints: bool = False,
    model_prefix: str = "pyomo",
) -> None:
    """Caller for cycle decomposition functions"""
    for amplicon_idx in range(len(bb.lr_graph)):
        bp_graph = bb.lr_graph[amplicon_idx]
        lseg = len(bp_graph.sequence_edges)
        lc = len(bp_graph.concordant_edges)
        ld = len(bp_graph.discordant_edges)
        lsrc = len(bp_graph.source_edges)

        total_weights = 0.0
        for sseg in bp_graph.sequence_edges:
            total_weights += sseg[7] * sseg[-1]  # type: ignore[operator]
        logger.info(f"Begin cycle decomposition for amplicon{amplicon_idx +1}.")
        logger.info(f"Total CN weights = {total_weights}.")

        bb.longest_path_constraints[amplicon_idx] = longest_path_dict(
            bb.path_constraints[amplicon_idx],
        )
        logger.info(
            f"Total num maximal subpath constraints = {len(bb.longest_path_constraints[amplicon_idx][0])}."
        )
        print(f"Solving {amplicon_idx}")
        for pathi in bb.longest_path_constraints[amplicon_idx][1]:
            logger.debug(
                f"Subpath constraint {pathi} = {bb.path_constraints[amplicon_idx][0][pathi]}"
            )

        k = max(10, ld // 2)  # Initial num cycles/paths
        logger.info(f"Initial num cycles/paths = {k}.")
        nnodes = len(bp_graph.nodes)  # Does not include s and t
        nedges = lseg + lc + ld + 2 * lsrc + 2 * len(bp_graph.endnodes)
        node_order = {}
        ni_ = 0
        for node in bb.lr_graph[amplicon_idx].nodes:
            node_order[node] = ni_
            ni_ += 1
        if nedges < k:
            k = nedges
            logger.info(f"Reset num cycles/paths to {k}.")
        sol_flag = 0
        while k <= nedges:
            if (
                nedges > 100
                or (
                    3 * k
                    + 3 * k * nedges
                    + 2 * k * nnodes
                    + k * len(bb.longest_path_constraints[amplicon_idx][0])
                )
                >= 10000
            ):
                (
                    total_cycle_weights_init,
                    total_path_satisfied_init,
                    cycles_init,
                    cycle_weights_init,
                    path_constraints_satisfied_init,
                ) = maximize_weights_greedy(
                    amplicon_id=amplicon_idx + 1,
                    bp_graph=bb.lr_graph[amplicon_idx],
                    total_weights=total_weights,
                    node_order=node_order,
                    pc_list=bb.longest_path_constraints[amplicon_idx][0],
                    cycle_id=0,  # TODO: FIX!
                    alpha=alpha,
                    p_total_weight=p_total_weight,
                    resolution=resolution,
                    cn_tol=0.005,
                    p_subpaths=0.9,
                    num_threads=num_threads,
                    postprocess=postprocess,
                    time_limit=time_limit,
                    model_prefix=model_prefix,
                    solver_to_use=solver_to_use,
                )
                logger.info("Completed greedy cycle decomposition.")
                logger.info(
                    f"Num cycles = {len(cycles_init[0])}; num paths = {len(cycles_init[1])}."
                )
                logger.info(
                    f"Total length weighted CN = {total_cycle_weights_init}/{total_weights}."
                )
                logger.info(
                    f"Total num subpath constraints satisfied = {total_path_satisfied_init}/{len(bb.longest_path_constraints[amplicon_idx][0])}."
                )
                if postprocess == 1:
                    lp_solution = minimize_cycles_post(
                        amplicon_idx + 1,
                        bb.lr_graph[amplicon_idx],
                        total_weights,
                        node_order,
                        bb.longest_path_constraints[amplicon_idx][0],
                        datatypes.InitialSolution(
                            cycles_init,
                            cycle_weights_init,
                            path_constraints_satisfied_init,
                        ),
                        min(
                            total_cycle_weights_init / total_weights * 0.9999,
                            p_total_weight,
                        ),
                        resolution,
                        num_threads,
                        time_limit,
                        model_prefix,
                        solver_to_use,
                    )

                    logger.info(
                        "Completed postprocessing of the greedy solution."
                    )
                    logger.info(
                        f"Num cycles = {len(lp_solution.walks[0])}; num paths = {len(lp_solution.walks[1])}."
                    )
                    logger.info(
                        f"Total length weighted CN = {lp_solution.total_weights_included}/{total_weights}."
                    )
                    logger.info(
                        f"Total num subpath constraints satisfied = {lp_solution.num_pc_satisfied}/{len(bb.longest_path_constraints[amplicon_idx][0])}."
                    )
                    bb.walks_by_amplicon[amplicon_idx] = lp_solution.walks
                    bb.walk_weights_by_amplicon[amplicon_idx] = (
                        lp_solution.walk_weights
                    )
                    bb.path_constraints_satisfied[amplicon_idx] = (
                        lp_solution.satisfied_pc
                    )
                else:
                    bb.walks_by_amplicon[amplicon_idx] = cycles_init
                    bb.walk_weights_by_amplicon[amplicon_idx] = (
                        cycle_weights_init
                    )
                    bb.path_constraints_satisfied[amplicon_idx] = (
                        path_constraints_satisfied_init
                    )
                sol_flag = 1
                break
            lp_solution = minimize_cycles(
                amplicon_idx + 1,
                bb.lr_graph[amplicon_idx],
                k,
                total_weights,
                node_order,
                bb.longest_path_constraints[amplicon_idx][0],
                p_total_weight,
                0.9,
                num_threads,
                time_limit,
                model_prefix,
                solver_to_use,
            )
            if (
                lp_solution.termination_condition
                == pyo.TerminationCondition.infeasible
            ):
                logger.info("Cycle decomposition is infeasible.")
                logger.info(f"Doubling k from {k} to {k * 2}.")
                k *= 2
            else:
                logger.info(f"Completed cycle decomposition with k = {k}.")
                logger.info(
                    f"Num cycles = {len(lp_solution.walks[0])}; num paths = {len(lp_solution.walks[1])}."
                )
                logger.info(
                    f"Total length weighted CN = {lp_solution.total_weights_included}/{total_weights}."
                )
                logger.info(
                    f"Total num subpath constraints satisfied = {len(lp_solution.satisfied_pc_set)}/{len(bb.longest_path_constraints[amplicon_idx][0])}."
                )

                bb.walks_by_amplicon[amplicon_idx] = lp_solution.walks
                bb.walk_weights_by_amplicon[amplicon_idx] = (
                    lp_solution.walk_weights
                )
                bb.path_constraints_satisfied[amplicon_idx] = (
                    lp_solution.satisfied_pc
                )
                sol_flag = 1
                break
        if sol_flag == 0:
            logger.info(
                "Cycle decomposition is infeasible, switch to greedy cycle decomposition."
            )
            (
                total_cycle_weights_init,
                total_path_satisfied_init,
                cycles_init,
                cycle_weights_init,
                path_constraints_satisfied_init,
            ) = maximize_weights_greedy(
                amplicon_id=amplicon_idx + 1,
                bp_graph=bb.lr_graph[amplicon_idx],
                total_weights=total_weights,
                node_order=node_order,
                pc_list=bb.longest_path_constraints[amplicon_idx][0],
                cycle_id=0,  # TODO: FIX!
                alpha=alpha,
                p_total_weight=p_total_weight,
                resolution=resolution,
                cn_tol=0.005,
                p_subpaths=0.9,
                num_threads=num_threads,
                postprocess=postprocess,
                time_limit=time_limit,
                model_prefix=model_prefix,
            )
            logger.info(
                "Completed greedy cycle decomposition.",
            )
            logger.info(
                "\tNum cycles = %d; num paths = %d."
                % (len(cycles_init[0]), len(cycles_init[1])),
            )
            logger.info(
                "\tTotal length weighted CN = %f/%f."
                % (total_cycle_weights_init, total_weights),
            )
            logger.info(
                "\tTotal num subpath constraints satisfied = %d/%d."
                % (
                    total_path_satisfied_init,
                    len(bb.longest_path_constraints[amplicon_idx][0]),
                ),
            )
            if postprocess == 1:
                lp_solution = minimize_cycles_post(
                    amplicon_idx + 1,
                    bb.lr_graph[amplicon_idx],
                    total_weights,
                    node_order,
                    bb.longest_path_constraints[amplicon_idx][0],
                    datatypes.InitialSolution(
                        cycles_init,
                        cycle_weights_init,
                        path_constraints_satisfied_init,
                    ),
                    min(
                        total_cycle_weights_init / total_weights * 0.9999,
                        p_total_weight,
                    ),
                    resolution,
                    num_threads,
                    time_limit,
                    model_prefix,
                    solver_to_use,
                )
                logger.info("Completed postprocessing of the greedy solution.")
                logger.info(
                    "Num cycles = %d; num paths = %d."
                    % (len(lp_solution.walks[0]), len(lp_solution.walks[1])),
                )
                logger.info(
                    "Total length weighted CN = %f/%f."
                    % (lp_solution.total_weights_included, total_weights),
                )
                logger.info(
                    "Total num subpath constraints satisfied = %d/%d."
                    % (
                        lp_solution.num_pc_satisfied,
                        len(bb.longest_path_constraints[amplicon_idx][0]),
                    ),
                )
                bb.walks_by_amplicon[amplicon_idx] = lp_solution.walks
                bb.walk_weights_by_amplicon[amplicon_idx] = (
                    lp_solution.walk_weights
                )
                bb.path_constraints_satisfied[amplicon_idx] = (
                    lp_solution.satisfied_pc
                )
            else:
                bb.walks_by_amplicon[amplicon_idx] = cycles_init
                bb.walk_weights_by_amplicon[amplicon_idx] = cycle_weights_init
                bb.path_constraints_satisfied[amplicon_idx] = (
                    path_constraints_satisfied_init
                )
        else:
            output.output_amplicon_cycles(
                amplicon_idx, bb, model_prefix, output_all_path_constraints
            )


def reconstruct_cycles(
    output_dir: str,
    output_all_path_constraints: bool,
    cycle_decomp_alpha: float,
    cycle_decomp_time_limit: int,
    cycle_decomp_threads: int,
    solver_to_use: datatypes.Solver,
    postprocess_greedy_sol: bool,
    bb: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
):
    logging.basicConfig(
        filename=f"{output_dir}/cycle_decomp.log",
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    bb.compute_path_constraints()
    logger.info("Computed all subpath constraints.")

    alpha_ = 0.01
    postprocess_ = 0
    nthreads = -1
    time_limit_ = 7200
    if cycle_decomp_alpha:
        alpha_ = cycle_decomp_alpha
    if postprocess_greedy_sol:
        postprocess_ = 1
    if cycle_decomp_threads:
        nthreads = cycle_decomp_threads
    if cycle_decomp_time_limit:
        time_limit_ = cycle_decomp_time_limit
    cycle_decomposition(
        bb,
        alpha=alpha_,
        num_threads=nthreads,
        postprocess=postprocess_,
        time_limit=time_limit_,
        solver_to_use=solver_to_use,
        output_all_path_constraints=output_all_path_constraints,
    )
    logger.info("Completed cycle decomposition for all amplicons.")
    logger.info(
        f"Wrote cycles for all complicons to {output_dir}/amplicon*_cycles.txt."
    )
