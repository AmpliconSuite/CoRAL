"""Functions used for cycle decomposition"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional

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
    solver_options: Optional[datatypes.SolverOptions] = None,
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
        output_dir: output directory for cycle decomposition model
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
    if not solver_options:
        solver_options = datatypes.SolverOptions()
    solver_options.time_limit_s = max(
        solver_options.time_limit_s, bp_graph.num_disc_edges * 300
    )  # each breakpoint edge is assigned >= 5 minutes)

    model_name = f"amplicon_{amplicon_id}_cycle_decomposition_{k=}"
    model_filepath = f"{solver_options.output_dir}/{solver_options.model_prefix}_{model_name}"

    model = models.concrete.get_model(
        bp_graph,
        k,
        total_weights,
        node_order,
        pc_list,
        model_name=model_name,
    )

    model.write(
        f"{model_filepath}.lp", io_options={"symbolic_solver_labels": True}
    )
    model.write(f"{model_filepath}_ampl.nl", format="nl")
    logger.debug(f"Completed model setup, wrote to {model_filepath}.lp.")

    solver = cycle_utils.get_solver(solver_options)
    results: pyomo.opt.SolverResults = solver.solve(model, tee=True)

    logger.debug(
        f"Completed optimization with status {results.solver.status}, condition {results.solver.termination_condition}."
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
    solver_options: Optional[datatypes.SolverOptions] = None,
) -> datatypes.CycleSolution:
    """Cycle decomposition by postprocessing the greedy MIQCP solution.

    Refines the initial solution generated by the greedy cycle model by
    minimizing the total # cycles out of the proposed cycles/paths.

    Args:
        g: breakpoint graph (object)
        total_weights: float, total length-weighted CN in breakpoint graph g
        node_order: dict maps each node in the input breakpoint graphg to a
            distinct integer, indicating a total order of the nodes in g
        pc_list: list of subpath constraints to be satisfied, each as a dict
            that maps an edge to its multiplicity
        init_sol: initial solution returned by maximize_weights_greedy
        p_total_weight: float between (0, 1), minimum proportion of
            length-weighted CN to be covered by the resulting cycles or paths,
            default value is 0.9
        resolution: float, minimum CN for each cycle or path, default value is
            0.1
        num_threads: integer, number of working threads for gurobipy, by default
            it tries to use up all available cores
        time_limit: integer, maximum allowed running time, in seconds,
            default is 7200 (2 hour)
        model_prefix: output prefix for gurobi *.lp model

    Returns:
        A solution container parsed from the MIQCP solver output with
        the following attributes:
            (1) Pyomo status of optimization model
            (2) Total length weighted CN in resulting cycles/paths
            (3) Total num subpath constraints satisfied by resulting
                cycles/paths
            (4) List of cycles, each as a dict which maps an edge to its
                multiplicity in the cycle
            (5) List of the corresponding CN of the above cycles
            (6) Subpath constraints (indices) satisfied by each cycle
            (7) List of paths, each as a dict which maps an edge to its
                multiplicity in the path
            (8) List of the corresponding CN of the above paths
            (9) Subpath constraints (indices) satisfied by each path
    """
    logger.debug(
        "Cycle decomposition with initial solution from the greedy strategy."
    )
    if not solver_options:
        solver_options = datatypes.SolverOptions()
    solver_options.time_limit_s = max(
        solver_options.time_limit_s, bp_graph.num_disc_edges * 300
    )  # each breakpoint edge is assigned >= 5 minutes)

    k = len(init_sol.walks[0]) + len(init_sol.walks[1])
    logger.debug(f"Reset k (num cycles) to {k}.")
    p_path_constraints = 0.0
    path_constraint_indices_ = []
    for paths in init_sol.satisfied_pc[0] + init_sol.satisfied_pc[1]:
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

    model_name = (
        f"amplicon_{amplicon_id}_cycle_decomposition_postprocessing_{k=}"
    )
    model_filepath = f"{solver_options.output_dir}/{solver_options.model_prefix}_{model_name}"

    model = models.concrete.get_model(
        bp_graph,
        k,
        total_weights,
        node_order,
        pc_list,
        model_name=model_name,
        is_post=True,
        resolution=resolution,
        p_total_weight=p_total_weight,
    )

    initialize_post_processing_solver(model, init_sol)
    model.write(
        f"{model_filepath}.lp", io_options={"symbolic_solver_labels": True}
    )
    model.write(f"{model_filepath}_ampl.nl", format="nl")

    logger.debug(f"Completed model setup, wrote to {model_filepath}.lp.")
    solver = cycle_utils.get_solver(solver_options)
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
    """Initialize solver based on values generated by greedy decomposition.

    Args:
        bp_graph: breakpoint graph (object)
        init_sol: initial solution returned by maximize_weights_greedy

    Returns:
        None - modifies the input model object directly with the provided
        solution's values for each variable.
    """
    for i in range(len(init_sol.walks[0])):
        model.z[i] = 1
        model.w[i] = init_sol.walk_weights[0][i]
        for var_name, var_idx in init_sol.walks[0][i]:
            if var_name == "x":
                model.x[var_idx, i] = init_sol.walks[0][i][(var_name, var_idx)]
            elif var_name == "c":
                model.c[var_idx, i] = init_sol.walks[0][i][(var_name, var_idx)]
            elif var_name == "d":
                model.d[var_idx, i] = init_sol.walks[0][i][(var_name, var_idx)]
            elif var_name == "y1":
                model.y1[var_idx, i] = init_sol.walks[0][i][(var_name, var_idx)]
            elif var_name == "y2":
                model.y2[var_idx, i] = init_sol.walks[0][i][(var_name, var_idx)]
    for i in range(len(init_sol.walks[1])):
        i_ = i + len(init_sol.walks[0])
        model.z[i_] = 1
        model.w[i_] = init_sol.walk_weights[1][i]
        for v, vi in init_sol.walks[1][i].keys():
            if v == "x":
                model.x[vi, i_] = init_sol.walks[1][i][(v, vi)]
            elif v == "c":
                model.c[vi, i_] = init_sol.walks[1][i][(v, vi)]
            elif v == "d":
                model.d[vi, i_] = init_sol.walks[1][i][(v, vi)]
            elif v == "y1":
                model.y1[vi, i_] = init_sol.walks[1][i][(v, vi)]
            elif v == "y2":
                model.y2[vi, i_] = init_sol.walks[1][i][(v, vi)]
    for i in range(len(init_sol.satisfied_pc[0])):
        for pi in init_sol.satisfied_pc[0][i]:
            model.r[pi, i] = 1
            model.R[pi] = 1
    for i in range(len(init_sol.satisfied_pc[1])):
        i_ = i + len(init_sol.satisfied_pc[0])
        for pi in init_sol.satisfied_pc[1][i]:
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
    solver_options: Optional[datatypes.SolverOptions] = None,
    postprocess: int = 0,
) -> datatypes.CycleSolution:
    """Greedy cycle decomposition by maximizing the total length-weighted CN.

    The basis of this model can essentially be considered a case of the standard
    model used in `minimize_cycles` with fixed k = 1. The objective function of
    the solver maximizes total CN rather than minimizing cycles, producing a
    single cycle at a time. We repeatedly call this model/solver until we reach
    one of the stopping conditions (involving total CN, satisfied path
    constraint, and cycle weight respectively).

    Returns:
        A solution container parsed from the MIQCP solver output with the
        following attributes:
            (1) Pyomo status of optimization model
            (2) Total length weighted CN in resulting cycles/paths
            (3) Total num subpath constraints satisfied by resulting
                cycles/paths
            (4) List of cycles, each as a dict which maps an edge to its
                multiplicity in the cycle
            (5) List of the corresponding CN of the above cycles
            (6) Subpath constraints (indices) satisfied by each cycle
            (7) List of paths, each as a dict which maps an edge to its
                multiplicity in the path
            (8) List of the corresponding CN of the above paths
            (9) Subpath constraints (indices) satisfied by each path
    """

    if not solver_options:
        solver_options = datatypes.SolverOptions()
    solver_options.time_limit_s = max(
        solver_options.time_limit_s, bp_graph.num_disc_edges * 300
    )  # each breakpoint edge is assigned >= 5 minutes)

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
            f"Iteration {cycle_id + 1} with remaining CN = {remaining_weights} and unsatisfied constraints = {num_unsatisfied_pc}/{len(pc_list)}."
        )
        logger.debug(f"Multiplication factor for subpath constraints = {pp}.")

        model_name = f"amplicon_{amplicon_id}_cycle_decomposition_greedy_{cycle_id + 1}_{alpha=}"
        model_filepath = f"{solver_options.output_dir}/{solver_options.model_prefix}_{model_name}"

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
            f"{model_filepath}.lp", io_options={"symbolic_solver_labels": True}
        )
        model.write(f"{model_filepath}_ampl.nl", format="nl")
        logger.debug(f"Completed model setup, wrote to {model_filepath}.lp.")

        solver = cycle_utils.get_solver(solver_options)
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
            logger.debug("Greedy cycle decomposition is infeasible, stopping.")
            break
        cycle_id += 1

        # Check if solver produced a cycle or path
        if curr_sol.num_cycles > 0:
            curr_walk_weight = curr_sol.walk_weights.cycles[0]
            full_solution.walks.cycles.append(curr_sol.walks.cycles[0])
            full_solution.walk_weights.cycles.append(curr_walk_weight)
            full_solution.satisfied_pc.cycles.append(
                curr_sol.satisfied_pc.cycles[0]
            )
        elif curr_sol.num_paths > 0:
            curr_walk_weight = curr_sol.walk_weights.paths[0]
            full_solution.walks.paths.append(curr_sol.walks.paths[0])
            full_solution.walk_weights.paths.append(curr_walk_weight)
            full_solution.satisfied_pc.paths.append(
                curr_sol.satisfied_pc.paths[0]
            )
        full_solution.satisfied_pc_set |= curr_sol.satisfied_pc_set
        for i in curr_sol.satisfied_pc_set:
            is_pc_unsatisfied[i] = False
        logger.debug(f"{is_pc_unsatisfied=}")

        # Update greedy stop conditions based on latest solution
        next_w = curr_walk_weight
        num_unsatisfied_pc = len(pc_list) - len(full_solution.satisfied_pc_set)
        for pi in curr_sol.satisfied_pc_set:
            is_pc_unsatisfied[pi] = False
        remaining_weights -= curr_sol.total_weights_included
        if curr_sol.total_weights_included < cn_tol * total_weights:
            logger.debug(
                f"Proportion of length-weighted CN less than {cn_tol=}, "
                "iteration terminated."
            )
            break
    return full_solution


def cycle_decomposition(
    bb: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
    solver_options: datatypes.SolverOptions,
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    postprocess: int = 0,
    *,
    output_all_path_constraints: bool = False,
) -> None:
    """Caller for cycle decomposition functions"""
    was_amplicon_solved: Dict[int, bool] = defaultdict(bool)  # default false
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
                lp_solution = maximize_weights_greedy(
                    amplicon_id=amplicon_idx + 1,
                    bp_graph=bb.lr_graph[amplicon_idx],
                    total_weights=total_weights,
                    node_order=node_order,
                    pc_list=bb.longest_path_constraints[amplicon_idx][0],
                    cycle_id=0,
                    alpha=alpha,
                    p_total_weight=p_total_weight,
                    resolution=resolution,
                    cn_tol=0.005,
                    p_subpaths=0.9,
                    solver_options=solver_options,
                    postprocess=postprocess,
                )
                if not lp_solution:
                    logger.info("Greedy cycle decomposition failed.")
                else:
                    logger.info("Completed greedy cycle decomposition.")
                    logger.info(
                        f"Num cycles = {lp_solution.num_cycles}; num paths = {lp_solution.num_paths}."
                    )
                    logger.info(
                        f"Total length weighted CN = {lp_solution.total_weights_included}/{total_weights}."
                    )
                    logger.info(
                        f"Total num subpath constraints satisfied = {lp_solution.num_pc_satisfied}/{len(bb.longest_path_constraints[amplicon_idx][0])}."
                    )
                if postprocess == 1:
                    lp_solution = postprocess_solution(
                        amplicon_idx,
                        bb,
                        total_weights,
                        node_order,
                        p_total_weight,
                        resolution,
                        lp_solution,
                        solver_options,
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
            lp_solution = minimize_cycles(
                amplicon_idx + 1,
                bb.lr_graph[amplicon_idx],
                k,
                total_weights,
                node_order,
                bb.longest_path_constraints[amplicon_idx][0],
                p_total_weight,
                0.9,
                solver_options,
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
            lp_solution = maximize_weights_greedy(
                amplicon_id=amplicon_idx + 1,
                bp_graph=bb.lr_graph[amplicon_idx],
                total_weights=total_weights,
                node_order=node_order,
                pc_list=bb.longest_path_constraints[amplicon_idx][0],
                cycle_id=0,
                alpha=alpha,
                p_total_weight=p_total_weight,
                resolution=resolution,
                cn_tol=0.005,
                p_subpaths=0.9,
                solver_options=solver_options,
                postprocess=postprocess,
            )
            logger.info("Completed greedy cycle decomposition.")
            logger.info(
                f"Num cycles = {lp_solution.num_cycles}; num paths = {lp_solution.num_paths}."
            )
            logger.info(
                f"Total length weighted CN = {lp_solution.total_weights_included}/{total_weights}."
            )
            logger.info(
                f"Total num subpath constraints satisfied = {lp_solution.num_pc_satisfied}/{len(bb.longest_path_constraints[amplicon_idx][0])}."
            )
            if postprocess == 1:
                lp_solution = postprocess_solution(
                    amplicon_idx,
                    bb,
                    total_weights,
                    node_order,
                    p_total_weight,
                    resolution,
                    lp_solution,
                    solver_options,
                )
            bb.walks_by_amplicon[amplicon_idx] = lp_solution.walks
            bb.walk_weights_by_amplicon[amplicon_idx] = lp_solution.walk_weights
            bb.path_constraints_satisfied[amplicon_idx] = (
                lp_solution.satisfied_pc
            )
        else:
            output.output_amplicon_cycles(
                amplicon_idx,
                bb,
                solver_options.output_dir,
                output_all_path_constraints,
            )
            was_amplicon_solved[amplicon_idx] = True
    output.output_summary_amplicon_stats(
        was_amplicon_solved,
        bb,
        output_dir=solver_options.output_dir,
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
        solver_options=datatypes.SolverOptions(
            nthreads, time_limit_, output_dir, "pyomo", solver_to_use
        ),
        alpha=alpha_,
        postprocess=postprocess_,
        output_all_path_constraints=output_all_path_constraints,
    )
    logger.info("Completed cycle decomposition for all amplicons.")
    logger.info(
        f"Wrote cycles for all complicons to {output_dir}/amplicon*_cycles.txt."
    )


def postprocess_solution(
    amplicon_idx: int,
    bam_to_bps: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
    total_weights: float,
    node_order: Dict[tuple[Any, Any, Any], int],
    p_total_weight: float,
    resolution: float,
    lp_solution: datatypes.CycleSolution,
    solver_options: datatypes.SolverOptions,
) -> datatypes.CycleSolution:
    """Postprocess the solution generated by the greedy cycle decomposition.

    Args:
        lp_solution: solution container parsed from the MIQCP solver output

    Returns:
        CycleLPSolution - modifies the input solution object directly with the postprocessed values.
    """
    lp_solution = minimize_cycles_post(
        amplicon_idx + 1,
        bam_to_bps.lr_graph[amplicon_idx],
        total_weights,
        node_order,
        bam_to_bps.longest_path_constraints[amplicon_idx][0],
        datatypes.InitialSolution(
            lp_solution.walks,
            lp_solution.walk_weights,
            lp_solution.satisfied_pc,
        ),
        min(
            lp_solution.total_weights_included / total_weights * 0.9999,
            p_total_weight,
        ),
        resolution,
        solver_options,
    )

    logger.info("Completed postprocessing of the greedy solution.")
    logger.info(
        f"Num cycles = {len(lp_solution.walks[0])}; num paths = {len(lp_solution.walks[1])}."
    )
    logger.info(
        f"Total length weighted CN = {lp_solution.total_weights_included}/{total_weights}."
    )
    logger.info(
        f"Total num subpath constraints satisfied = {lp_solution.num_pc_satisfied}/{len(bam_to_bps.longest_path_constraints[amplicon_idx][0])}."
    )
    return lp_solution
