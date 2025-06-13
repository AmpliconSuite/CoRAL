"""Functions used for cycle decomposition"""

from __future__ import annotations

import logging
import pathlib
from collections import defaultdict
from typing import Dict, List, Optional, cast

import memray
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

import coral.output.utils
import coral.summary
from coral import core_utils, datatypes, models
from coral.breakpoint import breakpoint_graph, infer_breakpoint_graph
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.core_utils import profile_fn_with_call_counter
from coral.datatypes import OutputPCOptions
from coral.models import cycle_utils
from coral.models.concrete import initialize_post_processing_solver
from coral.models.path_constraints import longest_path_dict
from coral.output import cycle_output

logger = logging.getLogger(__name__)


def minimize_cycles(
    bp_graph: BreakpointGraph,
    k: int,
    total_weights: float,
    node_order: Dict[datatypes.Node, int],
    p_total_weight: float = 0.9,
    p_bp_cn: float = 0.9,
    solver_options: datatypes.SolverOptions | None = None,
) -> datatypes.CycleSolution:
    """Cycle decomposition by minimizing the number of cycles/paths.

    Standard model for cycle detection used when generated BreakpointGraph
    is tractable.

    Args:
        amplicon_id: integer, amplicon ID
        g: breakpoint graph (object)
        k: integer, maximum mumber of cycles/paths allowed in cycle decomposition
        total_weights: float, total length-weighted CN in breakpoint graph g
        node_order: dict maps each node in the input breakpoint graph to a distinct integer, indicating a total order of the nodes in g
        pc_list: list of subpath constraints to be satisfied, each as a dict
            that maps an edge to its multiplicity
                *** note that all subpath constraints in this list are
                required to be satisfied ***
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
        A solution container parsed from the MIQCP solver output with the
        following attributes:
            (1) Pyomo status of optimization model
            (2) Total length weighted CN in resulting cycles/paths
            (3) Total num subpath constraints satisfied by resulting cycles/paths
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
        f"Regular cycle decomposition with at most {k} cycles/paths allowed."
    )
    if not solver_options:
        solver_options = datatypes.SolverOptions()
    # solver_options.time_limit_s = max(
    #     solver_options.time_limit_s, bp_graph.num_disc_edges * 300
    # )  # each breakpoint edge is assigned >= 5 minutes)

    model_name = (
        f"amplicon_{bp_graph.amplicon_idx+1}_cycle_decomposition_k_{k}"
    )
    prefixed_name = f"{solver_options.output_prefix}_{solver_options.model_prefix}_{model_name}"
    if solver_options.output_prefix == "":
        prefixed_name = f"{solver_options.model_prefix}_{model_name}"
    model_filepath = f"{solver_options.output_dir}/models/{prefixed_name}"

    model = models.concrete.generate_model(
        bp_graph,
        k,
        total_weights,
        node_order,
        model_name=model_name,
        model_filepath=model_filepath,
    )

    model.write(
        f"{model_filepath}.lp", io_options={"symbolic_solver_labels": True}
    )
    model.write(f"{model_filepath}_ampl.nl", format="nl")
    logger.debug(f"Completed model setup, wrote to {model_filepath}.lp.")

    solver = cycle_utils.get_solver(solver_options)
    results = solver.solve(model, model_filepath)

    logger.debug(
        f"Completed optimization with status {results.solver.status}, condition {results.solver.termination_condition}."
    )

    return cycle_utils.parse_solver_output(
        results,
        model,
        bp_graph,
        k,
        bp_graph.longest_path_constraints,
        total_weights,
        model_type=datatypes.ModelType.DEFAULT,
    )


def minimize_cycles_post(
    amplicon_id: int,
    bp_graph: BreakpointGraph,
    total_weights: float,
    node_order: Dict[datatypes.Node, int],
    pc_list: list,
    init_sol: datatypes.InitialSolution,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    solver_options: datatypes.SolverOptions | None = None,
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
    # solver_options.time_limit_s = max(
    #     solver_options.time_limit_s, bp_graph.num_disc_edges * 300
    # )  # each breakpoint edge is assigned >= 5 minutes)

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
        f"amplicon_{amplicon_id+1}_cycle_decomposition_postprocessing_k_{k}"
    )
    model_filepath = f"{solver_options.output_dir}/models/{solver_options.model_prefix}_{model_name}"
    if solver_options.model_prefix == "":
        model_filepath = f"{solver_options.output_dir}/models/{model_name}"

    model = models.concrete.generate_model(
        bp_graph,
        k,
        total_weights,
        node_order,
        model_name=model_name,
        model_filepath=model_filepath,
        is_post=True,
        resolution=resolution,
        p_total_weight=p_total_weight,
        init_sol=init_sol,
    )

    initialize_post_processing_solver(model, init_sol)

    solver = cycle_utils.get_solver(solver_options)
    results: cycle_utils.PyomoResults = solver.solve(
        model, model_filepath=model_filepath
    )
    return cycle_utils.parse_solver_output(
        results,
        model,
        bp_graph,
        k,
        pc_list,
        total_weights,
        model_type=datatypes.ModelType.DEFAULT,
    )


def greedy_solve(
    bp_graph: BreakpointGraph,
    total_weights: float,
    node_order: Dict[datatypes.Node, int],
    solver_options: datatypes.SolverOptions,
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    cn_tol: float = 0.005,
    p_subpaths: float = 0.9,
    *,
    should_postprocess: bool = False,
) -> datatypes.CycleSolution:
    lp_solution = maximize_weights_greedy(
        bp_graph=bp_graph,
        total_weights=total_weights,
        node_order=node_order,
        alpha=alpha,
        p_total_weight=p_total_weight,
        resolution=resolution,
        cn_tol=cn_tol,
        p_subpaths=p_subpaths,
        solver_options=solver_options,
    )
    if should_postprocess:
        lp_solution = postprocess_solution(
            bp_graph.amplicon_idx,
            bp_graph,
            total_weights,
            node_order,
            p_total_weight,
            resolution,
            lp_solution,
            solver_options,
        )

    if lp_solution.num_walks == 0:
        logger.warning(
            "Greedy cycle decomposition produced no cycles or paths."
        )
    elif lp_solution.termination_condition == pyo.TerminationCondition.maxTimeLimit:
        logger.warning(f"Completed partial greedy cycle decomposition after {lp_solution.model_metadata.k} iterations.")
    else:
        logger.info("Completed full greedy cycle decomposition.")
    return lp_solution


def maximize_weights_greedy(
    bp_graph: BreakpointGraph,
    total_weights: float,
    node_order: Dict[datatypes.Node, int],
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    cn_tol: float = 0.005,
    p_subpaths: float = 0.9,
    solver_options: datatypes.SolverOptions | None = None,
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
    # solver_options.time_limit_s = max(
    #     solver_options.time_limit_s, bp_graph.num_disc_edges * 300
    # )  # each breakpoint edge is assigned >= 5 minutes)

    remaining_weights = total_weights

    pc_list = bp_graph.longest_path_constraints
    is_pc_unsatisfied = [True for i in range(len(pc_list))]
    num_unsatisfied_pc = sum(is_pc_unsatisfied)
    remaining_cn = datatypes.EdgeToCN.from_graph(bp_graph)
    next_w = resolution * 1.1

    full_solution = datatypes.CycleSolution(
        solver_status=pyo.SolverStatus.unknown,
        termination_condition=pyo.TerminationCondition.unknown,
    )

    logger.debug(
        f"Greedy cycle decomposition with length weighted CN = "
        f"{remaining_weights} & num subpath constraints = {num_unsatisfied_pc}."
    )

    cycle_id = 1
    while next_w >= resolution and (
        remaining_weights > (1.0 - p_total_weight) * total_weights
        or num_unsatisfied_pc > np.floor((1.0 - p_subpaths) * len(pc_list))
    ):
        pp = 1.0
        if alpha > 0 and num_unsatisfied_pc > 0:
            # multi - objective optimization parameter
            pp = alpha * remaining_weights / num_unsatisfied_pc
        logger.debug(
            f"Iteration {cycle_id} with remaining CN = {remaining_weights} "
            f"& unsatisfied constraints = {num_unsatisfied_pc}/{len(pc_list)}."
        )
        logger.debug(f"Multiplication factor for subpath constraints = {pp}.")

        model_name = (
            f"amplicon_{bp_graph.amplicon_idx+1}_cycle_decomposition"
            f"_greedy_cycle{cycle_id}_alpha_{alpha}"
        )
        prefixed_name = f"{solver_options.output_prefix}_{solver_options.model_prefix}_{model_name}"
        if solver_options.model_prefix == "":
            prefixed_name = f"{solver_options.model_prefix}_{model_name}"
        model_filepath = f"{solver_options.output_dir}/models/{prefixed_name}"

        model = models.concrete.generate_model(
            bp_graph,
            k=1,
            total_weights=total_weights,
            node_order=node_order,
            model_name=model_name,
            model_filepath=model_filepath,
            is_greedy=True,
            pp=pp,
            is_pc_unsatisfied=is_pc_unsatisfied,
            remaining_cn=remaining_cn,
        )

        solver = cycle_utils.get_solver(solver_options)
        results = solver.solve(model, model_filepath=model_filepath)
        curr_sol = cycle_utils.parse_solver_output(
            results,
            model,
            bp_graph,
            1,
            pc_list,
            total_weights,
            model_type=datatypes.ModelType.GREEDY,
            remaining_cn=remaining_cn,
            resolution=resolution,
            is_pc_unsatisfied=is_pc_unsatisfied,
            alpha=alpha,
        )


        full_solution.termination_condition = curr_sol.termination_condition
        full_solution.solver_status = curr_sol.solver_status
        if (
            curr_sol.termination_condition
            == pyo.TerminationCondition.infeasible
        ):
            logger.debug("Greedy cycle decomposition is infeasible, stopping.")
            cycle_id -= 1 # Don't count failed iterations as cycles
            break
        if (
            curr_sol.termination_condition
            == pyo.TerminationCondition.maxTimeLimit
            and curr_sol.solver_status == pyo.SolverStatus.aborted
        ):
            logger.warning(f"Greedy cycle decomposition timed out during iteration {cycle_id}, stopping.")
            cycle_id -= 1 # Don't count failed iterations as cycles
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
        next_w = curr_walk_weight  # type: ignore[possibly-undefined]
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
    full_solution.model_metadata = datatypes.ModelMetadata(
        model_type=datatypes.ModelType.GREEDY,
        k=cycle_id,
        alpha=alpha,
        total_weights=total_weights,
        resolution=resolution,
        num_path_constraints=len(pc_list),
    )
    return full_solution


def does_graph_require_greedy_solve(bp_graph: BreakpointGraph, k: int) -> bool:
    # We avoid using the cycle-minimization ILP model for large graphs based on
    # this heuristic.
    return (
        bp_graph.num_edges > 100
        or (
            3 * k
            + 3 * k * bp_graph.num_edges
            + 2 * k * bp_graph.num_nodes
            + k * len(bp_graph.longest_path_constraints)
        )
        >= 10000
    )


@core_utils.profile_fn_with_call_counter
def solve_single_graph(
    bp_graph: BreakpointGraph,
    k: int,
    total_weights: float,
    node_order: Dict[datatypes.Node, int],
    solver_options: datatypes.SolverOptions,
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    p_bp_cn: float = 0.9,
    resolution: float = 0.1,
    *,
    should_force_greedy: bool = False,
    should_postprocess: bool = False,
) -> datatypes.CycleSolution:
    used_greedy = False
    while k <= bp_graph.num_edges:
        # When problem is too large, we begin with the greedy optimization
        if should_force_greedy or does_graph_require_greedy_solve(bp_graph, k):
            lp_solution = greedy_solve(
                bp_graph=bp_graph,
                total_weights=total_weights,
                node_order=node_order,
                solver_options=solver_options,
                alpha=alpha,
                p_total_weight=p_total_weight,
                resolution=resolution,
                cn_tol=0.005,
                p_subpaths=0.9,
                should_postprocess=should_postprocess,
            )
            used_greedy = True
            break

        lp_solution = minimize_cycles(
            bp_graph=bp_graph,
            k=k,
            total_weights=total_weights,
            node_order=node_order,
            p_total_weight=p_total_weight,
            p_bp_cn=p_bp_cn,
            solver_options=solver_options,
        )
        if (
            lp_solution.termination_condition
            != pyo.TerminationCondition.infeasible
        ):
            logger.info(f"Completed cycle decomposition with k = {k}.")
            break

        logger.info(
            f"Cycle decomposition is infeasible, doubling k from {k} to {k * 2}."
        )
        k *= 2
        continue

    was_unsolved = False
    if lp_solution.termination_condition == pyo.TerminationCondition.infeasible:
        logger.info(
            f"Cycle decomposition is infeasible despite k > {bp_graph.num_edges}"
        )
        was_unsolved = True
    elif (
        lp_solution.termination_condition
        == pyo.TerminationCondition.maxTimeLimit
        and lp_solution.solver_status == pyo.SolverStatus.aborted
    ):
        logger.info(
            "Cycle decomposition timed out without a valid solution, "
            "attempting greedy cycle decomposition."
        )
        was_unsolved = True

    if was_unsolved and not used_greedy:
        lp_solution = greedy_solve(
            bp_graph=bp_graph,
            total_weights=total_weights,
            node_order=node_order,
            solver_options=solver_options,
            alpha=alpha,
            p_total_weight=p_total_weight,
            resolution=resolution,
        )
        if (
            lp_solution.termination_condition
            == pyo.TerminationCondition.maxTimeLimit
            and lp_solution.solver_status == pyo.SolverStatus.aborted
        ):
            logger.warning(
                "Greedy cycle decomposition timed out without a valid solution."
            )
            return lp_solution

    logger.info(
        f"Num cycles = {lp_solution.num_cycles}; "
        f"num paths = {lp_solution.num_paths}."
    )
    logger.info(
        f"Total length weighted CN = "
        f"{lp_solution.total_weights_included}/{total_weights}."
    )
    logger.info(
        f"Total num subpath constraints satisfied = "
        f"{lp_solution.num_pc_satisfied}/"
        f"{len(bp_graph.longest_path_constraints)}."
    )
    logger.info(f"Completed cycle decomposition with k = {k}.")
    return lp_solution


def cycle_decomposition_single_graph(
    bp_graph: BreakpointGraph,
    solver_options: datatypes.SolverOptions,
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    *,
    should_postprocess: bool = False,
    should_force_greedy: bool = False,
    pc_output_option: OutputPCOptions = OutputPCOptions.LONGEST,
    ignore_path_constraints: bool = False,
) -> None:
    logger.info(
        f"Begin cycle decomposition for amplicon {bp_graph.amplicon_idx + 1}."
    )

    total_weights = sum(len(sseg) * sseg.cn for sseg in bp_graph.sequence_edges)

    logger.info(f"Total CN weights = {total_weights}.")

    if ignore_path_constraints:
        bp_graph.longest_path_constraints = []

    longest_path_constraints = bp_graph.longest_path_constraints

    logger.info(
        f"Total num maximal subpath constraints = "
        f"{len(longest_path_constraints)}."
    )
    for pc in longest_path_constraints:
        logger.debug(
            f"Subpath constraint {pc.pc_idx} (support={pc.support}) = "
            f"{bp_graph.path_constraints[pc.pc_idx]}"
        )

    k = max(10, len(bp_graph.discordant_edges) // 2)
    logger.info(f"Initial num cycles/paths = {k}.")

    node_order = {
        node: node_ordered_idx
        for node_ordered_idx, node in enumerate(bp_graph.node_adjacencies)
    }

    if bp_graph.num_edges < k:
        k = bp_graph.num_edges
        logger.info(f"Reset num cycles/paths to {k}.")

    lp_solution = solve_single_graph(
        bp_graph=bp_graph,
        k=k,
        total_weights=total_weights,
        node_order=node_order,
        solver_options=solver_options,
        alpha=alpha,
        p_total_weight=p_total_weight,
        resolution=resolution,
        should_force_greedy=should_force_greedy,
        should_postprocess=should_postprocess,
    )
    bp_graph.model_metadata = lp_solution.model_metadata
    if (
        lp_solution.termination_condition
        == pyo.TerminationCondition.maxTimeLimit
    ):
        if (not lp_solution.walks.paths and not lp_solution.walks.cycles):
            logger.warning(
                "Cycle decomposition timed out without a valid solution."
            )
            bp_graph.solution_status = datatypes.SolutionStatus.FAILURE
        else:
            logger.warning(
                "Cycle decomposition timed out with a partial solution."
            )
            bp_graph.solution_status = datatypes.SolutionStatus.PARTIAL
    else:
        bp_graph.solution_status = datatypes.SolutionStatus.SUCCESS
    bp_graph.walks = lp_solution.walks
    bp_graph.walk_weights = lp_solution.walk_weights
    bp_graph.path_constraints_satisfied = lp_solution.satisfied_pc
    bp_graph.mip_gap = lp_solution.mip_gap
    bp_graph.upper_bound = lp_solution.upper_bound
    cycle_output.output_amplicon_walks(
        bp_graph,
        str(solver_options.output_dir) + "/" + solver_options.output_prefix,
        pc_output_option=pc_output_option,
    )
    coral.summary.output.output_summary_amplicon_stats_single(bp_graph)


def cycle_decomposition_all_graphs(
    bp_graphs: list[BreakpointGraph],
    solver_options: datatypes.SolverOptions,
    alpha: float = 0.01,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,
    *,
    should_postprocess: bool = False,
    should_force_greedy: bool = False,
    pc_output_option: OutputPCOptions = OutputPCOptions.LONGEST,
    ignore_path_constraints: bool = False,
) -> None:
    """Caller for cycle decomposition functions"""
    was_amplicon_solved: Dict[int, bool] = defaultdict(bool)  # default false
    for bp_graph in bp_graphs:
        amplicon_idx = bp_graph.amplicon_idx
        cycle_decomposition_single_graph(
            bp_graph,
            solver_options,
            alpha,
            p_total_weight,
            resolution,
            should_postprocess=should_postprocess,
            should_force_greedy=should_force_greedy,
            pc_output_option=pc_output_option,
            ignore_path_constraints=ignore_path_constraints,
        )

        # Verify that cycle decomposition produced a valid solution
        if bp_graph.walks.paths or bp_graph.walks.cycles:
            was_amplicon_solved[amplicon_idx] = True
        else:
            was_amplicon_solved[amplicon_idx] = False

        coral.summary.output.output_summary_amplicon_stats_all(
            was_amplicon_solved, bp_graphs
        )


def reconstruct_cycles(
    bp_graphs: list[BreakpointGraph],
    solver_options: datatypes.SolverOptions,
    cycle_decomp_mode: datatypes.CycleDecompOptions,
    cycle_decomp_alpha: float = 0.01,
    *,
    should_postprocess_greedy_sol: bool = False,
    pc_output_option: OutputPCOptions = OutputPCOptions.LONGEST,
    ignore_path_constraints: bool = False,
) -> None:
    if solver_options.output_prefix == "":
        logging.basicConfig(
            filename=f"{solver_options.output_dir}/cycle_decomp.log",
            filemode="w",
            level=logging.INFO,
            format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
        )
    else:
        logging.basicConfig(
            filename=f"{solver_options.output_dir}/{solver_options.output_prefix}_cycle_decomp.log",
            filemode="w",
            level=logging.INFO,
            format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
        )
    
    if cycle_decomp_mode == datatypes.CycleDecompOptions.MAX_WEIGHT:
        cycle_decomposition_all_graphs(
            bp_graphs,
            solver_options=solver_options,
            alpha=cycle_decomp_alpha,
            should_postprocess=should_postprocess_greedy_sol,
            should_force_greedy=True,
            pc_output_option=pc_output_option,
            ignore_path_constraints=ignore_path_constraints,
        )
    else:
        cycle_decomposition_all_graphs(
            bp_graphs,
            solver_options=solver_options,
            alpha=cycle_decomp_alpha,
            should_postprocess=should_postprocess_greedy_sol,
            should_force_greedy=False,
            pc_output_option=pc_output_option,
            ignore_path_constraints=ignore_path_constraints,
        )
    logger.info("Completed cycle decomposition for all amplicons.")
    if solver_options.output_prefix == "":
        logger.info(
            "Wrote cycles for all complicons to "
            f"{solver_options.output_dir}/amplicon*_cycles.txt."
        )
    else:
        logger.info(
            "Wrote cycles for all complicons to "
            f"{solver_options.output_dir}/{solver_options.output_prefix}_amplicon*_cycles.txt."
        )


def postprocess_solution(
    amplicon_idx: int,
    bp_graph: BreakpointGraph,
    total_weights: float,
    node_order: Dict[datatypes.Node, int],
    p_total_weight: float,
    resolution: float,
    lp_solution: datatypes.CycleSolution,
    solver_options: datatypes.SolverOptions,
) -> datatypes.CycleSolution:
    """Postprocess the solution generated by the greedy cycle decomposition.

    Args:
        lp_solution: solution container parsed from the MIQCP solver output

    Returns:
        CycleLPSolution - modifies the input solution object directly with the
        postprocessed values.
    """
    lp_solution = minimize_cycles_post(
        amplicon_idx + 1,
        bp_graph,
        total_weights,
        node_order,
        bp_graph.longest_path_constraints,
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
        f"Total num subpath constraints satisfied = "
        f"{lp_solution.num_pc_satisfied}/"
        f"{len(bp_graph.longest_path_constraints)}."
    )
    return lp_solution
