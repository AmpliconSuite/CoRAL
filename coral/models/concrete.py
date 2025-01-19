from __future__ import annotations

import logging
from typing import List

import pyomo.environ as pyo

from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import EdgeToCN
from coral.models import constraints

logger = logging.getLogger(__name__)


def get_minimize_objective(
    model: pyo.Model,
    bp_graph: BreakpointGraph,
    k: int,
    total_weights: float,
    pc_list: List,
    is_post: bool,
) -> pyo.Objective:
    obj_value = 1.0 if not is_post else 2.0
    for i in range(k):
        obj_value += model.z[i]
        for seqi in range(bp_graph.num_seq_edges):
            obj_value -= (
                model.x[seqi, i]
                * model.w[i]
                * bp_graph.sequence_edges[seqi].gap
                / total_weights
            )
    if is_post:
        obj_value -= sum(
            model.R[pi] / len(pc_list) for pi in range(len(pc_list))
        )
    return pyo.Objective(sense=pyo.minimize, expr=obj_value)


def get_greedy_objective(
    model: pyo.Model,
    bp_graph: BreakpointGraph,
    pc_list: List,
    is_pc_unsatisfied: List[bool],
    pp: float,
) -> pyo.Objective:
    obj_value = 0.0
    for seqi in range(bp_graph.num_seq_edges):
        obj_value += (
            model.x[seqi, 0] * model.w[0] * bp_graph.sequence_edges[seqi].gap
        )
    for pi in range(len(pc_list)):
        if is_pc_unsatisfied[pi]:
            obj_value += model.r[pi, 0] * max(pp, 1.0)
    return pyo.Objective(sense=pyo.maximize, expr=obj_value)


class CycleLPModel(pyo.ConcreteModel):
    # TODO: add all constraints as decorators

    def __init__(
        self,
        name: str,
        bp_graph: BreakpointGraph,
        k: int,
        pc_list: List,
        is_post: bool = False,
    ) -> None:
        super().__init__(name=name)

        # region Indices
        self.k = pyo.RangeSet(0, k - 1)
        self.edge_idx = pyo.RangeSet(0, bp_graph.num_edges - 1)
        self.seq_edge_idx = pyo.RangeSet(0, len(bp_graph.sequence_edges) - 1)
        self.conc_edge_idx = pyo.RangeSet(0, len(bp_graph.concordant_edges) - 1)
        self.disc_edge_idx = pyo.RangeSet(0, len(bp_graph.discordant_edges) - 1)
        self.src_edge_idx = pyo.RangeSet(0, len(bp_graph.source_edges) - 1)
        self.node_idx = pyo.RangeSet(0, len(bp_graph.node_adjacencies) - 1)
        # endregion

        # region Variables
        # z[i]: indicating whether cycle or path i exists
        self.z = pyo.Var(self.k, domain=pyo.Binary)
        # w[i]: the weight of cycle or path i, continuous variable
        self.w = pyo.Var(
            self.k, domain=pyo.NonNegativeReals, bounds=(0.0, bp_graph.max_cn)
        )
        # x: the number of times an edge occur in cycle or path i
        self.x = pyo.Var(
            self.edge_idx, self.k, domain=pyo.Integers, bounds=(0.0, 10.0)
        )

        # Subpath constraints
        self.pc_idx = pyo.RangeSet(0, len(pc_list) - 1)
        self.r = pyo.Var(self.pc_idx, self.k, within=pyo.Binary)
        self.is_post = is_post
        if is_post:
            self.R = pyo.Var(self.pc_idx, within=pyo.Binary)

        # c: decomposition i is a cycle, and start at particular node
        self.c = pyo.Var(self.node_idx, self.k, within=pyo.Binary)

        # d: BFS/spanning tree order of the nodes in decomposition i (A4.15-A4.17)
        self.d = pyo.Var(
            self.node_idx,
            self.k,
            within=pyo.NonNegativeIntegers,
            bounds=(0.0, bp_graph.num_nodes + 2),
        )
        self.ds = pyo.Var(
            self.k,
            within=pyo.NonNegativeIntegers,
            bounds=(0.0, bp_graph.num_nodes + 2),
        )
        self.dt = pyo.Var(
            self.k,
            within=pyo.NonNegativeIntegers,
            bounds=(0.0, bp_graph.num_nodes + 2),
        )

        # y: spanning tree indicator (directed)
        self.y1 = pyo.Var(
            self.edge_idx, self.k, domain=pyo.Binary
        )  # Positive direction
        self.y2 = pyo.Var(
            self.edge_idx, self.k, domain=pyo.Binary
        )  # Negative direction
        # endregion


def get_model(
    bp_graph: BreakpointGraph,
    k: int,
    total_weights: float,
    node_order,
    pc_list,
    model_name: str,
    is_post: bool = False,
    is_greedy: bool = False,
    p_total_weight: float = 0.9,
    resolution: float = 0.1,  # Used only in post
    p_path_constraints: float | None = None,  # Used only in post
    p_bp_cn: float = 0.9,  # Used only in base
    pp: float = 1.0,  # Used only in greedy
    is_pc_unsatisfied: list[bool] | None = None,  # Used only in greedy
    remaining_cn: EdgeToCN | None = None,  # Used only in greedy
) -> pyo.ConcreteModel:
    logger.debug(f"Num nodes to be used in QP = {bp_graph.num_nodes}")
    logger.debug(f"Num edges to be used in QP = {bp_graph.num_edges}")

    model = CycleLPModel(
        name=f"{model_name}_{k}",
        bp_graph=bp_graph,
        k=k,
        pc_list=pc_list,
        is_post=is_post,
    )

    # Relationship between w[i] and z[i]
    # Below constraint is shared by `minimize_cycles` and `minimize_cycles_post`
    model.ConstraintWZ = pyo.Constraint(
        model.k,
        rule=lambda model, i: model.w[i] <= model.z[i] * bp_graph.max_cn,
    )

    if is_post or is_greedy:
        model.ConstraintWZResolution = pyo.Constraint(
            model.k, rule=lambda model, i: model.w[i] >= model.z[i] * resolution
        )

    if not is_greedy:
        model.ObjectiveMinCycles = get_minimize_objective(
            model, bp_graph, k, total_weights, pc_list, is_post
        )
    else:
        model.ObjectiveGreedy = get_greedy_objective(
            model, bp_graph, pc_list, is_pc_unsatisfied or [], pp
        )

    # Must include at least 0.9 * total CN weights (bilinear constraint)
    if not is_greedy:
        model.ConstraintTotalWeights = constraints.get_total_weight_constraint(
            model, k, bp_graph, p_total_weight * total_weights
        )

    model.ConstraintEulerianPath = pyo.Constraint(
        model.k, rule=constraints.get_eulerian_path_constraint(bp_graph)
    )
    model.ConstraintEulerianNodes = pyo.ConstraintList()
    for constraint in constraints.get_eulerian_node_constraints(
        model, k, bp_graph
    ):
        model.ConstraintEulerianNodes.add(constraint)

    constraints.set_copy_number_constraints(
        model,
        k,
        bp_graph,
        remaining_cn or EdgeToCN.from_graph(bp_graph),
        p_bp_cn=p_bp_cn,
        is_post=is_post,
        is_greedy=is_greedy,
    )

    model.ConstraintBPOccurrence = constraints.get_bp_occurrence_constraint(
        model, bp_graph
    )
    model.ConstraintCycleWeight = constraints.get_cycle_weight_constraint(
        model, bp_graph
    )
    model.ConstraintSingularBPEdge = constraints.get_single_bp_edge_constraint(
        model, bp_graph, node_order
    )

    constraints.set_cycle_tree_constraints(model, bp_graph)

    # Relationship between y and z:
    model.ConstraintY1Z = pyo.Constraint(
        model.k,
        model.edge_idx,
        rule=lambda model, i, j: model.y1[j, i] <= model.z[i],
    )
    model.ConstraintY2Z = pyo.Constraint(
        model.k,
        model.edge_idx,
        rule=lambda model, i, j: model.y2[j, i] <= model.z[i],
    )

    # Relationship between x, y and d
    model.ConstraintXY = pyo.Constraint(
        model.k,
        model.edge_idx,
        rule=lambda model, i, j: model.y1[j, i] + model.y2[j, i]
        <= model.x[j, i],
    )
    model.ConstraintXY1D = pyo.ConstraintList()
    model.ConstraintXY2D = pyo.ConstraintList()
    for i in range(k):
        for di in range(bp_graph.num_disc_edges):
            dedge = bp_graph.discordant_edges[di]
            if dedge.is_self_loop:  # exclude self loops
                model.ConstraintXY1D.add(
                    model.y1[
                        (bp_graph.num_seq_edges + bp_graph.num_conc_edges + di),
                        i,
                    ]
                    == 0
                )
                model.ConstraintXY2D.add(
                    model.y2[
                        (bp_graph.num_seq_edges + bp_graph.num_conc_edges + di),
                        i,
                    ]
                    == 0
                )

    model.ConstraintsXYD = pyo.ConstraintList()
    for constraint in constraints.get_spanning_tree_constraints(
        model, k, bp_graph, node_order
    ):
        model.ConstraintsXYD.add(constraint)
    for constraint in constraints.get_spanning_tree_constraints__non_endnode(
        model, k, bp_graph, node_order
    ):
        model.ConstraintsXYD.add(constraint)

    model.ConstraintSubpathEdges = pyo.ConstraintList()
    for constraint in constraints.get_shared_subpath_edge_constraints(
        model, k, pc_list, bp_graph
    ):
        model.ConstraintSubpathEdges.add(constraint)

    if pc_list and not is_greedy:
        model.ConstraintSubpathEdgesAddtl = pyo.ConstraintList()
        for constraint in constraints.get_addtl_subpath_edge_constraints(
            model, k, pc_list, is_post, p_path_constraints
        ):
            model.ConstraintSubpathEdgesAddtl.add(constraint)

    return model
