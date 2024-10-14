import logging
from typing import List, Optional, Tuple
from coral.breakpoint_graph import BreakpointGraph
import pyomo.environ as pyo
import pyomo.core
import pyomo.opt
import pyomo.util.infeasible

from coral.models.utils import EdgeToCN

logger = logging.getLogger(__name__)

def get_total_weight_constraint(model: pyo.ConcreteModel, k: int, bp_graph: BreakpointGraph, weight_threshold: float) -> pyo.Constraint:
    total_weight_expr = 0.0
    for i in range(k):
        for seqi in range(bp_graph.num_seq_edges):
            total_weight_expr += model.x[seqi, i] * model.w[i] * bp_graph.sequence_edges[seqi][-2]
    return pyo.Constraint(expr=total_weight_expr >= weight_threshold)

def get_eulerian_path_constraint(endnode_list, bp_graph: BreakpointGraph):
    def constrain_eulerian_path(model: pyo.Model, i: int) -> pyo.Expression:
        path_expr = 0.0
        for enodei in range(len(endnode_list)):
            edge_idx_offset = bp_graph.num_nonsrc_edges + 2 * bp_graph.num_src_edges + 2 * enodei
            path_expr += model.x[edge_idx_offset, i]  # (s, v)
            path_expr -= model.x[edge_idx_offset + 1, i]  # (v, t)
        for srci in range(bp_graph.num_src_edges):
            edge_idx_offset = bp_graph.num_nonsrc_edges + 2 * srci
            path_expr += model.x[edge_idx_offset, i]  # (s, v)
            path_expr -= model.x[edge_idx_offset + 1, i]  # (v, t)
        return path_expr == 0.0
    return constrain_eulerian_path

def get_eulerian_node_constraints(model: pyo.ConcreteModel, k: int, endnode_list, lseg, lc, ld, lsrc, bp_graph: BreakpointGraph) -> List[pyo.Constraint]:
    eulerian_node_constraints = []
    for node in bp_graph.nodes:
        if node in endnode_list:
            for i in range(k):
                edge_idx = bp_graph.num_nonsrc_edges + 2 * bp_graph.num_src_edges + 2 * endnode_list.index(node)
                eulerian_node_constraints.append(
                    model.x[edge_idx, i] + model.x[edge_idx + 1, i] == model.x[bp_graph.nodes[node][0][0], i]
                )
        else:
            for i in range(k):
                ec_expr = 0.0
                for seqi in bp_graph.nodes[node][0]:
                    ec_expr += model.x[seqi, i]
                for ci in bp_graph.nodes[node][1]:
                    ec_expr -= model.x[(lseg + ci), i]
                for di in bp_graph.nodes[node][2]:
                    ec_expr -= model.x[(lseg + lc + di), i]
                for srci in bp_graph.nodes[node][3]:
                    ec_expr -= model.x[(lseg + lc + ld + 2 * srci), i]  # connected to s
                    ec_expr -= model.x[(lseg + lc + ld + 2 * srci) + 1, +i]  # connected to t
                eulerian_node_constraints.append(ec_expr == 0.0)
    return eulerian_node_constraints

# Same base/post/greedy
def get_spanning_tree_constraints(model: pyo.ConcreteModel, k :int,  bp_graph: BreakpointGraph, endnode_list, node_order) -> List[pyo.Constraint]:
    constraints = []
    for i in range(k):
        t_expr_x = 0.0  # linear
        t_expr_y = 0.0  # linear
        t_expr_yd = 0.0  # quad
        for node in endnode_list:
            idx_offset = bp_graph.num_nonsrc_edges + 2 * bp_graph.num_src_edges + 2 * endnode_list.index(node)
            t_expr_x += model.x[idx_offset + 1, i]
            t_expr_y += model.y1[idx_offset + 1, i]
            t_expr_yd += model.y1[idx_offset + 1, i] * (
                model.dt[i] - model.d[node_order[node], i]
            )  # node -> t
            expr_x = 0.0  # linear
            expr_y = 0.0  # linear
            expr_xc = 0.0  # quad
            expr_yd = 0.0  # quad
            for seqi in bp_graph.nodes[node][0]:
                sseg = bp_graph.sequence_edges[seqi]
                node_ = (sseg[0], sseg[1], "-")
                if node_ == node:
                    node_ = (sseg[0], sseg[2], "+")
                expr_x += model.x[seqi, i]
                expr_xc += model.x[seqi, i] * model.c[node_order[node], i]
                if node_order[node_] <= node_order[node]:
                    expr_y += model.y1[seqi, i]
                    expr_yd += model.y1[seqi, i] * (model.d[node_order[node], i] - model.d[node_order[node_], i])
                else:
                    expr_y += model.y2[seqi, i]
                    expr_yd += model.y2[seqi, i] * (model.d[node_order[node], i] - model.d[node_order[node_], i])

            expr_x += model.x[idx_offset, i]  # from s
            expr_xc += (
                model.x[idx_offset, i] * model.c[node_order[node], i]
            )
            expr_x += model.x[idx_offset + 1, i]  # to t
            expr_xc += (
                model.x[idx_offset + 1, i]
                * model.c[node_order[node], i]
            )
            expr_y += model.y1[idx_offset, i]  # from s
            expr_yd += model.y1[idx_offset, i] * (
                model.d[node_order[node], i] - model.ds[i]
            )
            constraints.append(expr_x * (bp_graph.num_nodes + 2) >= model.d[node_order[node], i])
            constraints.append(expr_y <= 1.0)
            constraints.append(expr_y * bp_graph.num_edges * k + expr_xc >= expr_x)
            constraints.append(expr_yd * bp_graph.num_edges * k + expr_xc >= expr_x)
        
        for srci in range(bp_graph.num_src_edges):
            srce = bp_graph.source_edges[srci]
            idx_offset = bp_graph.num_nonsrc_edges + 2 * srci
            t_expr_x += model.x[idx_offset + 1, i]
            t_expr_y += model.y1[idx_offset + 1, i]
            t_expr_yd += model.y1[idx_offset + 1, i] * (
                model.dt[i] - model.d[node_order[(srce[3], srce[4], srce[5])], i]
            )
        constraints.append(t_expr_x * (bp_graph.num_nodes + 2) >= model.dt[i])
        constraints.append(t_expr_y <= 1.0)
        constraints.append(t_expr_y * bp_graph.num_edges * k >= t_expr_x)
        constraints.append(t_expr_yd >= t_expr_x)
    return constraints


# Same base + post
def get_spanning_tree_constraints__non_endnode(model: pyo.ConcreteModel, k: int,  bp_graph: BreakpointGraph, endnode_list, node_order) -> List[pyo.Constraint]:
    constraints = []
    for i in range(k):
        for node in bp_graph.nodes.keys():
            if node not in endnode_list:
                expr_x = 0.0
                expr_y = 0.0
                expr_xc = 0.0
                expr_yd = 0.0
                for seqi in bp_graph.nodes[node][0]:
                    sseg = bp_graph.sequence_edges[seqi]
                    node_ = (sseg[0], sseg[1], "-")
                    if node_ == node:
                        node_ = (sseg[0], sseg[2], "+")
                    expr_x += model.x[seqi, i]
                    expr_xc += model.x[seqi, i] * model.c[node_order[node], i]
                    if node_order[node_] <= node_order[node]:
                        expr_y += model.y1[seqi, i]
                        expr_yd += model.y1[seqi, i] * (model.d[node_order[node], i] - model.d[node_order[node_], i])
                    else:
                        expr_y += model.y2[seqi, i]
                        expr_yd += model.y2[seqi, i] * (model.d[node_order[node], i] - model.d[node_order[node_], i])
                for ci in bp_graph.nodes[node][1]:
                    cedge = bp_graph.concordant_edges[ci]
                    node_ = (cedge[0], cedge[1], cedge[2])
                    if node_ == node:
                        node_ = (cedge[3], cedge[4], cedge[5])
                    expr_x += model.x[(bp_graph.num_seq_edges + ci), i]
                    expr_xc += model.x[(bp_graph.num_seq_edges + ci), i] * model.c[node_order[node], i]
                    if node_order[node_] <= node_order[node]:
                        expr_y += model.y1[(bp_graph.num_seq_edges + ci), i]
                        expr_yd += model.y1[(bp_graph.num_seq_edges + ci), i] * (
                            model.d[node_order[node], i] - model.d[node_order[node_], i]
                        )
                    else:
                        expr_y += model.y2[(bp_graph.num_seq_edges + ci), i]
                        expr_yd += model.y2[(bp_graph.num_seq_edges + ci), i] * (
                            model.d[node_order[node], i] - model.d[node_order[node_], i]
                        )
                for di in bp_graph.nodes[node][2]:
                    dedge = bp_graph.discordant_edges[di]
                    node_ = (dedge[0], dedge[1], dedge[2])
                    if node_ == node:
                        node_ = (dedge[3], dedge[4], dedge[5])
                    idx_offset = bp_graph.num_seq_edges + bp_graph.num_conc_edges + di
                    expr_x += model.x[idx_offset, i]
                    expr_xc += model.x[idx_offset, i] * model.c[node_order[node], i]
                    if node_order[node_] <= node_order[node]:
                        expr_y += model.y1[idx_offset, i]
                        expr_yd += model.y1[idx_offset, i] * (
                            model.d[node_order[node], i] - model.d[node_order[node_], i]
                        )
                    else:
                        expr_y += model.y2[idx_offset, i]
                        expr_yd += model.y2[idx_offset, i] * (
                            model.d[node_order[node], i] - model.d[node_order[node_], i]
                        )
                for srci in bp_graph.nodes[node][3]:
                    idx_offset = bp_graph.num_seq_edges + bp_graph.num_conc_edges + bp_graph.num_disc_edges + srci * 2
                    expr_x += model.x[idx_offset, i]
                    expr_x += model.x[idx_offset + 1, i]
                    expr_xc += model.x[idx_offset, i] * model.c[node_order[node], i]
                    expr_xc += model.x[idx_offset + 1, i] * model.c[node_order[node], i]
                    expr_y += model.y1[idx_offset, i]
                    expr_yd += model.y1[idx_offset, i] * (model.d[node_order[node], i] - model.ds[i])
                
                constraints.append(expr_x * (bp_graph.num_nodes + 2) >= model.d[node_order[node], i])
                constraints.append(expr_y <= 1.0)
                constraints.append(expr_y * bp_graph.num_edges * k + expr_xc >= expr_x)
                constraints.append(expr_yd * bp_graph.num_edges * k + expr_xc >= expr_x)
    return constraints

# Consistent across `minimize_cycles` and `minimize_cycles_post`
def get_shared_subpath_edge_constraints(model: pyo.ConcreteModel, k: int, pc_list: List, bp_graph: BreakpointGraph) -> List[pyo.Constraint]:
    subpath_constraints = []
    for pi, path_constraint_ in enumerate(pc_list):
        for edge in path_constraint_:
            for i in range(k):
                if edge[0] == "s":
                    subpath_constraints.append(model.x[edge[1], i] >= model.r[pi, i] * path_constraint_[edge])
                elif edge[0] == "c":
                    subpath_constraints.append(
                        model.x[(bp_graph.num_seq_edges + edge[1]), i] >= model.r[pi, i] * path_constraint_[edge]
                    )
                else:
                    subpath_constraints.append(
                        model.x[(bp_graph.num_seq_edges + bp_graph.num_conc_edges + edge[1]), i] >= model.r[pi, i] * path_constraint_[edge]
                    )

    return subpath_constraints

def get_addtl_subpath_edge_constraints(model: pyo.ConcreteModel, k: int, pc_list: List, is_post: bool = False, p_path_constraints: Optional[float] = None) -> List[pyo.Constraint]:
    if not is_post:
        return  [pyo.Constraint(model.pc_idx, rule=lambda model, pi: sum(model.r[pi, i] for i in range(k)) >= 1.0)]
    assert p_path_constraints, f"Need to pass p_path_constraints when generating post-processing model"
    constraints = []
    sum_R = 0.0
    for pi in range(len(pc_list)):
        sum_R += model.R[pi]
        sum_r = 0.0
        for i in range(k):
            sum_r += model.r[i]
            constraints.append(pyo.Constraint(model.R[pi] >= model.r[i]))
        constraints.append(pyo.Constraint(sum_r >= model.R[pi]))
    constraints.append(pyo.Constraint(sum_R >= p_path_constraints * len(pc_list)))
    return constraints
        
def set_copy_number_constraints(model: pyo.ConcreteModel, k: int, bp_graph: BreakpointGraph, edge_to_cn: EdgeToCN, 
                                is_post: bool = False, p_bp_cn: float = 0.9, is_greedy: bool = False) -> None:
    # CN constraints (all quadratic)
    # Use EdgeToCN for constraints rather than directly accessing BreakpointGraph in order to support greedy model
    model.ConstraintCNSource = pyo.Constraint(
        model.src_edge_idx,
        rule=lambda model, srci: sum(
            model.w[i] * (model.x[(bp_graph.num_nonsrc_edges + 2 * srci), i] + model.x[(bp_graph.num_nonsrc_edges + 2 * srci) + 1, i])
            for i in range(k)
        )
        <= edge_to_cn.source[srci],
    )
    model.ConstraintCNSequence = pyo.Constraint(
        model.seq_edge_idx,
        rule=lambda model, seqi: sum(model.w[i] * model.x[seqi, i] for i in range(k))
        <= edge_to_cn.sequence[seqi],
    )
    model.ConstraintCNConcordant = pyo.Constraint(
        model.conc_edge_idx,
        rule=lambda model, ci: sum(model.w[i] * model.x[(bp_graph.num_seq_edges + ci), i] for i in range(k))
        <= edge_to_cn.concordant[ci],
    )
    model.ConstraintCNDiscordant = pyo.Constraint(
        model.disc_edge_idx,
        rule=lambda model, di: sum(model.w[i] * model.x[(bp_graph.num_seq_edges + bp_graph.num_conc_edges + di), i] for i in range(k))
        <= edge_to_cn.discordant[di],
    )
    if not is_post and not is_greedy:
        # Base minimize cycles model has an additional constraint
        model.ConstraintCNDiscordantNonPost = pyo.Constraint(
            model.disc_edge_idx,
            rule=lambda model, di: sum(model.w[i] * model.x[(bp_graph.num_seq_edges + bp_graph.num_conc_edges + di), i] for i in range(k))
            >= p_bp_cn * edge_to_cn.discordant[di], # type: ignore[operator]
        )
    if is_greedy:
        model.ConstraintDiscordantGreedy = pyo.ConstraintList()
        for di in range(bp_graph.num_disc_edges):
            logger.debug(f"Ignored discordant edge {di} in greedy model")
            model.ConstraintDiscordantGreedy.add(model.x[(bp_graph.num_seq_edges + bp_graph.num_conc_edges + di)] == 0.0)


def get_bp_occurrence_constraint(model: pyo.ConcreteModel, bp_graph: BreakpointGraph) -> pyo.Constraint:
    # Occurrence of breakpoints in each cycle/path
    discordant_multiplicities = bp_graph.infer_discordant_edge_multiplicities()
    def constraint_bp_occurence(model: pyo.Model, i: int, di: int) -> pyo.Expression:
        return model.x[(bp_graph.num_seq_edges + bp_graph.num_conc_edges + di), i] <= discordant_multiplicities[di]
    return  pyo.Constraint(model.k, model.disc_edge_idx, rule=constraint_bp_occurence)


# Constraint identical for base + post
def get_cycle_weight_constraint(model: pyo.Model, bp_graph: BreakpointGraph, endnode_list) -> pyo.Constraint:
    # Relationship between c and x
    def constrain_c_x(model: pyo.Model, i: int) -> pyo.Expression:
        cycle_expr = 0.0
        cycle_expr += sum(model.c[ni, i] for ni in range(bp_graph.num_nodes))
        cycle_expr += sum(model.x[(bp_graph.num_nonsrc_edges + 2 * bp_graph.num_src_edges + 2 * enodei), i] for enodei in range(len(endnode_list)))   # (s, v)
        cycle_expr += sum(model.x[(bp_graph.num_nonsrc_edges + 2 * srci), i] for srci in range(bp_graph.num_src_edges))
        return cycle_expr <= 1.0

    return pyo.Constraint(model.k, rule=constrain_c_x)


# Constraint identical for base + post
def get_single_bp_edge_constraint(model: pyo.Model, bp_graph: BreakpointGraph, node_order) -> pyo.Constraint:
    # There must be a concordant/discordant edge occuring one time
    def constrain_singular_bp_edge(model: pyo.Model, i: int) -> pyo.Expression:
        expr_xc = 0.0
        if not bp_graph.nodes:
            return pyo.Constraint.Skip  #  need to skip if no nodes to avoid trivial constraint error
        for node in bp_graph.nodes:
            for ci in set(bp_graph.nodes[node][1]):
                expr_xc += model.c[node_order[node], i] * model.x[(bp_graph.num_seq_edges + ci), i]
            for di in set(bp_graph.nodes[node][2]):
                expr_xc += model.c[node_order[node], i] * model.x[(bp_graph.num_seq_edges + bp_graph.num_conc_edges + di), i]
        # Skip trivial constraints when all components have 0 coeff
        return expr_xc <= 1.0 if expr_xc else pyo.Constraint.Skip
    return pyo.Constraint(model.k, rule=constrain_singular_bp_edge)

# Constraint identical for base + post
def set_cycle_tree_constraints(model: pyo.Model, bp_graph: BreakpointGraph) -> None:
    # Relationship between c and d
    model.ConstraintCD1 = pyo.Constraint(model.k, model.node_idx, rule=lambda model, i, ni: model.d[ni, i] >= model.c[ni, i])
    model.ConstraintCD2 = pyo.Constraint(model.k, rule=lambda model, i: (sum(model.c[ni, i] for ni in range(bp_graph.num_nodes)) + model.ds[i]) <= 1.0)
    
