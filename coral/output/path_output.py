from __future__ import annotations

import logging
import random

from coral import constants
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    ConcordantEdge,
    DirectedEdge,
    DirectedWalk,
    DiscordantEdge,
    EdgeId,
    EdgeType,
    Node,
    OptimizationWalk,
    PathMetric,
    Strand,
    Walk,
)

logger = logging.getLogger(__name__)


def eulerian_cycle_t(
    g: BreakpointGraph,
    edge_counts_next_cycle: dict[EdgeId, int],
    path_constraints_next_cycle: list[Walk],
    path_constraints_support: list[int],
) -> DirectedWalk:
    """Return an eulerian traversal of a cycle, represented by a list of
     directed edges.

    g: breakpoint graph (object)
    edges_next_cycle: subgraph induced by the cycle, as a dict that maps an edge
        to its multiplicity
    path_constraints_next_cycle: list of subpath constraints to be satisfied,
        each as a list of alternating nodes and edges
        ***
        Note: the traversal may not satisfy all subpath constraints in case not
        all subpath constraints are satisfied, return the eulerian traversal
        satisfying the maximum number of subpath constraints
        ***
    path_constraints_support: num long reads supporting each subpath constraint
    """
    lseg = len(g.sequence_edges)

    # A cycle is edge - node list starting and ending with the same edge
    eulerian_cycle: Walk = []

    # Since Eulerian, there could be subcycles in the middle of a cycle
    eulerian_subcycle: DirectedWalk = []  # Cycle in AA cycle format
    best_cycle: DirectedWalk = []  # Cycle in AA cycle format
    valid = 0
    num_trials = 0
    num_pc = len(path_constraints_next_cycle)
    unsatisfied_path_metric = PathMetric(
        path_idxs=list(range(num_pc)),
        path_length=100 * num_pc,
        path_support=100 * max(path_constraints_support + [0]),
    )
    while valid <= 0 and num_trials < 1000:
        valid = 1
        num_trials += 1
        eulerian_cycle = []
        eulerian_subcycle = []
        edges_cur = edge_counts_next_cycle.copy()

        # Start with the edge with smallest index and on the positive strand
        last_seq_edge_idx = lseg
        for edge in edges_cur:
            if edge[0] == "e":
                last_seq_edge_idx = min(last_seq_edge_idx, edge[1])
        last_edge_dir = "+"
        eulerian_cycle.append(EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx))
        eulerian_subcycle.append(
            DirectedEdge(last_seq_edge_idx + 1, Strand.FORWARD)
        )
        while len(edges_cur) > 0:
            seq_edge = g.sequence_edges[last_seq_edge_idx]
            node = Node(seq_edge.chr, seq_edge.end, Strand.FORWARD)
            if last_edge_dir == "-":
                node = Node(seq_edge.chr, seq_edge.start, Strand.REVERSE)
            eulerian_cycle.append(node)
            next_bp_edges = []  # Since cycle, only consider discordant edges and concordant edges
            for ci in g.node_adjacencies[node].concordant:
                next_bp_edges.append(EdgeId(EdgeType.CONCORDANT, ci))
            for di in g.node_adjacencies[node].discordant:
                next_bp_edges.append(EdgeId(EdgeType.DISCORDANT, di))
            del_list = [
                i
                for i in range(len(next_bp_edges))
                if next_bp_edges[i] not in edges_cur
            ]
            for i in del_list[::-1]:
                del next_bp_edges[i]
            if len(next_bp_edges) == 0:
                valid = 0
                break
            bp_edge: ConcordantEdge | DiscordantEdge
            if len(next_bp_edges) == 1:  # No branching on the path
                eulerian_cycle.append(next_bp_edges[0])
                edges_cur[next_bp_edges[0]] = (
                    int(edges_cur[next_bp_edges[0]]) - 1
                )
                if edges_cur[next_bp_edges[0]] == 0:
                    del edges_cur[next_bp_edges[0]]
                if next_bp_edges[0][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[0][1]]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[0][1]]
                node_ = bp_edge.node1
                if node == node_:
                    node_ = bp_edge.node2
                eulerian_cycle.append(node_)
                last_seq_edge_idx = g.node_adjacencies[node_].sequence[0]
                eulerian_cycle.append(
                    EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)
                )
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_subcycle.append(
                        DirectedEdge(last_seq_edge_idx + 1, Strand.FORWARD)
                    )
                else:
                    last_edge_dir = "-"
                    eulerian_subcycle.append(
                        DirectedEdge(last_seq_edge_idx + 1, Strand.REVERSE)
                    )
                edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)] = (
                    int(edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)])
                    - 1
                )
                if edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)] == 0:
                    del edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)]
            else:
                r = random.randint(0, len(next_bp_edges) - 1)
                eulerian_cycle.append(next_bp_edges[r])
                edges_cur[next_bp_edges[r]] = (
                    int(edges_cur[next_bp_edges[r]]) - 1
                )
                if edges_cur[next_bp_edges[r]] == 0:
                    del edges_cur[next_bp_edges[r]]
                if next_bp_edges[r][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[r][1]]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[r][1]]
                node_ = bp_edge.node1
                if node == node_:
                    node_ = bp_edge.node2
                eulerian_cycle.append(node_)
                last_seq_edge_idx = g.node_adjacencies[node_].sequence[0]
                eulerian_cycle.append(
                    EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)
                )
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_subcycle.append(
                        DirectedEdge(last_seq_edge_idx + 1, Strand.FORWARD)
                    )
                else:
                    last_edge_dir = "-"
                    eulerian_subcycle.append(
                        DirectedEdge(last_seq_edge_idx + 1, Strand.REVERSE)
                    )
                edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)] = (
                    int(edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)])
                    - 1
                )
                if edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)] == 0:
                    del edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge_idx)]
        if valid == 1 and len(best_cycle) == 0:
            best_cycle = eulerian_subcycle
        path_metric = PathMetric(path_idxs=[], path_length=0, path_support=0)
        # check if the remaining path constraints are satisfied
        for pathi in range(len(path_constraints_next_cycle)):
            path_ = path_constraints_next_cycle[pathi]
            path0 = path_[0]
            s = 0
            for ei in range(len(eulerian_cycle) - 1):
                obj = eulerian_cycle[ei]
                if obj == path0:
                    s_ = 1
                    for i in range(len(path_)):
                        if (
                            eulerian_cycle[:-1][
                                (ei + i) % (len(eulerian_cycle) - 1)
                            ]
                            != path_[i]
                        ):
                            s_ = 0
                            break
                    if s_ == 1:
                        s = 1
                        break
                    s_ = 1
                    for i in range(len(path_)):
                        if eulerian_cycle[:-1][ei - i] != path_[i]:
                            s_ = 0
                            break
                    if s_ == 1:
                        s = 1
                        break
            if s == 0 and valid == 1:
                path_metric.path_idxs.append(pathi)
                path_metric.path_length += len(path_)
                path_metric.path_support += path_constraints_support[pathi]
        if valid == 1 and len(path_metric.path_idxs) > 0:
            valid = -1
        if (
            valid != 0
            and (
                len(path_metric.path_idxs)
                < len(unsatisfied_path_metric.path_idxs)
            )
            or (
                len(path_metric.path_idxs)
                == len(unsatisfied_path_metric.path_idxs)
                and path_metric.path_length
                < unsatisfied_path_metric.path_length
            )
            or (
                len(path_metric.path_idxs)
                == len(unsatisfied_path_metric.path_idxs)
                and path_metric.path_length
                == unsatisfied_path_metric.path_length
                and path_metric.path_support
                < unsatisfied_path_metric.path_support
            )
        ):
            unsatisfied_path_metric.path_idxs = path_metric.path_idxs
            unsatisfied_path_metric.path_length = path_metric.path_length
            unsatisfied_path_metric.path_support = path_metric.path_support
            best_cycle = eulerian_subcycle
    if len(unsatisfied_path_metric.path_idxs) == 0:
        logger.debug("Cycle satisfies all subpath constraints.")
    else:
        logger.debug("The following path constraints are not satisfied:")
        for pathi in unsatisfied_path_metric.path_idxs:
            logger.debug(f"{path_constraints_next_cycle[pathi]}")
    return best_cycle


def eulerian_path_t(
    g: BreakpointGraph,
    edges_next_path: OptimizationWalk,
    path_constraints_next_path: list[Walk],
    path_constraints_support: list[int],
) -> DirectedWalk:
    """Return an eulerian traversal of an s-t walk, represented by a list of
     directed edges.

    g: breakpoint graph (object)
    edges_next_path: subgraph induced by the s-t walk, as a dict that maps an edge to its multiplicity
            ***
            must include s and t in the dict
            ***
    path_constraints_next_path: list of subpath constraints to be satisfied,
            each as a list of alternating nodes and edges
            ***
            Note: the traversal may not satisfy all subpath constraints
            in case not all subpath constraints are satisfied, return the eulerian traversal satisfying the
            maximum number of subpath constraints
            ***
    path_constraints_support: num long reads supporting each subpath constraint
    """
    lseg = len(g.sequence_edges)
    endnode_list = list(g.endnode_adjacencies.keys())

    eulerian_path: Walk = []  # A path is edge - node list starting and ending with edges
    # Since Eulerian, there could be subcycles in the middle of a path
    eulerian_subpath: DirectedWalk = []  # Path in AA cycle format
    best_path: DirectedWalk = []  # Path in AA cycle format
    valid = 0
    num_trials = 0
    l = len(path_constraints_next_path)
    unsatisfied_path_metric = PathMetric(
        path_idxs=list(range(l)),
        path_length=100 * l,
        path_support=100 * max(path_constraints_support + [0]),
    )
    while valid <= 0 and num_trials < 1000:
        valid = 1
        num_trials += 1
        eulerian_path = []
        eulerian_subpath = []
        edges_cur = edges_next_path.copy()
        last_seq_edge = lseg
        last_edge_dir = "+"

        src_edge: EdgeId
        edge: EdgeId
        for edge in edges_cur:  # Start with the edge with smallest index
            if edge[0] == "s" or edge[0] == "t":
                src_edge = edge
                node = g.source_edges[edge[1]].node
                if len(eulerian_path) == 0:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path.append(EdgeId(EdgeType.TERMINAL, -1))
                    eulerian_path.append(node)
                    last_seq_edge = g.node_adjacencies[node].sequence[0]
                elif g.node_adjacencies[node].sequence[0] < last_seq_edge:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path[-1] = node
                    last_seq_edge = g.node_adjacencies[node].sequence[0]
            elif edge[0] == "ns" or edge[0] == "nt":
                src_edge = edge
                node = endnode_list[edge[1]]
                if len(eulerian_path) == 0:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path.append(EdgeId(EdgeType.TERMINAL, -1))
                    eulerian_path.append(node)
                    last_seq_edge = g.node_adjacencies[node].sequence[0]
                elif g.node_adjacencies[node].sequence[0] < last_seq_edge:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path[-1] = node
                    last_seq_edge = g.node_adjacencies[node].sequence[0]
        del edges_cur[src_edge]
        eulerian_path.append(EdgeId(EdgeType.SEQUENCE, last_seq_edge))
        if last_edge_dir == "+":
            eulerian_subpath.append(
                DirectedEdge(last_seq_edge + 1, Strand.FORWARD)
            )
        else:
            eulerian_subpath.append(
                DirectedEdge(last_seq_edge + 1, Strand.REVERSE)
            )
        edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)] = (
            int(edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)]) - 1
        )
        if edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)] == 0:
            del edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)]
        while len(edges_cur) > 0:
            seq_edge = g.sequence_edges[last_seq_edge]
            node = Node(seq_edge.chr, seq_edge.end, Strand.FORWARD)
            if last_edge_dir == "-":
                node = Node(seq_edge.chr, seq_edge.start, Strand.REVERSE)
            eulerian_path.append(node)
            if len(edges_cur) == 1 and (
                list(edges_cur.keys())[0][0] == "s"
                or list(edges_cur.keys())[0][0] == "ns"
                or list(edges_cur.keys())[0][0] == "t"
                or list(edges_cur.keys())[0][0] == "nt"
            ):
                eulerian_path.append(EdgeId(EdgeType.TERMINAL, -1))
                break
            # Since cycle, only consider discordant/concordant edges
            next_bp_edges: list[EdgeId] = []
            for ci in g.node_adjacencies[node].concordant:
                next_bp_edges.append(EdgeId(EdgeType.CONCORDANT, ci))
            for di in g.node_adjacencies[node].discordant:
                next_bp_edges.append(EdgeId(EdgeType.DISCORDANT, di))
            del_list = [
                i
                for i in range(len(next_bp_edges))
                if next_bp_edges[i] not in edges_cur
            ]
            for i in del_list[::-1]:
                del next_bp_edges[i]
            if len(next_bp_edges) == 0:
                valid = 0
                break

            bp_edge: ConcordantEdge | DiscordantEdge
            if len(next_bp_edges) == 1:  # No branching on the path
                eulerian_path.append(next_bp_edges[0])
                edges_cur[next_bp_edges[0]] = (
                    int(edges_cur[next_bp_edges[0]]) - 1
                )
                if edges_cur[next_bp_edges[0]] == 0:
                    del edges_cur[next_bp_edges[0]]
                if next_bp_edges[0][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[0][1]]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[0][1]]
                node_ = bp_edge.node1
                if node == node_:
                    node_ = bp_edge.node2
                eulerian_path.append(node_)
                last_seq_edge = g.node_adjacencies[node_].sequence[0]
                eulerian_path.append(EdgeId(EdgeType.SEQUENCE, last_seq_edge))
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_subpath.append(
                        DirectedEdge(last_seq_edge + 1, Strand.FORWARD)
                    )
                else:
                    last_edge_dir = "-"
                    eulerian_subpath.append(
                        DirectedEdge(last_seq_edge + 1, Strand.REVERSE)
                    )
                edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)] = (
                    int(edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)]) - 1
                )
                if edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)] == 0:
                    del edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)]
            else:
                r = random.randint(0, len(next_bp_edges) - 1)
                eulerian_path.append(next_bp_edges[r])
                edges_cur[next_bp_edges[r]] = (
                    int(edges_cur[next_bp_edges[r]]) - 1
                )
                if edges_cur[next_bp_edges[r]] == 0:
                    del edges_cur[next_bp_edges[r]]
                if next_bp_edges[r][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[r][1]]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[r][1]]
                node_ = bp_edge.node1
                if node == node_:
                    node_ = bp_edge.node2
                eulerian_path.append(node_)
                last_seq_edge = g.node_adjacencies[node_].sequence[0]
                eulerian_path.append(EdgeId(EdgeType.SEQUENCE, last_seq_edge))
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_subpath.append(
                        DirectedEdge(last_seq_edge + 1, Strand.FORWARD)
                    )
                else:
                    last_edge_dir = "-"
                    eulerian_subpath.append(
                        DirectedEdge(last_seq_edge + 1, Strand.REVERSE)
                    )
                edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)] = (
                    int(edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)]) - 1
                )
                if edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)] == 0:
                    del edges_cur[EdgeId(EdgeType.SEQUENCE, last_seq_edge)]
        if valid == 1 and len(best_path) == 0:
            best_path = eulerian_subpath
        path_metric = PathMetric(
            path_idxs=[],
            path_length=0,
            path_support=0,
        )
        # check if the remaining path constraints are satisfied
        for pathi in range(len(path_constraints_next_path)):
            path_ = path_constraints_next_path[pathi]
            s = 0
            for ei in range(2, len(eulerian_path) - 1 - len(path_)):
                if (
                    eulerian_path[ei : ei + len(path_)] == path_[:]
                    or eulerian_path[ei : ei + len(path_)] == path_[::-1]
                ):
                    s = 1
                    break
            if s == 0 and valid == 1:
                path_metric.path_idxs.append(pathi)
                path_metric.path_length += len(path_)
                path_metric.path_support += path_constraints_support[pathi]
        if valid == 1 and len(path_metric.path_idxs) > 0:
            valid = -1
        if (
            valid != 0
            and (
                len(path_metric.path_idxs)
                < len(unsatisfied_path_metric.path_idxs)
            )
            or (
                len(path_metric.path_idxs)
                == len(unsatisfied_path_metric.path_idxs)
                and path_metric.path_length
                < unsatisfied_path_metric.path_length
            )
            or (
                len(path_metric.path_idxs)
                == len(unsatisfied_path_metric.path_idxs)
                and path_metric.path_length
                == unsatisfied_path_metric.path_length
                and path_metric.path_support
                < unsatisfied_path_metric.path_support
            )
        ):
            unsatisfied_path_metric.path_idxs = path_metric.path_idxs
            unsatisfied_path_metric.path_length = path_metric.path_length
            unsatisfied_path_metric.path_support = path_metric.path_support
            best_path = eulerian_subpath
    if len(unsatisfied_path_metric.path_idxs) == 0:
        logger.debug("Path satisfies all subpath constraints.")
    else:
        logger.debug("The following path constraints are not satisfied:")
        for pathi in unsatisfied_path_metric.path_idxs:
            logger.debug(f"{path_constraints_next_path[pathi]}")
    return best_path
