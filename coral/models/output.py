from __future__ import annotations

import io
import logging
import random
from typing import Dict

from coral import constants
from coral.breakpoint import infer_breakpoint_graph
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.constants import CHR_TAG_TO_IDX
from coral.datatypes import (
    AmpliconWalk,
    ConcordantEdge,
    DirectedEdge,
    DirectedWalk,
    DiscordantEdge,
    EdgeId,
    EdgeType,
    Node,
    PathConstraint,
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
    l = len(path_constraints_next_cycle)
    unsatisfied_path_metric = PathMetric(
        path_idxs=list(range(l)),
        path_length=100 * l,
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
    edges_next_path: AmpliconWalk,
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


def output_all_cycles(
    bb: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
    output_dir: str,
    output_all_paths: bool = False,
) -> None:
    """Write the result from cycle decomposition into *.cycles files"""
    for amplicon_idx in range(len(bb.lr_graph)):
        output_amplicon_cycles(amplicon_idx, bb, output_dir, output_all_paths)


def output_amplicon_cycles(
    amplicon_idx: int,
    bb: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
    output_dir: str,
    output_all_paths: bool = False,
) -> None:
    """Write the result from cycle decomposition into *.cycles files"""
    logger.info(f"Output cycles for amplicon {amplicon_idx+1}.")
    fp = open(f"{output_dir}/amplicon{amplicon_idx + 1}_cycles.txt", "w")
    interval_num = 1
    ai_amplicon = [
        ai
        for ai in bb.amplicon_intervals
        if bb.ccid2id[ai.amplicon_id] == amplicon_idx + 1
    ]
    ai_amplicon = sorted(
        ai_amplicon, key=lambda ai: (CHR_TAG_TO_IDX[ai.chr], ai.start)
    )
    for ai in ai_amplicon:
        fp.write(f"Interval\t{interval_num}\t{ai.chr}\t{ai.start}\t{ai.end}\n")
        interval_num += 1

    fp.write("List of cycle segments\n")
    for seqi in range(len(bb.lr_graph[amplicon_idx].sequence_edges)):
        sseg = bb.lr_graph[amplicon_idx].sequence_edges[seqi]
        fp.write(f"Segment\t{seqi + 1}\t{sseg.chr}\t{sseg.start}\t{sseg.end}\n")

    satisfied_path_constraints = bb.path_constraints_satisfied[amplicon_idx]
    longest_path_constraints = bb.longest_path_constraints[amplicon_idx]

    if output_all_paths:
        fp.write("List of all subpath constraints\n")
        for pathi in range(len(bb.path_constraints[amplicon_idx])):
            write_path_constraint_to_file(
                pathi, bb.path_constraints[amplicon_idx][pathi], fp
            )
    else:
        fp.write("List of longest subpath constraints\n")
        path_constraint_indices_ = []
        for paths in (
            satisfied_path_constraints.cycles + satisfied_path_constraints.paths
        ):
            for pathi in paths:
                if pathi not in path_constraint_indices_:
                    path_constraint_indices_.append(pathi)
        for constraint_i, finalized_pc in enumerate(longest_path_constraints):
            write_path_constraint_to_file(
                constraint_i,
                bb.path_constraints[amplicon_idx][finalized_pc.pc_idx],
                fp,
            )
            if constraint_i in path_constraint_indices_:
                fp.write("Satisfied\n")
            else:
                fp.write("Unsatisfied\n")

    walk_weights = bb.walk_weights_by_amplicon[amplicon_idx]

    # sort walks according to weights
    walk_indices = sorted(
        [(0, i) for i in range(len(walk_weights.cycles))]
        + [(1, i) for i in range(len(walk_weights.paths))],
        key=lambda item: walk_weights[item[0]][item[1]],
        reverse=True,
    )

    for cycle_i in walk_indices:
        if cycle_i[0] == 0:  # cycles
            logger.debug(
                f"Traversing next cycle, CN = {walk_weights.cycles[cycle_i[1]]}"
            )
            path_constraints_satisfied_cycle = []
            path_constraints_support_cycle = []
            for pathi in satisfied_path_constraints.cycles[cycle_i[1]]:
                pathi_ = bb.longest_path_constraints[amplicon_idx][pathi].pc_idx
                path_constraints_satisfied_cycle.append(
                    bb.path_constraints[amplicon_idx][pathi_].path
                )
                path_constraints_support_cycle.append(
                    bb.longest_path_constraints[amplicon_idx][pathi].support,
                )
            cycle_seg_list = eulerian_cycle_t(
                bb.lr_graph[amplicon_idx],
                bb.walks_by_amplicon[amplicon_idx][cycle_i[0]][cycle_i[1]],
                path_constraints_satisfied_cycle,
                path_constraints_support_cycle,
            )
            assert cycle_seg_list[0] == cycle_seg_list[-1]
            fp.write(f"Cycle={walk_indices.index(cycle_i) + 1};")
            fp.write(f"Copy_count={walk_weights.cycles[cycle_i[1]]};")
            fp.write("Segments=")
            for segi in range(len(cycle_seg_list) - 2):
                fp.write(f"{cycle_seg_list[segi][0]}{cycle_seg_list[segi][1]},")
            fp.write(f"{cycle_seg_list[-2][0]}{cycle_seg_list[-2][1]}")
            if not output_all_paths:
                fp.write(";Path_constraints_satisfied=")
                for pathi in range(
                    len(
                        bb.path_constraints_satisfied[amplicon_idx][cycle_i[0]][
                            cycle_i[1]
                        ]
                    )
                    - 1,
                ):
                    fp.write(
                        "%d,"
                        % (
                            bb.path_constraints_satisfied[amplicon_idx][
                                cycle_i[0]
                            ][cycle_i[1]][pathi]
                            + 1
                        ),
                    )
                if (
                    len(
                        bb.path_constraints_satisfied[amplicon_idx][cycle_i[0]][
                            cycle_i[1]
                        ]
                    )
                    > 0
                ):
                    fp.write(
                        "%d\n"
                        % (
                            bb.path_constraints_satisfied[amplicon_idx][
                                cycle_i[0]
                            ][cycle_i[1]][-1]
                            + 1
                        ),
                    )
                else:
                    fp.write("\n")
            else:
                fp.write("\n")
        else:  # paths
            logger.debug(
                f"Traversing next path, CN = {bb.walk_weights_by_amplicon[amplicon_idx][cycle_i[0]][cycle_i[1]]}"
            )
            path_constraints_satisfied_path = []
            path_constraints_support_path = []
            for pathi in bb.path_constraints_satisfied[amplicon_idx][
                cycle_i[0]
            ][cycle_i[1]]:
                pathi_ = bb.longest_path_constraints[amplicon_idx][pathi].pc_idx
                path_constraints_satisfied_path.append(
                    bb.path_constraints[amplicon_idx][pathi_].path,
                )
                path_constraints_support_path.append(
                    bb.longest_path_constraints[amplicon_idx][pathi].support
                )
            cycle_seg_list = eulerian_path_t(
                bb.lr_graph[amplicon_idx],
                bb.walks_by_amplicon[amplicon_idx][cycle_i[0]][cycle_i[1]],
                path_constraints_satisfied_path,
                path_constraints_support_path,
            )
            print(
                cycle_seg_list,
                bb.walks_by_amplicon[amplicon_idx][cycle_i[0]][cycle_i[1]],
            )
            fp.write("Cycle=%d;" % (walk_indices.index(cycle_i) + 1))
            fp.write(
                f"Copy_count={bb.walk_weights_by_amplicon[amplicon_idx][cycle_i[0]][cycle_i[1]]};"
            )
            fp.write("Segments=0+,")
            for segi in range(len(cycle_seg_list) - 1):
                fp.write(f"{cycle_seg_list[segi][0]}{cycle_seg_list[segi][1]},")

            fp.write(f"{cycle_seg_list[-1][0]}{cycle_seg_list[-1][1]},0-")
            if not output_all_paths:
                fp.write(";Path_constraints_satisfied=")
                for pathi in range(
                    len(
                        bb.path_constraints_satisfied[amplicon_idx][cycle_i[0]][
                            cycle_i[1]
                        ]
                    )
                    - 1,
                ):
                    fp.write(
                        "%d,"
                        % (
                            bb.path_constraints_satisfied[amplicon_idx][
                                cycle_i[0]
                            ][cycle_i[1]][pathi]
                            + 1
                        ),
                    )
                if (
                    len(
                        bb.path_constraints_satisfied[amplicon_idx][cycle_i[0]][
                            cycle_i[1]
                        ]
                    )
                    > 0
                ):
                    fp.write(
                        "%d\n"
                        % (
                            bb.path_constraints_satisfied[amplicon_idx][
                                cycle_i[0]
                            ][cycle_i[1]][-1]
                            + 1
                        ),
                    )
                else:
                    fp.write("\n")
            else:
                fp.write("\n")
    fp.close()


def output_summary_amplicon_stats(
    was_amplicon_solved: Dict[int, bool],
    bb: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
    output_dir: str,
) -> None:
    logger.info("Outputting solution info for all amplicons.")

    fp = open(f"{output_dir}/amplicon_summary.txt", "w")
    fp.write(
        f"{sum(was_amplicon_solved.values())}/{len(bb.lr_graph)} amplicons "
        "solved.\n"
    )
    fp.write("--------------------------------------------------------------\n")
    for amplicon_idx in range(len(bb.lr_graph)):
        fp.write(f"Amplicon {amplicon_idx + 1}: ")
        if not was_amplicon_solved[amplicon_idx]:
            fp.write("UNSOLVED.\n")
            continue
        fp.write("Solved.\n")
    # TODO: mirror AA summary output
    fp.close()


def write_path_constraint_to_file(
    path_idx: int, path_constraint: PathConstraint, fp: io.TextIOWrapper
) -> None:
    fp.write(f"Path constraint\t{path_idx+1}\t")
    fp.write(path_constraint.to_file_str())
