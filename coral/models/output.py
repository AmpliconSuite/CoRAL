from __future__ import annotations

import logging
import random
from typing import Any, Dict

from coral import constants
from coral.breakpoint import infer_breakpoint_graph
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.constants import CHR_TAG_TO_IDX

logger = logging.getLogger(__name__)


def eulerian_cycle_t(
    g: BreakpointGraph,
    edges_next_cycle,
    path_constraints_next_cycle,
    path_constraints_support,
):
    """Return an eulerian traversal of a cycle, represented by a dict of edges

    g: breakpoint graph (object)
    edges_next_cycle: subgraph induced by the cycle, as a dict that maps an edge to its multiplicity
    path_constraints_next_cycle: list of subpath constraints to be satisfied,
            each as a list of alternating nodes and edges
            ***
            Note: the traversal may not satisfy all subpath constraints
            in case not all subpath constraints are satisfied, return the eulerian traversal satisfying the
            maximum number of subpath constraints
            ***
    path_constraints_support: num long reads supporting each subpath constraint
    """
    lseg = len(g.sequence_edges)

    eulerian_cycle: list[
        Any
    ] = []  # A cycle is edge - node list starting and ending with the same edge
    # Since Eulerian, there could be subcycles in the middle of a cycle
    eulerian_cycle_: list[Any] = []  # Cycle in AA cycle format
    best_cycle: list[Any] = []  # Cycle in AA cycle format
    valid = 0
    num_trials = 0
    l = len(path_constraints_next_cycle)
    unsatisfied_path_metric = [
        range(l),
        100 * l,
        100 * max(path_constraints_support + [0]),
    ]
    while valid <= 0 and num_trials < 1000:
        valid = 1
        num_trials += 1
        eulerian_cycle = []
        eulerian_cycle_ = []
        edges_cur = edges_next_cycle.copy()
        last_seq_edge = lseg  # Start with the edge with smallest index and on the positive strand
        for edge in edges_cur.keys():
            if edge[0] == "e":
                last_seq_edge = min(last_seq_edge, edge[1])
        last_edge_dir = "+"
        eulerian_cycle.append(("s", last_seq_edge))
        eulerian_cycle_.append(str(last_seq_edge + 1) + "+")
        while len(edges_cur) > 0:
            seq_edge = g.sequence_edges[last_seq_edge]
            node = (seq_edge[0], seq_edge[2], "+")
            if last_edge_dir == "-":
                node = (seq_edge[0], seq_edge[1], "-")
            eulerian_cycle.append(node)
            next_bp_edges = []  # Since cycle, only consider discordant edges and concordant edges
            for ci in g.nodes[node][1]:
                next_bp_edges.append(("c", ci))
            for di in g.nodes[node][2]:
                next_bp_edges.append(("d", di))
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
            if len(next_bp_edges) == 1:  # No branching on the path
                eulerian_cycle.append(next_bp_edges[0])
                edges_cur[next_bp_edges[0]] = (
                    int(edges_cur[next_bp_edges[0]]) - 1
                )
                if edges_cur[next_bp_edges[0]] == 0:
                    del edges_cur[next_bp_edges[0]]
                bp_edge = []
                if next_bp_edges[0][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[0][1]][:6]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[0][1]][:6]
                node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
                if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
                    node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
                eulerian_cycle.append(node_)
                last_seq_edge = g.nodes[node_][0][0]
                eulerian_cycle.append(("s", last_seq_edge))
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_cycle_.append(str(last_seq_edge + 1) + "+")
                else:
                    last_edge_dir = "-"
                    eulerian_cycle_.append(str(last_seq_edge + 1) + "-")
                edges_cur[("e", last_seq_edge)] = (
                    int(edges_cur[("e", last_seq_edge)]) - 1
                )
                if edges_cur[("e", last_seq_edge)] == 0:
                    del edges_cur[("e", last_seq_edge)]
            else:
                r = random.randint(0, len(next_bp_edges) - 1)
                eulerian_cycle.append(next_bp_edges[r])
                edges_cur[next_bp_edges[r]] = (
                    int(edges_cur[next_bp_edges[r]]) - 1
                )
                if edges_cur[next_bp_edges[r]] == 0:
                    del edges_cur[next_bp_edges[r]]
                bp_edge = []
                if next_bp_edges[r][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[r][1]][:6]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[r][1]][:6]
                node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
                if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
                    node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
                eulerian_cycle.append(node_)
                last_seq_edge = g.nodes[node_][0][0]
                eulerian_cycle.append(("s", last_seq_edge))
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_cycle_.append(str(last_seq_edge + 1) + "+")
                else:
                    last_edge_dir = "-"
                    eulerian_cycle_.append(str(last_seq_edge + 1) + "-")
                edges_cur[("e", last_seq_edge)] = (
                    int(edges_cur[("e", last_seq_edge)]) - 1
                )
                if edges_cur[("e", last_seq_edge)] == 0:
                    del edges_cur[("e", last_seq_edge)]
        if valid == 1 and len(best_cycle) == 0:
            best_cycle = eulerian_cycle_
        path_metric: list[Any] = [[], 0, 0]
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
                path_metric[0].append(pathi)
                path_metric[1] += len(path_)
                path_metric[2] += path_constraints_support[pathi]
        if valid == 1 and len(path_metric[0]) > 0:
            valid = -1
        if (
            valid != 0
            and (len(path_metric[0]) < len(unsatisfied_path_metric[0]))
            or (
                len(path_metric[0]) == len(unsatisfied_path_metric[0])
                and path_metric[1] < unsatisfied_path_metric[1]
            )
            or (
                len(path_metric[0]) == len(unsatisfied_path_metric[0])
                and path_metric[1] == unsatisfied_path_metric[1]
                and path_metric[2] < unsatisfied_path_metric[2]
            )
        ):
            unsatisfied_path_metric[0] = path_metric[0]
            unsatisfied_path_metric[1] = path_metric[1]
            unsatisfied_path_metric[2] = path_metric[2]
            best_cycle = eulerian_cycle_
    if len(unsatisfied_path_metric[0]) == 0:
        logger.debug(
            "Cycle satisfies all subpath constraints.",
        )
    else:
        logger.debug("The following path constraints are not satisfied:")
        for pathi in unsatisfied_path_metric[0]:
            logger.debug(f"{path_constraints_next_cycle[pathi]}")
    return best_cycle


def eulerian_path_t(
    g: BreakpointGraph,
    edges_next_path,
    path_constraints_next_path,
    path_constraints_support,
):
    """Return an eulerian traversal of an s-t walk, represented by a dict of edges

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
    endnode_list = [node for node in g.endnodes.keys()]

    eulerian_path: list[
        Any
    ] = []  # A path is edge - node list starting and ending with edges
    # Since Eulerian, there could be subcycles in the middle of a path
    eulerian_path_: list[Any] = []  # Path in AA cycle format
    best_path: list[Any] = []  # Path in AA cycle format
    valid = 0
    num_trials = 0
    l = len(path_constraints_next_path)
    unsatisfied_path_metric = [
        range(l),
        100 * l,
        100 * max(path_constraints_support + [0]),
    ]
    while valid <= 0 and num_trials < 1000:
        valid = 1
        num_trials += 1
        eulerian_path = []
        eulerian_path_ = []
        edges_cur = edges_next_path.copy()
        src_edge = ()
        last_seq_edge = lseg
        last_edge_dir = "+"
        for edge in edges_cur.keys():  # Start with the edge with smallest index
            if edge[0] == "s" or edge[0] == "t":
                src_edge = edge
                node = (
                    g.source_edges[edge[1]][3],
                    g.source_edges[edge[1]][4],
                    g.source_edges[edge[1]][5],
                )
                if len(eulerian_path) == 0:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path.append(("$", -1))
                    eulerian_path.append(node)
                    last_seq_edge = g.nodes[node][0][0]
                elif g.nodes[node][0][0] < last_seq_edge:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path[-1] = node
                    last_seq_edge = g.nodes[node][0][0]
            elif edge[0] == "ns" or edge[0] == "nt":
                src_edge = edge
                node = endnode_list[edge[1]]
                if len(eulerian_path) == 0:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path.append(("$", -1))
                    eulerian_path.append(node)
                    last_seq_edge = g.nodes[node][0][0]
                elif g.nodes[node][0][0] < last_seq_edge:
                    last_edge_dir = constants.INVERT_STRAND_DIRECTION[node[2]]
                    eulerian_path[-1] = node
                    last_seq_edge = g.nodes[node][0][0]
        del edges_cur[src_edge]
        eulerian_path.append(("s", last_seq_edge))
        if last_edge_dir == "+":
            eulerian_path_.append(str(last_seq_edge + 1) + "+")
        else:
            eulerian_path_.append(str(last_seq_edge + 1) + "-")
        edges_cur[("e", last_seq_edge)] = (
            int(edges_cur[("e", last_seq_edge)]) - 1
        )
        if edges_cur[("e", last_seq_edge)] == 0:
            del edges_cur[("e", last_seq_edge)]
        while len(edges_cur) > 0:
            seq_edge = g.sequence_edges[last_seq_edge]
            node = (seq_edge[0], seq_edge[2], "+")
            if last_edge_dir == "-":
                node = (seq_edge[0], seq_edge[1], "-")
            eulerian_path.append(node)
            if len(edges_cur) == 1 and (
                list(edges_cur.keys())[0][0] == "s"
                or list(edges_cur.keys())[0][0] == "ns"
                or list(edges_cur.keys())[0][0] == "t"
                or list(edges_cur.keys())[0][0] == "nt"
            ):
                eulerian_path.append(("$", -1))
                break
            next_bp_edges = []  # Since cycle, only consider discordant edges and concordant edges
            for ci in g.nodes[node][1]:
                next_bp_edges.append(("c", ci))
            for di in g.nodes[node][2]:
                next_bp_edges.append(("d", di))
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
            if len(next_bp_edges) == 1:  # No branching on the path
                eulerian_path.append(next_bp_edges[0])
                edges_cur[next_bp_edges[0]] = (
                    int(edges_cur[next_bp_edges[0]]) - 1
                )
                if edges_cur[next_bp_edges[0]] == 0:
                    del edges_cur[next_bp_edges[0]]
                bp_edge = []
                if next_bp_edges[0][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[0][1]][:6]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[0][1]][:6]
                node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
                if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
                    node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
                eulerian_path.append(node_)
                last_seq_edge = g.nodes[node_][0][0]
                eulerian_path.append(("s", last_seq_edge))
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_path_.append(str(last_seq_edge + 1) + "+")
                else:
                    last_edge_dir = "-"
                    eulerian_path_.append(str(last_seq_edge + 1) + "-")
                edges_cur[("e", last_seq_edge)] = (
                    int(edges_cur[("e", last_seq_edge)]) - 1
                )
                if edges_cur[("e", last_seq_edge)] == 0:
                    del edges_cur[("e", last_seq_edge)]
            else:
                r = random.randint(0, len(next_bp_edges) - 1)
                eulerian_path.append(next_bp_edges[r])
                edges_cur[next_bp_edges[r]] = (
                    int(edges_cur[next_bp_edges[r]]) - 1
                )
                if edges_cur[next_bp_edges[r]] == 0:
                    del edges_cur[next_bp_edges[r]]
                bp_edge = []
                if next_bp_edges[r][0] == "c":
                    bp_edge = g.concordant_edges[next_bp_edges[r][1]][:6]
                else:
                    bp_edge = g.discordant_edges[next_bp_edges[r][1]][:6]
                node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
                if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
                    node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
                eulerian_path.append(node_)
                last_seq_edge = g.nodes[node_][0][0]
                eulerian_path.append(("s", last_seq_edge))
                if node_[2] == "-":
                    last_edge_dir = "+"
                    eulerian_path_.append(str(last_seq_edge + 1) + "+")
                else:
                    last_edge_dir = "-"
                    eulerian_path_.append(str(last_seq_edge + 1) + "-")
                edges_cur[("e", last_seq_edge)] = (
                    int(edges_cur[("e", last_seq_edge)]) - 1
                )
                if edges_cur[("e", last_seq_edge)] == 0:
                    del edges_cur[("e", last_seq_edge)]
        if valid == 1 and len(best_path) == 0:
            best_path = eulerian_path_
        path_metric: list[Any] = [[], 0, 0]
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
                path_metric[0].append(pathi)
                path_metric[1] += len(path_)
                path_metric[2] += path_constraints_support[pathi]
        if valid == 1 and len(path_metric[0]) > 0:
            valid = -1
        if (
            valid != 0
            and (len(path_metric[0]) < len(unsatisfied_path_metric[0]))
            or (
                len(path_metric[0]) == len(unsatisfied_path_metric[0])
                and path_metric[1] < unsatisfied_path_metric[1]
            )
            or (
                len(path_metric[0]) == len(unsatisfied_path_metric[0])
                and path_metric[1] == unsatisfied_path_metric[1]
                and path_metric[2] < unsatisfied_path_metric[2]
            )
        ):
            unsatisfied_path_metric[0] = path_metric[0]
            unsatisfied_path_metric[1] = path_metric[1]
            unsatisfied_path_metric[2] = path_metric[2]
            best_path = eulerian_path_
    if len(unsatisfied_path_metric[0]) == 0:
        logger.debug("Path satisfies all subpath constraints.")
    else:
        logger.debug("The following path constraints are not satisfied:")
        for pathi in unsatisfied_path_metric[0]:
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
        ai_amplicon, key=lambda ai: (CHR_TAG_TO_IDX[ai.chr_tag], ai.start)
    )
    for ai in ai_amplicon:
        fp.write(
            f"Interval\t{interval_num}\t{ai.chr_tag}\t{ai.start}\t{ai.end}\n"
        )
        interval_num += 1

    fp.write("List of cycle segments\n")
    for seqi in range(len(bb.lr_graph[amplicon_idx].sequence_edges)):
        sseg = bb.lr_graph[amplicon_idx].sequence_edges[seqi]
        fp.write(f"Segment\t{seqi + 1}\t{sseg[0]}\t{sseg[1]}\t{sseg[2]}\n")
    if output_all_paths:
        fp.write("List of all subpath constraints\n")
        for pathi in range(len(bb.path_constraints[amplicon_idx][0])):
            fp.write("Path constraint\t%d\t" % (pathi + 1))
            path_ = bb.path_constraints[amplicon_idx][0][pathi]
            if path_[0][1] > path_[-1][1]:
                path_ = path_[::-1]
            for i in range(len(path_)):
                if i % 4 == 0:
                    if i < len(path_) - 1:
                        if path_[i + 1][2] == "+":
                            fp.write("%d+," % (path_[i][1] + 1))
                        else:
                            fp.write("%d-," % (path_[i][1] + 1))
                    elif path_[i - 1][2] == "+":
                        fp.write("%d-\t" % (path_[i][1] + 1))
                    else:
                        fp.write("%d+\t" % (path_[i][1] + 1))
            fp.write(
                "Support=%d\n" % (bb.path_constraints[amplicon_idx][1][pathi])
            )
    else:
        fp.write("List of longest subpath constraints\n")
        path_constraint_indices_ = []
        for paths in (
            bb.path_constraints_satisfied[amplicon_idx][0]
            + bb.path_constraints_satisfied[amplicon_idx][1]
        ):
            for pathi in paths:
                if pathi not in path_constraint_indices_:
                    path_constraint_indices_.append(pathi)
        for constraint_i in range(
            len(bb.longest_path_constraints[amplicon_idx][1])
        ):
            fp.write("Path constraint\t%d\t" % (constraint_i + 1))
            pathi = bb.longest_path_constraints[amplicon_idx][1][constraint_i]
            path_ = bb.path_constraints[amplicon_idx][0][pathi]
            if path_[0][1] > path_[-1][1]:
                path_ = path_[::-1]
            for i in range(len(path_)):
                if i % 4 == 0:
                    if i < len(path_) - 1:
                        if path_[i + 1][2] == "+":
                            fp.write("%d+," % (path_[i][1] + 1))
                        else:
                            fp.write("%d-," % (path_[i][1] + 1))
                    elif path_[i - 1][2] == "+":
                        fp.write("%d-\t" % (path_[i][1] + 1))
                    else:
                        fp.write("%d+\t" % (path_[i][1] + 1))
            fp.write(
                "Support=%d\t"
                % (bb.longest_path_constraints[amplicon_idx][2][constraint_i]),
            )
            if constraint_i in path_constraint_indices_:
                fp.write("Satisfied\n")
            else:
                fp.write("Unsatisfied\n")

    # sort cycles according to weights
    cycle_indices = sorted(
        [
            (0, i)
            for i in range(len(bb.walk_weights_by_amplicon[amplicon_idx][0]))
        ]
        + [
            (1, i)
            for i in range(len(bb.walk_weights_by_amplicon[amplicon_idx][1]))
        ],
        key=lambda item: bb.walk_weights_by_amplicon[amplicon_idx][item[0]][
            item[1]
        ],
        reverse=True,
    )

    print(cycle_indices)
    for cycle_i in cycle_indices:
        if cycle_i[0] == 0:  # cycles
            logger.debug(
                f"Traversing next cycle, CN = {bb.walk_weights_by_amplicon[amplicon_idx][cycle_i[0]][cycle_i[1]]}"
            )
            path_constraints_satisfied_cycle = []
            path_constraints_support_cycle = []
            for pathi in bb.path_constraints_satisfied[amplicon_idx][
                cycle_i[0]
            ][cycle_i[1]]:
                pathi_ = bb.longest_path_constraints[amplicon_idx][1][pathi]
                path_constraints_satisfied_cycle.append(
                    bb.path_constraints[amplicon_idx][0][pathi_],
                )
                path_constraints_support_cycle.append(
                    bb.longest_path_constraints[amplicon_idx][2][pathi],
                )
            cycle_seg_list = eulerian_cycle_t(
                bb.lr_graph[amplicon_idx],
                bb.walks_by_amplicon[amplicon_idx][cycle_i[0]][cycle_i[1]],
                path_constraints_satisfied_cycle,
                path_constraints_support_cycle,
            )
            assert cycle_seg_list[0] == cycle_seg_list[-1]
            fp.write("Cycle=%d;" % (cycle_indices.index(cycle_i) + 1))
            fp.write(
                "Copy_count=%s;"
                % str(
                    bb.walk_weights_by_amplicon[amplicon_idx][cycle_i[0]][
                        cycle_i[1]
                    ]
                ),
            )
            fp.write("Segments=")
            for segi in range(len(cycle_seg_list) - 2):
                fp.write(
                    f"{int(cycle_seg_list[segi][:-1])}{cycle_seg_list[segi][-1]},"
                )
            fp.write(
                "%d%s" % (int(cycle_seg_list[-2][:-1]), cycle_seg_list[-2][-1])
            )
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
                pathi_ = bb.longest_path_constraints[amplicon_idx][1][pathi]
                path_constraints_satisfied_path.append(
                    bb.path_constraints[amplicon_idx][0][pathi_],
                )
                path_constraints_support_path.append(
                    bb.longest_path_constraints[amplicon_idx][2][pathi],
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
            fp.write("Cycle=%d;" % (cycle_indices.index(cycle_i) + 1))
            fp.write(
                "Copy_count=%s;"
                % str(
                    bb.walk_weights_by_amplicon[amplicon_idx][cycle_i[0]][
                        cycle_i[1]
                    ]
                ),
            )
            fp.write("Segments=0+,")
            for segi in range(len(cycle_seg_list) - 1):
                fp.write(
                    "%d%s,"
                    % (int(cycle_seg_list[segi][:-1]), cycle_seg_list[segi][-1])
                )

            fp.write(
                "%d%s,0-"
                % (int(cycle_seg_list[-1][:-1]), cycle_seg_list[-1][-1])
            )
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
):
    logger.info(f"Outputting solution info for all amplicons.")

    fp = open(f"{output_dir}/amplicon_summary.txt", "w")
    fp.write(
        f"{sum(was_amplicon_solved.values())}/{len(bb.lr_graph)} amplicons solved.\n"
    )
    fp.write("--------------------------------------------------------------\n")
    for amplicon_idx in range(len(bb.lr_graph)):
        fp.write(f"Amplicon {amplicon_idx + 1}: ")
        if not was_amplicon_solved[amplicon_idx]:
            fp.write("UNSOLVED.\n")
            continue
        else:
            fp.write("Solved.\n")
    # TODO: mirror AA summary output
    fp.close()
