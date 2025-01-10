"""Functions used for constructing subpath constraints"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from coral.breakpoint import breakpoint_utilities

edge_type_to_index = {"s": 0, "c": 1, "d": 2}


def valid_path(g, path):
    """Check if the input subpath constraint is valid in graph g
    A valid path must be an alternating sequence of nodes and edges;
            consist of >= 3 sequence edges;
            start and end with sequence edges.
            edges must alternate between sequence edges and breakpoint (concordant or discordant) edges

    g: breakpoint graph (object)
    path: a list of alternating nodes and edges

    Returns: True/False
    """
    if len(path) <= 3 or len(path) % 2 == 0:
        return False
    if path[0][0] != "s" or path[-1][0] != "s":
        return False
    for i in range(len(path)):
        if i % 2 == 0:
            if len(path[i]) != 2:
                return False
        else:
            if len(path[i]) != 3:
                return False
            e1 = path[i - 1]
            e2 = path[i + 1]
            try:
                if (e1[0] == "s" and e2[0] == "s") or (
                    e1[0] != "s" and e2[0] != "s"
                ):
                    return False
                if e1[1] not in g.nodes[path[i]][edge_type_to_index[e1[0]]]:
                    return False
                if e2[1] not in g.nodes[path[i]][edge_type_to_index[e2[0]]]:
                    return False
            except:
                return False
    return True


def alignment_to_path(g, rint, min_overlap=500):
    """Traverse through the input breakpoint graph to convert a single alignment to a path

    g: breakpoint graph (object)
    rint: alignment interval on the reference genome
    min_overlap: required overlap (in bp) between the alignment and the first/last sequence edge in the resulting path,
            default value is 500bp

    Returns: the resulting path as a list of alternating nodes and edges
    """
    seq_edge_list = []
    for segi in range(len(g.sequence_edges)):
        sseg = g.sequence_edges[segi]
        if breakpoint_utilities.interval_overlap(rint, sseg):
            seq_edge_list.append(segi)
    if len(seq_edge_list) == 0:
        return []
    seq_edge_list = sorted(
        seq_edge_list, key=lambda item: g.sequence_edges[item][1]
    )
    segi0 = seq_edge_list[0]
    if (
        len(seq_edge_list) > 1
        and min(g.sequence_edges[segi0][2], rint[2])
        - max(g.sequence_edges[segi0][1], rint[1])
        < min_overlap
    ):
        del seq_edge_list[0]
    segi0 = seq_edge_list[0]
    while len(seq_edge_list) > 1 and g.sequence_edges[segi0][7] < min_overlap:
        del seq_edge_list[0]
        segi0 = seq_edge_list[0]
    segi0 = seq_edge_list[-1]
    if (
        len(seq_edge_list) > 1
        and min(g.sequence_edges[segi0][2], rint[2])
        - max(g.sequence_edges[segi0][1], rint[1])
        < min_overlap
    ):  # need to parameterize this
        del seq_edge_list[-1]
    segi0 = seq_edge_list[-1]
    while len(seq_edge_list) > 1 and g.sequence_edges[segi0][7] < min_overlap:
        del seq_edge_list[-1]
        segi0 = seq_edge_list[-1]
    if len(seq_edge_list) <= 2:
        return []
    segi0 = seq_edge_list[0]
    node1 = (g.sequence_edges[segi0][0], g.sequence_edges[segi0][1], "-")
    segi0 = seq_edge_list[-1]
    node2 = (g.sequence_edges[segi0][0], g.sequence_edges[segi0][2], "+")
    path_ = traverse_through_sequence_edge(g, node1, node2)[1:-1]
    return path_


def chimeric_alignment_to_path_l(g, rints, ai, bp_node, min_overlap=500):
    """Given a breakpoint graph, a list of consecutive alignments, and an end node,
            return a traversal from the alignment indexed at ai to the end node

    g: breakpoint graph (object)
    rints: alignment intervals on the reference genome
    ai: index of the starting alignment
    bp_node: end node
    min_overlap: required overlap (in bp) between the alignment and the first/last sequence edge in the resulting path,
            default value is 500bp

    Returns: the resulting path as a list of alternating nodes and edges
            note that the resulting path (additionally) starts with a node
    """
    al = rints[ai]
    seq_edge_list = []
    for segi in range(len(g.sequence_edges)):
        sseg = g.sequence_edges[segi]
        if al[-1] == "+":
            if breakpoint_utilities.interval_overlap(al, sseg):
                seq_edge_list.append([segi, "+"])
        elif breakpoint_utilities.interval_overlap([al[0], al[2], al[1]], sseg):
            seq_edge_list.append([segi, "-"])
    if len(seq_edge_list) == 0:
        return []
    if seq_edge_list[0][1] == "+":
        seq_edge_list = sorted(
            seq_edge_list, key=lambda item: g.sequence_edges[item[0]][1]
        )
        segi0 = seq_edge_list[0][0]
        if (
            len(seq_edge_list) > 1
            and min(g.sequence_edges[segi0][2], al[2])
            - max(g.sequence_edges[segi0][1], al[1])
            < min_overlap
        ):
            del seq_edge_list[0]
        segi0 = seq_edge_list[0][0]
        while (
            len(seq_edge_list) > 0 and g.sequence_edges[segi0][7] < min_overlap
        ):
            del seq_edge_list[0]
            if len(seq_edge_list) > 0:
                segi0 = seq_edge_list[0][0]
        # check if the rightmost node connects to the breakpoint edge at index edi
        while len(seq_edge_list) > 0:
            segi_last = seq_edge_list[-1][0]
            rnode = (
                g.sequence_edges[segi_last][0],
                g.sequence_edges[segi_last][2],
                "+",
            )
            if rnode != bp_node:
                del seq_edge_list[-1]
            else:
                break
    else:
        seq_edge_list = sorted(
            seq_edge_list,
            key=lambda item: g.sequence_edges[item[0]][1],
            reverse=True,
        )
        segi0 = seq_edge_list[0][0]
        if (
            len(seq_edge_list) > 1
            and min(g.sequence_edges[segi0][2], al[1])
            - max(g.sequence_edges[segi0][1], al[2])
            < min_overlap
        ):
            del seq_edge_list[0]
        segi0 = seq_edge_list[0][0]
        while (
            len(seq_edge_list) > 0 and g.sequence_edges[segi0][7] < min_overlap
        ):
            del seq_edge_list[0]
            if len(seq_edge_list) > 0:
                segi0 = seq_edge_list[0][0]
        # check if the rightmost node connects to the breakpoint edge at index edi
        while len(seq_edge_list) > 0:
            segi_last = seq_edge_list[-1][0]
            rnode = (
                g.sequence_edges[segi_last][0],
                g.sequence_edges[segi_last][1],
                "-",
            )
            if rnode != bp_node:
                del seq_edge_list[-1]
            else:
                break
    if len(seq_edge_list) == 0:
        return []
    path_l = []
    for si in range(len(seq_edge_list)):
        path_l.append(("s", seq_edge_list[si][0]))
        if seq_edge_list[si][1] == "+":
            path_l.append(
                (
                    g.sequence_edges[seq_edge_list[si][0]][0],
                    g.sequence_edges[seq_edge_list[si][0]][2],
                    "+",
                ),
            )
        else:
            path_l.append(
                (
                    g.sequence_edges[seq_edge_list[si][0]][0],
                    g.sequence_edges[seq_edge_list[si][0]][1],
                    "-",
                ),
            )
        if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == "+":
            if (
                g.sequence_edges[seq_edge_list[si][0]][2] + 1
                == g.sequence_edges[seq_edge_list[si + 1][0]][1]
            ):
                for ci in range(len(g.concordant_edges)):
                    if (
                        g.concordant_edges[ci][0]
                        == g.sequence_edges[seq_edge_list[si][0]][0]
                        and g.sequence_edges[seq_edge_list[si][0]][2]
                        == g.concordant_edges[ci][1]
                        and g.sequence_edges[seq_edge_list[si + 1][0]][1]
                        == g.concordant_edges[ci][4]
                    ):
                        path_l.append(("c", ci))
                        path_l.append(
                            (
                                g.sequence_edges[seq_edge_list[si][0]][0],
                                g.sequence_edges[seq_edge_list[si + 1][0]][1],
                                "-",
                            ),
                        )
                        break
        if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == "-":
            if (
                g.sequence_edges[seq_edge_list[si][0]][1] - 1
                == g.sequence_edges[seq_edge_list[si + 1][0]][2]
            ):
                for ci in range(len(g.concordant_edges)):
                    if (
                        g.concordant_edges[ci][0]
                        == g.sequence_edges[seq_edge_list[si][0]][0]
                        and g.sequence_edges[seq_edge_list[si + 1][0]][2]
                        == g.concordant_edges[ci][1]
                        and g.sequence_edges[seq_edge_list[si][0]][1]
                        == g.concordant_edges[ci][4]
                    ):
                        path_l.append(("c", ci))
                        path_l.append(
                            (
                                g.sequence_edges[seq_edge_list[si][0]][0],
                                g.sequence_edges[seq_edge_list[si + 1][0]][2],
                                "+",
                            ),
                        )
                        break
    return path_l


def chimeric_alignment_to_path_r(g, rints, ai, bp_node, min_overlap=500):
    """Given a breakpoint graph, a list of consecutive alignments, and a start node,
            return a traversal from the starting node to the alignment indexed at ai

    g: breakpoint graph (object)
    rints: alignment intervals on the reference genome
    ai: index of the end alignment
    bp_node: start node
    min_overlap: required overlap (in bp) between the alignment and the first/last sequence edge in the resulting path,
            default value is 500bp

    Returns: the resulting path as a list of alternating nodes and edges
            note that the resulting path (additionally) ends with a node
    """
    ar = rints[ai]
    seq_edge_list = []
    for segi in range(len(g.sequence_edges)):
        sseg = g.sequence_edges[segi]
        if ar[-1] == "+":
            if breakpoint_utilities.interval_overlap(ar, sseg):
                seq_edge_list.append([segi, "+"])
        elif breakpoint_utilities.interval_overlap([ar[0], ar[2], ar[1]], sseg):
            seq_edge_list.append([segi, "-"])
    if len(seq_edge_list) == 0:
        return []
    if seq_edge_list[0][1] == "+":
        seq_edge_list = sorted(
            seq_edge_list, key=lambda item: g.sequence_edges[item[0]][1]
        )
        segi1 = seq_edge_list[-1][0]
        if (
            min(g.sequence_edges[segi1][2], ar[2])
            - max(g.sequence_edges[segi1][1], ar[1])
            < 500
        ):  # need to parameterize this
            del seq_edge_list[-1]
        if len(seq_edge_list) == 0:
            return []
        segi1 = seq_edge_list[-1][0]
        while len(seq_edge_list) > 0 and g.sequence_edges[segi1][7] < 500:
            del seq_edge_list[-1]
            if len(seq_edge_list) > 0:
                segi1 = seq_edge_list[-1][0]
        # check if the leftmost node connects to the breakpoint edge at index edi
        while len(seq_edge_list) > 0:
            segi_last = seq_edge_list[0][0]
            lnode = (
                g.sequence_edges[segi_last][0],
                g.sequence_edges[segi_last][1],
                "-",
            )
            if lnode != bp_node:
                del seq_edge_list[0]
            else:
                break
    else:
        seq_edge_list = sorted(
            seq_edge_list,
            key=lambda item: g.sequence_edges[item[0]][1],
            reverse=True,
        )
        segi1 = seq_edge_list[-1][0]
        if (
            min(g.sequence_edges[segi1][2], ar[1])
            - max(g.sequence_edges[segi1][1], ar[2])
            < 500
        ):  # need to parameterize this
            del seq_edge_list[-1]
        if len(seq_edge_list) == 0:
            return []
        segi1 = seq_edge_list[-1][0]
        while len(seq_edge_list) > 0 and g.sequence_edges[segi1][7] < 500:
            del seq_edge_list[-1]
            if len(seq_edge_list) > 0:
                segi1 = seq_edge_list[-1][0]
        while len(seq_edge_list) > 0:
            segi_last = seq_edge_list[0][0]
            lnode = (
                g.sequence_edges[segi_last][0],
                g.sequence_edges[segi_last][2],
                "+",
            )
            if lnode != bp_node:
                del seq_edge_list[0]
            else:
                break
    if len(seq_edge_list) == 0:
        return []
    path_r = []
    for si in range(len(seq_edge_list)):
        if seq_edge_list[si][1] == "+":
            path_r.append(
                (
                    g.sequence_edges[seq_edge_list[si][0]][0],
                    g.sequence_edges[seq_edge_list[si][0]][1],
                    "-",
                ),
            )
        else:
            path_r.append(
                (
                    g.sequence_edges[seq_edge_list[si][0]][0],
                    g.sequence_edges[seq_edge_list[si][0]][2],
                    "+",
                ),
            )
        path_r.append(("s", seq_edge_list[si][0]))
        if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == "+":
            if (
                g.sequence_edges[seq_edge_list[si][0]][2] + 1
                == g.sequence_edges[seq_edge_list[si + 1][0]][1]
            ):
                for ci in range(len(g.concordant_edges)):
                    if (
                        g.concordant_edges[ci][0]
                        == g.sequence_edges[seq_edge_list[si][0]][0]
                        and g.sequence_edges[seq_edge_list[si][0]][2]
                        == g.concordant_edges[ci][1]
                        and g.sequence_edges[seq_edge_list[si + 1][0]][1]
                        == g.concordant_edges[ci][4]
                    ):
                        path_r.append(
                            (
                                g.sequence_edges[seq_edge_list[si][0]][0],
                                g.sequence_edges[seq_edge_list[si][0]][2],
                                "+",
                            ),
                        )
                        path_r.append(("c", ci))
                        break
        if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == "-":
            if (
                g.sequence_edges[seq_edge_list[si][0]][1] - 1
                == g.sequence_edges[seq_edge_list[si + 1][0]][2]
            ):
                for ci in range(len(g.concordant_edges)):
                    if (
                        g.concordant_edges[ci][0]
                        == g.sequence_edges[seq_edge_list[si][0]][0]
                        and g.sequence_edges[seq_edge_list[si + 1][0]][2]
                        == g.concordant_edges[ci][1]
                        and g.sequence_edges[seq_edge_list[si][0]][1]
                        == g.concordant_edges[ci][4]
                    ):
                        path_r.append(
                            (
                                g.sequence_edges[seq_edge_list[si][0]][0],
                                g.sequence_edges[seq_edge_list[si][0]][1],
                                "-",
                            ),
                        )
                        path_r.append(("c", ci))
                        break
    return path_r


def chimeric_alignment_to_path_i(g, rints, ai1, ai2, di):
    """Given a breakpoint graph and a list of consecutive alignments,
            return a traversal from the alignment indexed at ai1 to the alignment indexed at ai2,
            through discordant edge indexed at di

    g: breakpoint graph (object)
    rints: alignment intervals on the reference genome
    ai1: index of the start alignment
    ai2: index of the end alignment
    di: index of discordant edge in g

    Returns: the resulting path as a list of alternating nodes and edges
    """
    path_ = [("d", di)]
    node1 = (
        g.discordant_edges[di][0],
        g.discordant_edges[di][1],
        g.discordant_edges[di][2],
    )
    node2 = (
        g.discordant_edges[di][3],
        g.discordant_edges[di][4],
        g.discordant_edges[di][5],
    )
    if ai1 > ai2:
        path_ = (
            chimeric_alignment_to_path_l(g, rints, ai2, node2)
            + path_
            + chimeric_alignment_to_path_r(g, rints, ai1, node1)
        )
    else:
        path_ = (
            chimeric_alignment_to_path_l(g, rints, ai1, node1)
            + path_
            + chimeric_alignment_to_path_r(g, rints, ai2, node2)
        )
    return path_


def traverse_through_sequence_edge(g, start_node, end_node):
    """Given a breakpoint graph and two nodes, return a traversal through sequence and concordant edges between the two nodes

    g: breakpoint graph (object)
    start_node: start node
    end_node: end node - start and end node must locate at different (left/right) ends on the corresponding sequence edges

    Returns: the resulting path as a list of alternating nodes and edges
            note that the resulting path (additionally) starts and ends with the given nodes
    """
    assert start_node[2] != end_node[2]
    path_ = [start_node]
    seqi = g.nodes[start_node][0][0]
    seq_edge = g.sequence_edges[seqi]
    next_end = (seq_edge[0], seq_edge[1], "-")
    if start_node[2] == "-":
        next_end = (seq_edge[0], seq_edge[2], "+")
    path_.append(("s", seqi))
    path_.append(next_end)
    while next_end != end_node:
        try:
            ci = g.nodes[next_end][1][0]
        except:
            return (
                path_  # ignore the alignments spanning two amplicon intervals
            )
        path_.append(("c", ci))
        cedge = g.concordant_edges[ci]
        next_start = (cedge[0], cedge[1], cedge[2])
        if next_start == next_end:
            next_start = (cedge[3], cedge[4], cedge[5])
        path_.append(next_start)
        seqi = g.nodes[next_start][0][0]
        seq_edge = g.sequence_edges[seqi]
        next_end = (seq_edge[0], seq_edge[1], "-")
        if next_start[2] == "-":
            next_end = (seq_edge[0], seq_edge[2], "+")
        path_.append(("s", seqi))
        path_.append(next_end)
    return path_


def chimeric_alignment_to_path(g, rints, ai_list, bp_list):
    """Convert chimeric alignments to path"""
    path_ = []
    lastnode = ()
    for i in range(len(bp_list)):
        di = bp_list[i]
        node1 = (
            g.discordant_edges[di][0],
            g.discordant_edges[di][1],
            g.discordant_edges[di][2],
        )
        node2 = (
            g.discordant_edges[di][3],
            g.discordant_edges[di][4],
            g.discordant_edges[di][5],
        )
        if ai_list[i][0] > ai_list[i][1]:
            if i == 0:
                path_ = chimeric_alignment_to_path_l(
                    g, rints, ai_list[i][1], node2
                ) + [
                    ("d", bp_list[i]),
                ]
                lastnode = node1
            else:
                path_ += traverse_through_sequence_edge(g, lastnode, node2)
                path_.append(("d", bp_list[i]))
                lastnode = node1
                if i == len(bp_list) - 1:
                    path_ += chimeric_alignment_to_path_r(
                        g, rints, ai_list[i][0], node1
                    )
        elif i == 0:
            path_ = chimeric_alignment_to_path_l(
                g, rints, ai_list[i][0], node1
            ) + [
                ("d", bp_list[i]),
            ]
            lastnode = node2
        else:
            path_ += traverse_through_sequence_edge(g, lastnode, node1)
            path_.append(("d", bp_list[i]))
            lastnode = node2
            if i == len(bp_list) - 1:
                path_ += chimeric_alignment_to_path_r(
                    g, rints, ai_list[i][1], node2
                )
    return path_


def longest_path_dict(
    path_constraints_: list[list[Any]],
) -> list[list[Any]]:
    """Convert paths from a list of alternating nodes and edges into a dict of edges.
    Only keep the longest paths, i.e., those which are not a subpath of any other path
    """
    res_paths: list[list[Any]] = [[], [], []]
    for pathi in range(len(path_constraints_[0])):
        path = path_constraints_[0][pathi]
        path_constraint: dict[int, int] = defaultdict(int)
        for ei in range(len(path)):
            if ei % 2 == 0:
                path_constraint[path[ei]] += 1
        res_paths[0].append(path_constraint)
        res_paths[1].append(pathi)
        res_paths[2].append(path_constraints_[1][pathi])
    for pathi in range(len(res_paths[0]))[::-1]:
        path_constraint = res_paths[0][pathi]
        subpath_flag = -1
        for pathi_ in range(len(res_paths[0])):
            path_constraint_ = res_paths[0][pathi_]
            s1 = 1
            for edge in path_constraint.keys():
                if (
                    edge not in path_constraint_
                    or path_constraint_[edge] < path_constraint[edge]
                ):
                    s1 = 0
                    break
            if s1 == 1 and pathi_ != pathi:
                subpath_flag = pathi_
                break
        if subpath_flag >= 0:
            del res_paths[0][pathi]
            del res_paths[1][pathi]
            res_paths[2][subpath_flag] = max(
                res_paths[2][subpath_flag], res_paths[2][pathi]
            )
            del res_paths[2][pathi]
    return res_paths
