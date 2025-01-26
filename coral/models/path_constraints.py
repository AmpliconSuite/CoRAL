"""Functions used for constructing subpath constraints"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from coral.breakpoint import breakpoint_utilities
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    BPIndexedAlignments,
    Edge,
    EdgeId,
    FinalizedPathConstraint,
    Node,
    PathConstraint,
    ReferenceInterval,
    SequenceEdge,
    Strand,
    Walk,
)

edge_type_to_index = {"s": 0, "c": 1, "d": 2}


def valid_path(g: BreakpointGraph, path: Walk) -> bool:
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
                if (
                    e1[1]
                    not in g.node_adjacencies[path[i]][
                        edge_type_to_index[e1[0]]
                    ]
                ):
                    return False
                if (
                    e2[1]
                    not in g.node_adjacencies[path[i]][
                        edge_type_to_index[e2[0]]
                    ]
                ):
                    return False
            except:
                return False
    return True


def alignment_to_path(
    g: BreakpointGraph, ref_intv: ReferenceInterval, min_overlap: int = 500
) -> Walk:
    """Traverse through the input breakpoint graph to convert a single alignment
      to a path.

    g: breakpoint graph (object)
    ref_intv: alignment interval on the reference genome
    min_overlap: required overlap (in bp) between the alignment and the
        first/last sequence edge in the resulting path, default value is 500bp

    Returns: the resulting path as a list of alternating nodes and edges
    """
    overlapping_seq_edges = sorted(
        [
            seq_edge
            for seq_edge in g.sequence_edges
            if ref_intv.does_overlap(seq_edge.interval)
            and (
                min(seq_edge.end, ref_intv.end)
                - max(seq_edge.start, ref_intv.start)
                >= min_overlap
            )
            and (seq_edge.gap >= min_overlap)
        ]
    )

    if len(overlapping_seq_edges) <= 2:
        return []

    node1 = overlapping_seq_edges[0].start_node
    node2 = overlapping_seq_edges[-1].end_node
    path_ = traverse_through_sequence_edge(g, node1, node2)[1:-1]
    return path_


def chimeric_alignment_to_path_l(
    g: BreakpointGraph,
    ref_intv: ReferenceInterval,
    bp_node: Node,
    min_overlap: int = 500,
) -> Walk:
    """Given a breakpoint graph, a reference alignment, and an end node,
            return a traversal from the alignment to the end node

    g: breakpoint graph (object)
    ref_intv: alignment interval on the reference genome
    bp_node: end node
    min_overlap: required overlap (in bp) between the alignment and the
        first/last sequence edge in the resulting path, default value is 500bp

    Returns: the resulting path as a list of alternating nodes and edges
            note that the resulting path (additionally) starts with a node
    """
    seq_edge_idxs: list[tuple[int, Strand]] = []
    for segi in range(len(g.sequence_edges)):
        sseg = g.sequence_edges[segi]
        if ref_intv.strand == Strand.FORWARD:
            if ref_intv.does_overlap(sseg.interval):
                seq_edge_idxs.append((segi, Strand.FORWARD))
        elif ref_intv.reverse.does_overlap(sseg.interval):
            seq_edge_idxs.append((segi, Strand.REVERSE))
    if len(seq_edge_idxs) == 0:
        return []
    # TODO: clean up redundancy below
    if seq_edge_idxs[0][1] == Strand.FORWARD:
        seq_edge_idxs = sorted(
            seq_edge_idxs, key=lambda item: g.sequence_edges[item[0]].start
        )
        segi0 = seq_edge_idxs[0][0]
        if (
            len(seq_edge_idxs) > 1
            and min(g.sequence_edges[segi0].end, ref_intv.end)
            - max(g.sequence_edges[segi0].start, ref_intv.start)
            < min_overlap
        ):
            del seq_edge_idxs[0]
        segi0 = seq_edge_idxs[0][0]
        while (
            len(seq_edge_idxs) > 0 and g.sequence_edges[segi0].gap < min_overlap
        ):
            del seq_edge_idxs[0]
            if len(seq_edge_idxs) > 0:
                segi0 = seq_edge_idxs[0][0]
        # check if the rightmost node connects to the breakpoint edge at index edi
        while len(seq_edge_idxs) > 0:
            segi_last = seq_edge_idxs[-1][0]
            if g.sequence_edges[segi_last].end_node != bp_node:
                del seq_edge_idxs[-1]
            else:
                break
    else:
        seq_edge_idxs = sorted(
            seq_edge_idxs,
            key=lambda item: g.sequence_edges[item[0]].start,
            reverse=True,
        )
        segi0 = seq_edge_idxs[0][0]
        if (
            len(seq_edge_idxs) > 1
            and min(g.sequence_edges[segi0].end, ref_intv.start)
            - max(g.sequence_edges[segi0].start, ref_intv.end)
            < min_overlap
        ):
            del seq_edge_idxs[0]
        segi0 = seq_edge_idxs[0][0]
        while (
            len(seq_edge_idxs) > 0 and g.sequence_edges[segi0].gap < min_overlap
        ):
            del seq_edge_idxs[0]
            if len(seq_edge_idxs) > 0:
                segi0 = seq_edge_idxs[0][0]
        # check if the rightmost node connects to the breakpoint edge at index edi
        while len(seq_edge_idxs) > 0:
            segi_last = seq_edge_idxs[-1][0]
            if g.sequence_edges[segi_last].start_node != bp_node:
                del seq_edge_idxs[-1]
            else:
                break
    if len(seq_edge_idxs) == 0:
        return []
    path_l: Walk = []
    for si, (seq_edge_idx, seq_strand) in enumerate(seq_edge_idxs):
        seq_edge = g.sequence_edges[seq_edge_idx]
        path_l.append(EdgeId("s", seq_edge_idx))
        if seq_strand == Strand.FORWARD:
            path_l.append(seq_edge.end_node)
        else:
            path_l.append(seq_edge.start_node)
        if not si >= len(seq_edge_idxs) - 1:
            continue
        next_seq_edge = g.sequence_edges[seq_edge_idxs[si + 1][0]]
        if (
            seq_strand == Strand.FORWARD
            and seq_edge.end + 1 == next_seq_edge.start
        ):
            for ci in range(len(g.concordant_edges)):
                if breakpoint_utilities.does_bp_edge_join_sequence_edges(
                    g.concordant_edges[ci],
                    seq_edge,
                    next_seq_edge,
                ):
                    path_l.append(EdgeId("c", ci))
                    path_l.append(
                        Node(
                            seq_edge.chr,
                            next_seq_edge.start,
                            Strand.REVERSE,
                        ),
                    )
                    break
        elif (
            seq_strand == Strand.REVERSE
            and seq_edge.start - 1 == next_seq_edge.end
        ):
            for ci in range(len(g.concordant_edges)):
                if breakpoint_utilities.does_bp_edge_join_sequence_edges(
                    g.concordant_edges[ci],
                    next_seq_edge,
                    seq_edge,
                ):
                    path_l.append(EdgeId("c", ci))
                    path_l.append(
                        Node(
                            seq_edge.chr,
                            next_seq_edge.end,
                            Strand.FORWARD,
                        ),
                    )
                    break
    return path_l


# TODO: unify this with path_l method
def chimeric_alignment_to_path_r(
    g: BreakpointGraph,
    ref_intv: ReferenceInterval,
    bp_node: Node,
    min_overlap: int = 500,
) -> Walk:
    """Given a breakpoint graph, a reference alignment, and a start
    node, return a traversal from the starting node to the alignment.

    g: breakpoint graph (object)
    ref_intv: end alignment interval on the reference genome
    bp_node: start node
    min_overlap: required overlap (in bp) between the alignment and the
        first/last sequence edge in the resulting path, default value is 500bp

    Returns: the resulting path as a list of alternating nodes and edges
            note that the resulting path (additionally) ends with a node
    """
    seq_edge_idxs: list[tuple[int, Strand]] = []
    for segi in range(len(g.sequence_edges)):
        sseg = g.sequence_edges[segi]
        if ref_intv.strand == Strand.FORWARD:
            if ref_intv.does_overlap(sseg.interval):
                seq_edge_idxs.append((segi, Strand.FORWARD))
        elif ref_intv.reverse.does_overlap(sseg.interval):
            seq_edge_idxs.append((segi, Strand.REVERSE))
    if len(seq_edge_idxs) == 0:
        return []
    if seq_edge_idxs[0][1] == Strand.FORWARD:
        seq_edge_idxs = sorted(
            seq_edge_idxs, key=lambda item: g.sequence_edges[item[0]].start
        )
        segi1 = seq_edge_idxs[-1][0]
        last_seq_edge = g.sequence_edges[segi1]
        if (
            min(last_seq_edge.end, ref_intv.end)
            - max(last_seq_edge.start, ref_intv.start)
            < min_overlap
        ):
            del seq_edge_idxs[-1]
        if len(seq_edge_idxs) == 0:
            return []
        segi1 = seq_edge_idxs[-1][0]
        while len(seq_edge_idxs) > 0 and g.sequence_edges[segi1].gap < 500:
            del seq_edge_idxs[-1]
            if len(seq_edge_idxs) > 0:
                segi1 = seq_edge_idxs[-1][0]
        # check if the leftmost node connects to the breakpoint edge at index edi
        while len(seq_edge_idxs) > 0:
            segi_last = seq_edge_idxs[0][0]
            if g.sequence_edges[segi_last].start_node != bp_node:
                del seq_edge_idxs[0]
            else:
                break
    else:
        seq_edge_idxs = sorted(
            seq_edge_idxs,
            key=lambda item: g.sequence_edges[item[0]].start,
            reverse=True,
        )
        segi1 = seq_edge_idxs[-1][0]
        if (
            min(g.sequence_edges[segi1].end, ref_intv.start)
            - max(g.sequence_edges[segi1].start, ref_intv.end)
            < min_overlap
        ):
            del seq_edge_idxs[-1]
        if len(seq_edge_idxs) == 0:
            return []
        segi1 = seq_edge_idxs[-1][0]
        while (
            len(seq_edge_idxs) > 0 and g.sequence_edges[segi1].gap < min_overlap
        ):
            del seq_edge_idxs[-1]
            if len(seq_edge_idxs) > 0:
                segi1 = seq_edge_idxs[-1][0]
        while len(seq_edge_idxs) > 0:
            segi_last = seq_edge_idxs[0][0]
            if g.sequence_edges[segi_last].end_node != bp_node:
                del seq_edge_idxs[0]
            else:
                break
    if len(seq_edge_idxs) == 0:
        return []
    path_r: Walk = []
    for si, (seq_edge_idx, seq_strand) in enumerate(seq_edge_idxs):
        seq_edge = g.sequence_edges[seq_edge_idx]
        if seq_edge_idxs[si][1] == Strand.FORWARD:
            path_r.append(seq_edge.start_node)
        else:
            path_r.append(seq_edge.end_node)
        path_r.append(EdgeId("s", seq_edge_idx))

        if not si >= len(seq_edge_idxs) - 1:
            continue
        next_seq_edge = g.sequence_edges[seq_edge_idxs[si + 1][0]]
        if seq_strand == Strand.FORWARD and (
            seq_edge.end + 1 == next_seq_edge.start
        ):
            for ci in range(len(g.concordant_edges)):
                if breakpoint_utilities.does_bp_edge_join_sequence_edges(
                    g.concordant_edges[ci],
                    seq_edge,
                    next_seq_edge,
                ):
                    path_r.append(seq_edge.end_node)
                    path_r.append(EdgeId("c", ci))
                    break
        if seq_strand == Strand.REVERSE and (
            seq_edge.start - 1 == next_seq_edge.end
        ):
            for ci in range(len(g.concordant_edges)):
                if breakpoint_utilities.does_bp_edge_join_sequence_edges(
                    g.concordant_edges[ci],
                    next_seq_edge,
                    seq_edge,
                ):
                    path_r.append(seq_edge.start_node)
                    path_r.append(EdgeId("c", ci))
                    break
    return path_r


def chimeric_alignment_to_path_i(
    g: BreakpointGraph,
    ref_intvs: list[ReferenceInterval],
    bp_alignment: BPIndexedAlignments,
) -> Walk:
    """Given a breakpoint graph and a list of consecutive alignments,
    return a traversal from the alignment indexed at ai1 to the alignment
    indexed at ai2, through discordant edge indexed at di.

    g: breakpoint graph (object)
    ref_intvs: alignment intervals on the reference genome
    bp_alignment: tuple of (index of the start alignment, index of the end
        alignment, index of discordant edge)

    Returns: the resulting path as a list of alternating nodes and edges
    """
    path_: Walk = [EdgeId("d", bp_alignment.discordant_idx)]
    discordant_edge = g.discordant_edges[bp_alignment.discordant_idx]
    node1 = discordant_edge.node1
    node2 = discordant_edge.node2
    if bp_alignment.alignment1 > bp_alignment.alignment2:
        path_ = (
            chimeric_alignment_to_path_l(
                g, ref_intvs[bp_alignment.alignment2], node2
            )
            + path_
            + chimeric_alignment_to_path_r(
                g, ref_intvs[bp_alignment.alignment1], node1
            )
        )
    else:
        path_ = (
            chimeric_alignment_to_path_l(
                g, ref_intvs[bp_alignment.alignment1], node1
            )
            + path_
            + chimeric_alignment_to_path_r(
                g, ref_intvs[bp_alignment.alignment2], node2
            )
        )
    return path_


def traverse_through_sequence_edge(
    g: BreakpointGraph, start_node: Node, end_node: Node
) -> Walk:
    """Given a breakpoint graph and two nodes, return a traversal through sequence and concordant edges between the two nodes

    g: breakpoint graph (object)
    start_node: start node
    end_node: end node - start and end node must locate at different (left/right) ends on the corresponding sequence edges

    Returns: the resulting path as a list of alternating nodes and edges
            note that the resulting path (additionally) starts and ends with the given nodes
    """
    assert start_node.strand != end_node.strand
    path_: Walk = [start_node]
    seqi = g.node_adjacencies[start_node].sequence[0]
    seq_edge = g.sequence_edges[seqi]
    next_end = seq_edge.start_node
    if start_node.strand == Strand.REVERSE:
        next_end = seq_edge.end_node
    path_.append(EdgeId("s", seqi))
    path_.append(next_end)
    while next_end != end_node:
        # ignore the alignments spanning two amplicon intervals
        if not (conc_edges := g.node_adjacencies[next_end].concordant):
            return path_

        ci = conc_edges[0]
        path_.append(EdgeId("c", ci))
        cedge = g.concordant_edges[ci]
        next_start = cedge.node1
        if next_start == next_end:
            next_start = cedge.node2
        path_.append(next_start)

        seqi = g.node_adjacencies[next_start].sequence[0]
        seq_edge = g.sequence_edges[seqi]
        next_end = seq_edge.start_node
        if next_start.strand == Strand.REVERSE:
            next_end = seq_edge.end_node
        path_.append(EdgeId("s", seqi))
        path_.append(next_end)
    return path_


def chimeric_alignment_to_path(
    g: BreakpointGraph,
    ref_intvs: list[ReferenceInterval],
    ai_list: list[tuple[int, int]],
    bp_list: list[int],
) -> Walk:
    """Convert chimeric alignments to path"""
    path_: Walk = []
    lastnode: Node
    for i in range(len(bp_list)):
        di = bp_list[i]
        disc_edge = g.discordant_edges[di]
        node1 = disc_edge.node1
        node2 = disc_edge.node2

        if ai_list[i][0] > ai_list[i][1]:
            if i == 0:
                path_ = chimeric_alignment_to_path_l(
                    g, ref_intvs[ai_list[i][1]], node2
                ) + [
                    EdgeId("d", bp_list[i]),
                ]
                lastnode = node1
            else:
                path_ += traverse_through_sequence_edge(g, lastnode, node2)
                path_.append(EdgeId("d", bp_list[i]))
                lastnode = node1
                if i == len(bp_list) - 1:
                    path_ += chimeric_alignment_to_path_r(
                        g, ref_intvs[ai_list[i][0]], node1
                    )
        elif i == 0:
            path_ = chimeric_alignment_to_path_l(
                g, ref_intvs[ai_list[i][0]], node1
            ) + [
                EdgeId("d", bp_list[i]),
            ]
            lastnode = node2
        else:
            path_ += traverse_through_sequence_edge(g, lastnode, node1)
            path_.append(EdgeId("d", bp_list[i]))
            lastnode = node2
            if i == len(bp_list) - 1:
                path_ += chimeric_alignment_to_path_r(
                    g, ref_intvs[ai_list[i][1]], node2
                )
    return path_


def longest_path_dict(
    path_constraints_: list[PathConstraint],
) -> list[FinalizedPathConstraint]:
    """Convert paths from a list of alternating nodes and edges into a dict of
    edges. Only keep the longest paths, i.e., those which are not a subpath of
    any other path.
    """
    final_paths: list[FinalizedPathConstraint] = []
    for path_i in range(len(path_constraints_)):
        path = path_constraints_[path_i].path
        edge_counts: dict[EdgeId, int] = defaultdict(int)
        for edge_id in range(len(path)):
            # Only track edge counts, not nodes
            if edge_id % 2 == 0:
                edge_counts[path[edge_id]] += 1  # type: ignore[index]
        final_paths.append(
            FinalizedPathConstraint(
                edge_counts=edge_counts,
                pc_idx=path_i,
                support=path_constraints_[path_i].support,
            )
        )
    for path_i, final_pc in enumerate(reversed(final_paths)):
        subpath_idx: int | None = None
        path_edge_counts = final_pc.edge_counts
        for subpath_i, subpath_constraint in enumerate(final_paths):
            subpath_edge_counts = subpath_constraint.edge_counts
            is_subpath = True
            for edge in path_edge_counts:
                if (
                    edge not in subpath_edge_counts
                    or subpath_edge_counts[edge] < path_edge_counts[edge]
                ):
                    is_subpath = False
                    break
            if is_subpath and subpath_i != path_i:
                subpath_idx = subpath_i
                break
        if subpath_idx is not None:
            del final_paths[path_i]
            del final_paths[path_i]
            subpath_constraint.support = max(
                subpath_constraint.support, final_paths[path_i].support
            )
            del final_paths[path_i]
    return final_paths
