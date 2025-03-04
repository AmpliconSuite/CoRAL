"""
Utilities for scoring reconstructions against simulated amplicons.
"""

from __future__ import annotations

import sys
import warnings
from itertools import combinations, product
from typing import Optional

from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import Node, SequenceEdge, Strand

warnings.simplefilter(action="ignore", category=FutureWarning)

import networkx as nx
import numpy as np
import pandas as pd
import pyranges
import scipy

from coral.scoring import io_utils, scoring_utils


def is_breakpoint_match(
    node1: Node,
    node2: Node,
    other_node1: Node,
    other_node2: Node,
    tolerance: int = 100,
) -> bool:
    # Need to check both orientations
    return is_oriented_breakpoint_match(
        node1, node2, other_node1, other_node2, tolerance
    ) or is_oriented_breakpoint_match(
        node1, node2, other_node2, other_node1, tolerance
    )


def is_oriented_breakpoint_match(
    node1: Node,
    node2: Node,
    other_node1: Node,
    other_node2: Node,
    tolerance: int = 100,
) -> bool:
    if node1.chr != other_node1.chr or node2.chr != other_node2.chr:
        return False

    if node1.strand != other_node1.strand or node2.strand != other_node2.strand:
        return False

    return (
        abs(node1.pos - other_node1.pos) < tolerance
        and abs(node2.pos - other_node2.pos) < tolerance
    )


def repartition_intervals(cycle, intervals):
    new_intervals = []
    for _, interval in cycle.iterrows():
        segments = pyranges.PyRanges(intervals).intersect(
            pyranges.PyRanges(pd.DataFrame(interval).T)
        )

        if len(segments) > 0 and interval.orientation == "-":
            segments_df = segments.df.iloc[::-1].copy()
        else:
            segments_df = segments.df.copy()

        orientation = interval.orientation
        for _, segment in segments_df.iterrows():
            new_intervals.append(
                [
                    segment.Chromosome,
                    segment.Start,
                    segment.End,
                    orientation,
                    segment.segment_label,
                ]
            )

    return pd.DataFrame(
        new_intervals,
        columns=["Chromosome", "Start", "End", "orientation", "segment_label"],
    )


def find_longest_subsequence(sequence1, sequence2, segment_lengths):
    m = len(sequence1)
    n = len(sequence2)

    # instantiate DP array
    L = [[None] * (n + 1) for i in range(m + 1)]

    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif sequence1[i - 1] == sequence2[j - 1]:
                L[i][j] = L[i - 1][j - 1] + segment_lengths[sequence1[i - 1][0]]
            else:
                L[i][j] = max(L[i - 1][j], L[i][j - 1])

    return L[m][n]


def compute_normalized_longest_cycle_subsequence(
    true_cycle: pd.DataFrame,
    reconstructed_cycles: pd.DataFrame,
    binned_genome: pyranges.PyRanges,
) -> tuple[float, float]:
    flip_orientation = lambda x: "+" if x == "-" else "-"

    overall_max_lcs, overall_max_normalized_lcs = 0, 0

    for cycle in reconstructed_cycles["cycle_id"].unique():
        _true_cycle = true_cycle.copy()
        _true_cycle["segment"] = _true_cycle.apply(
            lambda x: (x.Chromosome, x.Start, x.End), axis=1
        )

        reconstructed_cycle = reconstructed_cycles[
            reconstructed_cycles["cycle_id"] == cycle
        ].copy()
        reconstructed_cycle["segment"] = reconstructed_cycle.apply(
            lambda x: (x.Chromosome, x.Start, x.End), axis=1
        )

        binned_genome, _ = io_utils.bin_genome(_true_cycle, reconstructed_cycle)

        grs = {
            n: s
            for n, s in zip(
                ["bins", "e1", "e2"],
                [
                    pyranges.PyRanges(binned_genome),
                    pyranges.PyRanges(_true_cycle),
                    pyranges.PyRanges(reconstructed_cycle),
                ],
            )
        }

        overlaps = pyranges.count_overlaps(grs)
        overlaps = overlaps.as_df()
        overlaps["interval_length"] = abs(overlaps["End"] - overlaps["Start"])
        overlaps["total_interval_length"] = overlaps.apply(
            lambda x: max(x.e1, x.e2) * x.interval_length, axis=1
        )
        total_length = overlaps["total_interval_length"].sum()

        # find bins that are covered by both reconstruction and true
        overlaps = overlaps[(overlaps["e1"] > 0) & (overlaps["e2"] > 0)]
        overlaps["segment_label"] = [i + 1 for i in range(len(overlaps))]

        if len(overlaps) == 0:
            continue

        # repartition cycles with respect to overlaps
        true_cycle_repartitioned = repartition_intervals(_true_cycle, overlaps)
        reconstructed_cycle_repartitioned = repartition_intervals(
            reconstructed_cycle, overlaps
        )

        # _true_cycle = _true_cycle[_true_cycle['segment_label'] != 0]
        # reconstructed_cycle = reconstructed_cycle[reconstructed_cycle['segment_label'] != 0]

        true_cycle_string = true_cycle_repartitioned.apply(
            lambda x: (x.segment_label, x.orientation), axis=1
        ).values
        reconstructed_cycle_string_fw = reconstructed_cycle_repartitioned.apply(
            lambda x: (x.segment_label, x.orientation), axis=1
        ).values
        reconstructed_cycle_string_rev = [
            (x[0], flip_orientation(x[1]))
            for x in reconstructed_cycle_string_fw[::-1]
        ]

        max_lcs = 0

        for i in range(len(reconstructed_cycle_string_fw)):
            rotated_reconstruction = (
                list(reconstructed_cycle_string_fw)[i:]
                + list(reconstructed_cycle_string_fw)[:i]
            )
            common_subsequence_length = find_longest_subsequence(
                true_cycle_string,
                rotated_reconstruction,
                overlaps[["segment_label", "interval_length"]]
                .set_index("segment_label")["interval_length"]
                .to_dict(),
            )
            if common_subsequence_length > max_lcs:
                max_lcs = common_subsequence_length

        for i in range(len(reconstructed_cycle_string_rev)):
            rotated_reconstruction = (
                list(reconstructed_cycle_string_rev)[i:]
                + list(reconstructed_cycle_string_rev)[:i]
            )
            common_subsequence_length = find_longest_subsequence(
                true_cycle_string,
                rotated_reconstruction,
                overlaps[["segment_label", "interval_length"]]
                .set_index("segment_label")["interval_length"]
                .to_dict(),
            )
            if common_subsequence_length > max_lcs:
                max_lcs = common_subsequence_length

        if (max_lcs / total_length) > overall_max_normalized_lcs:
            overall_max_lcs, overall_max_normalized_lcs = (
                max_lcs,
                (max_lcs / total_length),
            )

    return overall_max_lcs, overall_max_normalized_lcs


def get_fragments_similarity_unweighted(
    true_fragments: pyranges.PyRanges,
    reconstructed_fragments: pyranges.PyRanges,
    bins: pyranges.PyRanges,
) -> tuple[float, float]:
    """
    Compute the overlap distance between the fragment counts of two reconstructions which overlap every genomic bin.

    Args:
        e1 (pd.DataFrame): First reconstruction
        e2 (pd.DataFrame): Second reconstruction
        bins (pd.DataFrame): Bins of the intervals union e1 and e2
    """

    # find minimum distance in cycles
    overlap_distance = np.inf
    length_ratio = 0

    for cycle in reconstructed_fragments.df["cycle_id"].unique():
        reconstructed_cycle = reconstructed_fragments[
            reconstructed_fragments.df["cycle_id"] == cycle
        ]
        grs = {
            n: s
            for n, s in zip(
                ["bins", "e1", "e2"],
                [bins, true_fragments, reconstructed_cycle],
            )
        }

        overlaps = pyranges.count_overlaps(grs)

        overlaps = overlaps.as_df()
        overlaps["overlapping_score"] = overlaps.apply(
            lambda x: abs(x.e1 - x.e2), axis=1
        )
        overlaps["len"] = abs(overlaps["End"] - overlaps["Start"])
        # account for duplicated fragments
        overlaps["totallen"] = overlaps.apply(
            lambda x: max(x.e1, x.e2) * x.len, axis=1
        )
        overlaps["prod"] = overlaps["overlapping_score"] * overlaps["len"]

        true_length = overlaps.apply(lambda x: x.e1 * x.len, axis=1).sum()
        reconstructed_length = overlaps.apply(
            lambda x: x.e2 * x.len, axis=1
        ).sum()

        _overlap_distance = (overlaps["prod"].sum()) / (
            overlaps["totallen"].sum()
        )

        if _overlap_distance < overlap_distance:
            overlap_distance = _overlap_distance
            length_ratio = (reconstructed_length - true_length) / true_length

    return (1.0 - overlap_distance), length_ratio


def find_overlapping_sequence_edges(
    true_graph: BreakpointGraph,
    reconstructed_graph: BreakpointGraph,
    tolerance: int = 1000,
) -> int:
    # naively test if we find an edge that looks similar
    # for all edges in edges1
    num_edges_found = 0
    for edge in true_graph.sequence_edges:
        node1 = Node(edge.chr, edge.start, Strand.REVERSE)
        node2 = Node(edge.chr, edge.end, Strand.FORWARD)

        for other_edge in reconstructed_graph.sequence_edges:
            other_node1 = Node(other_edge.chr, other_edge.start, Strand.REVERSE)
            other_node2 = Node(other_edge.chr, other_edge.end, Strand.FORWARD)
            if is_breakpoint_match(
                node1, node2, other_node1, other_node2, tolerance
            ):
                num_edges_found += 1
                break

    return num_edges_found  # , num_edges_found / len(true_graph.sequence_edges)


def find_overlapping_bp_edges(
    true_graph: BreakpointGraph,
    reconstructed_graph: BreakpointGraph,
    tolerance: int = 1000,
) -> int:
    # naively test if we find an edge that looks similar
    # for all edges in edges1
    num_edges_found = 0
    for edge in true_graph.discordant_edges:
        node1 = Node(edge.node1.chr, edge.node1.pos, edge.node1.strand)
        node2 = Node(edge.node2.chr, edge.node2.pos, edge.node2.strand)

        for other_edge in reconstructed_graph.discordant_edges:
            other_node1 = Node(
                other_edge.node1.chr,
                other_edge.node1.pos,
                other_edge.node1.strand,
            )
            other_node2 = Node(
                other_edge.node2.chr,
                other_edge.node2.pos,
                other_edge.node2.strand,
            )
            if is_breakpoint_match(
                node1, node2, other_node1, other_node2, tolerance
            ):
                num_edges_found += 1
                break

    return (
        num_edges_found  # , num_edges_found / len(true_graph.breakpoint_edges)
    )


def create_cycle_graph(cycle):
    cycle_graph_fw = nx.DiGraph()
    cycle_graph_rev = nx.DiGraph()

    flip_orientation = lambda x: "-" if x == "+" else "+"

    for segment_i in range(len(cycle) - 1):
        segment_1, segment_2 = (
            cycle.loc[segment_i, ["Chromosome", "Start", "End", "orientation"]],
            cycle.loc[
                segment_i + 1, ["Chromosome", "Start", "End", "orientation"]
            ],
        )

        cycle_graph_fw.add_node(segment_i, label=tuple(segment_1.values))
        cycle_graph_fw.add_node(segment_i + 1, label=tuple(segment_2.values))

        cycle_graph_fw.add_edge(segment_i, segment_i + 1)

        rev_segment_1 = segment_1.values
        rev_segment_2 = segment_2.values

        rev_segment_1[3] = flip_orientation(rev_segment_1[3])
        rev_segment_2[3] = flip_orientation(rev_segment_2[3])

        cycle_graph_rev.add_node(segment_i, label=tuple(rev_segment_1))
        cycle_graph_rev.add_node(segment_i + 1, label=tuple(rev_segment_2))

        cycle_graph_rev.add_edge(segment_i + 1, segment_i)

    last_segment = cycle.loc[
        len(cycle) - 1, ["Chromosome", "Start", "End", "orientation"]
    ]

    cycle_graph_fw.add_node(len(cycle) - 1, label=tuple(last_segment.values))
    cycle_graph_fw.add_edge(len(cycle) - 1, 0)

    rev_last_segment = last_segment.values
    rev_last_segment[3] = flip_orientation(rev_last_segment[3])

    cycle_graph_rev.add_node(len(cycle) - 1, label=tuple(rev_last_segment))
    cycle_graph_rev.add_edge(
        0,
        len(cycle) - 1,
    )

    return cycle_graph_fw, cycle_graph_rev


def match_segments(
    true_cycle, true_graph, reconstructed_cycle, reconstructed_graph
):
    matching_graph = nx.Graph()
    for segment_i in range(len(reconstructed_cycle)):
        r_chr, r_start, r_end, r_strand = reconstructed_cycle.loc[
            segment_i, ["Chromosome", "Start", "End", "orientation"]
        ].values

        for segment_j in range(len(true_cycle)):
            weight = 3e9
            true_chr, true_start, true_end, true_strand = true_cycle.loc[
                segment_j, ["Chromosome", "Start", "End", "orientation"]
            ].values

            if true_chr == r_chr:
                weight = abs(true_start - r_start) + abs(true_end - r_end)

            if (
                true_chr,
                true_start,
                true_end,
                true_strand,
            ) not in reconstructed_graph.nodes:
                matching_graph.add_edge(
                    ("true", segment_j),
                    segment_i,
                    weight=weight,
                )

    segment_mapping = nx.bipartite.minimum_weight_full_matching(matching_graph)

    remapped_cycle_graph = nx.relabel_nodes(
        reconstructed_graph,
        segment_mapping,
        copy=True,
    )

    # remap back by removing the identifier 'true'
    clean_segment_mapping = {
        n: n[1] if type(n) == tuple else n for n in remapped_cycle_graph.nodes
    }
    clean_remapped_cycle_graph = nx.relabel_nodes(
        remapped_cycle_graph,
        clean_segment_mapping,
        copy=True,
    )

    for n in clean_remapped_cycle_graph:
        if n not in true_graph.nodes:
            continue

        true_graph_label = true_graph.nodes(data=True)[n]["label"]
        reconstructed_label = clean_remapped_cycle_graph.nodes(data=True)[n][
            "label"
        ]
        relabeled = (
            true_graph_label[0],
            true_graph_label[1],
            true_graph_label[2],
            reconstructed_label[3],
        )  # preserve orientation

        clean_remapped_cycle_graph.nodes(data=True)[n]["label"] = relabeled

    return clean_remapped_cycle_graph


def rank_triplet(triplet, true_path, reconstructed_path):
    def find_index(lst, value):
        return [i for i in range(len(lst)) if lst[i] == value]

    triplet_i, triplet_j, triplet_k = [triplet[i] for i in range(3)]

    if (
        triplet_i in reconstructed_path
        and triplet_j in reconstructed_path
        and triplet_k in reconstructed_path
    ):
        true_i, true_j, true_k = (
            find_index(true_path, triplet_i),
            find_index(true_path, triplet_j),
            find_index(true_path, triplet_k),
        )
        reconstructed_i, reconstructed_j, reconstructed_k = (
            find_index(reconstructed_path, triplet_i),
            find_index(reconstructed_path, triplet_j),
            find_index(reconstructed_path, triplet_k),
        )

        all_acceptable_orderings = product(*[true_i, true_j, true_k])
        all_reconstructed_orderings = product(
            *[reconstructed_i, reconstructed_j, reconstructed_k]
        )

        for ti, tj, tk in all_acceptable_orderings:
            for ri, rj, rk in all_reconstructed_orderings:
                true_ranks = tuple(scipy.stats.rankdata([ti, tj, tk]))

                reconstructed_ranks = tuple(scipy.stats.rankdata([ri, rj, rk]))

                if true_ranks == reconstructed_ranks:
                    return 1

    return 0


def score_triplets_correct(true_cycle, reconstructed_cycles):
    def triplet_not_in_path(path, triplet):
        return np.any([triplet[i] not in path for i in range(len(triplet))])

    true_cycle_graph, _ = create_cycle_graph(true_cycle)

    best_triplets_correct = np.nan

    for cycle_id in reconstructed_cycles["cycle_id"].unique():
        reconstructed_cycle = reconstructed_cycles[
            reconstructed_cycles["cycle_id"] == cycle_id
        ].reset_index()

        (
            reconstructed_cycle_graph_fw,
            reconstructed_cycle_graph_rev,
        ) = create_cycle_graph(reconstructed_cycle)

        remapped_reconstructed_cycle_graph_fw = match_segments(
            true_cycle,
            true_cycle_graph,
            reconstructed_cycle,
            reconstructed_cycle_graph_fw,
        )
        remapped_reconstructed_cycle_graph_rev = match_segments(
            true_cycle,
            true_cycle_graph,
            reconstructed_cycle,
            reconstructed_cycle_graph_rev,
        )

        # find best source possible
        source = np.intersect1d(
            list(true_cycle_graph.nodes),
            list(remapped_reconstructed_cycle_graph_fw.nodes),
        )
        if len(source) == 0:
            source = np.intersect1d(
                list(true_cycle_graph.nodes),
                list(remapped_reconstructed_cycle_graph_rev.nodes),
            )
            if len(source) == 0:
                continue
        source = source[0]

        # find paths
        true_path = []
        for node in nx.depth_first_search.dfs_preorder_nodes(
            true_cycle_graph, source=source
        ):
            true_path.append(true_cycle_graph.nodes(data=True)[node]["label"])

        reconstructed_fw_path = []
        if source in remapped_reconstructed_cycle_graph_fw.nodes:
            for node in nx.depth_first_search.dfs_preorder_nodes(
                remapped_reconstructed_cycle_graph_fw, source=source
            ):
                reconstructed_fw_path.append(
                    remapped_reconstructed_cycle_graph_fw.nodes(data=True)[
                        node
                    ]["label"]
                )

        reconstructed_rev_path = []
        if source in remapped_reconstructed_cycle_graph_rev.nodes:
            for node in nx.depth_first_search.dfs_preorder_nodes(
                remapped_reconstructed_cycle_graph_rev, source=source
            ):
                reconstructed_rev_path.append(
                    remapped_reconstructed_cycle_graph_rev.nodes(data=True)[
                        node
                    ]["label"]
                )

        triplets_correct = 0
        num_triplets = 0

        all_triplets = combinations(true_path, 3)
        for triplet_i, triplet_j, triplet_k in set(all_triplets):
            if (
                triplet_i == triplet_j
                or triplet_i == triplet_k
                or triplet_j == triplet_k
            ):
                continue

            num_triplets += 1

            if triplet_not_in_path(
                reconstructed_fw_path, (triplet_i, triplet_j, triplet_k)
            ) and triplet_not_in_path(
                reconstructed_rev_path, (triplet_i, triplet_j, triplet_k)
            ):
                continue

            triplet_score_fw = rank_triplet(
                (triplet_i, triplet_j, triplet_k),
                true_path,
                reconstructed_fw_path,
            )
            triplet_score_rev = rank_triplet(
                (triplet_i, triplet_j, triplet_k),
                true_path,
                reconstructed_rev_path,
            )

            triplets_correct += max(triplet_score_fw, triplet_score_rev)

        if np.isnan(best_triplets_correct) and triplets_correct > 0:
            best_triplets_correct = triplets_correct / max(1, num_triplets)

        elif (triplets_correct / max(1, num_triplets)) > best_triplets_correct:
            best_triplets_correct = triplets_correct / max(1, num_triplets)

    return best_triplets_correct


def get_seq_edge_df(sequence_edges: list[SequenceEdge]) -> pd.DataFrame:
    seq_edge_df = scoring_utils.get_seq_edges_from_decoil(
        sequence_edges,
        base_coverage=13.0,
    )
    return seq_edge_df


def get_cycle_copy_number_ratio(
    cycles_bed_df: pd.DataFrame, sequence_edges: list[SequenceEdge], n: int = 1
) -> float:
    """Computes the ratio of the best cycle weighted copy-number to the total."""
    # sequence_edges["CopyCount"] = sequence_edges["CopyCount"].astype(float)
    # # sequence_edges['ident'] = sequence_edges.apply(lambda x: (x.Chrom1, x.Start, x.Chrom2, x.End), axis=1)
    # # to_drop = sequence_edges.loc[sequence_edges['CopyCount'] < 5.0, 'ident'].values

    # # sequence_edges = sequence_edges[~sequence_edges['ident'].isin(to_drop)]

    # sequence_edges = sequence_edges[sequence_edges["CopyCount"] >= 5.0]

    seq_edges = [edge for edge in sequence_edges if edge.cn >= 5.0]
    total_weight = sum(
        [edge.cn * abs(edge.end - edge.start) for edge in seq_edges]
    )

    # total_weight = sequence_edges.apply(
    #     lambda x: x.CopyCount * abs(x.End - x.Start), axis=1
    # ).sum()

    cycle_to_ratio = []
    for _, cycle in cycles_bed_df.groupby("cycle_id"):
        # cycle['ident'] = cycle.apply(lambda x: (x.Chromosome, x.Start, x.Chromosome, x.End), axis=1)
        # cycle = cycle[~cycle['ident'].isin(to_drop)]

        if len(cycle) == 0:
            continue

        if np.any(cycle["iscyclic"] == False):
            continue

        cycle_weight = (
            cycle.apply(lambda x: x.weight * abs(x.End - x.Start), axis=1).sum()
        ) / total_weight
        cycle_to_ratio.append(cycle_weight)

    if n >= len(cycle_to_ratio):
        return 1.0
    return np.sum(np.partition(cycle_to_ratio, -n)[-n:])
