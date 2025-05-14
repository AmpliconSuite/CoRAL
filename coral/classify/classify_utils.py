from __future__ import annotations

import importlib.resources
from collections import defaultdict
from typing import TYPE_CHECKING

from intervaltree import IntervalTree

from coral import constants, datatypes, supplemental_data
from coral.breakpoint import breakpoint_utils

if TYPE_CHECKING:
    from coral.breakpoint.breakpoint_graph import BreakpointGraph


def find_closest_true_graphs(
    true_graphs: dict[str, BreakpointGraph],
    recon_graphs: dict[str, BreakpointGraph],
) -> dict[str, str | None]:
    # Build a dict of the form {true_graph_id: recon_graph_id}
    # If no matching true graph is found, the value is None
    matching_true_graphs = defaultdict(None)
    for true_id, true_graph in true_graphs.items():
        max_overlapping_bp_edges = 0
        for recon_id, recon_graph in recon_graphs.items():
            num_overlapping_bp_edges = (
                breakpoint_utils.find_overlapping_bp_edges(
                    true_graph, recon_graph
                )
            )
            if num_overlapping_bp_edges > max_overlapping_bp_edges:
                max_overlapping_bp_edges = num_overlapping_bp_edges
                matching_true_graphs[true_id] = recon_id
    return matching_true_graphs


def build_low_complexity_region_map() -> dict[str, IntervalTree]:
    regions_file = (
        importlib.resources.files(supplemental_data)
        / "exclude.cnvnator_100bp.GRCh38.20170403.bed"
    )
    regions_dict = defaultdict(IntervalTree)
    with regions_file.open("r") as infile:
        _ = next(infile)
        for line in infile:
            if not (fields := line.rstrip().rsplit(sep="\t")):
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            # Ignore decoy contigs
            if chrom not in constants.CHR_TAG_TO_IDX:
                continue
            # Ignore regions that are too short
            if not (end - start > 7500):
                continue
            regions_dict[chrom].addi(start, end)
    return regions_dict


def get_sequence_edge_properties(
    bp_graph: BreakpointGraph,
    low_complexity_region_map: dict[str, IntervalTree],
    min_cn_gain: float = 4.5,
) -> tuple[float, int]:
    """
    Returns the max CN over seq edges, as well as the total size of amplified
    intervals with CN over the given gain threshold.
    """
    max_cn = 0
    tot_amplified_seq_len = 0
    tot_weight = 0.0
    for edge in bp_graph.sequence_edges:
        if low_complexity_region_map[edge.chr].overlaps(edge.start, edge.end):
            continue
        if len(edge) < 1000:
            continue
        tot_weight += len(edge) * edge.cn
        max_cn = max(max_cn, edge.cn)
        if edge.cn > min_cn_gain:
            tot_amplified_seq_len += len(edge)
    return max_cn, tot_amplified_seq_len


def get_discordant_edge_properties(
    bp_graph: BreakpointGraph,
    low_complexity_region_map: dict[str, IntervalTree],
    foldback_dist_threshold: int = 25_000,
) -> tuple[int, int, float]:
    """
    Returns the number of total rearranged edges, the number of foldback edges,
    and the fraction of reads that support foldback edges.
    """
    fb_reads, non_fb_reads, fb_edges, rearranged_edges = 0, 0, 0, 0
    for edge in bp_graph.discordant_edges:
        if (
            edge.node1.pos in low_complexity_region_map[edge.node1.chr]
            or edge.node2.pos in low_complexity_region_map[edge.node2.chr]
        ):
            continue
        if (
            edge.node1.chr == edge.node2.chr
            and abs(edge.node2.pos - edge.node1.pos)
            <= 2000  # TODO: check with Jens if this should also be foldback_dist_threshold, isntead of matching 2000 from `compute_f_from_AA_graph`
            and edge.node1.strand == datatypes.Strand.FORWARD
            and edge.node2.strand == datatypes.Strand.REVERSE
        ):
            continue
        if edge.node1.strand != edge.node2.strand:
            non_fb_reads += edge.lr_count
            if abs(edge.node2.pos - edge.node1.pos) > foldback_dist_threshold:
                rearranged_edges += 1
            continue
        rearranged_edges += 1
        if (
            edge.node1.chr != edge.node2.chr
            or (edge.node2.pos - edge.node1.pos) > foldback_dist_threshold
        ):
            non_fb_reads += edge.lr_count
            continue
        # Satisfies all conditions for a foldback edge
        fb_reads += edge.lr_count
        fb_edges += 1
    return (
        rearranged_edges,
        fb_edges,
        fb_reads / max(1, fb_reads + non_fb_reads),
    )
