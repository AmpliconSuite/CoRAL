"""
Utilities for breakpoint graph inference.
"""

from __future__ import annotations

import io
import logging
from collections import Counter, defaultdict
from typing import (
    TYPE_CHECKING,
    Generator,
    List,
    Optional,
    Sequence,
    Tuple,
    cast,
)

import numpy as np
import pysam

from coral import bam_types, cigar_parsing, core_types, core_utils, datatypes
from coral.constants import CHR_TAG_TO_IDX
from coral.datatypes import (
    AmpliconInterval,
    BPAlignments,
    Breakpoint,
    BreakpointStats,
    ChimericAlignment,
    ChrPairOrientation,
    CNSInterval,
    Interval,
    Node,
    Strand,
)

if TYPE_CHECKING:
    from coral.breakpoint.breakpoint_graph import BreakpointGraph

logger = logging.getLogger(__name__)


def interval_overlap(int1: list[int], int2: list[int]) -> bool:
    """Check if two chromosome intervals overlap (share a subsequence).

    Intervals are given in the form [chr, start, end], where:
        chr: chromosome number
        s: start position/index
        e: end position/index
    """
    return (
        int1[0] == int2[0]
        and int(int1[1]) <= int(int2[2])
        and int(int2[1]) <= int(int1[2])
    )


def interval_include(int1: list[int], int2: list[int]) -> bool:
    """
    Check if an interval in the form of [chr, s, e] is fully included in another
    """
    return (
        int1[0] == int2[0]
        and int(int1[1]) >= int(int2[1])
        and int(int1[2]) <= int(int2[2])
    )


def interval_adjacent(int1: list[int], int2: list[int]) -> bool:
    """
    Check if two intervals in the form of [chr, s, e] are adjacent
    """
    if int1[0] != int2[0]:
        return False
    if int1[1] <= int2[1]:
        return int2[1] == int1[2] + 1
    return int1[1] == int2[2] + 1


def interval_overlap_l(
    intv: Interval, other_intvs: Sequence[Interval]
) -> int | None:
    """
    Check if an interval in the form of [chr, s, e] overlaps with any member of
    a list of intervals.
    """
    for other_intv_idx, other_intv in enumerate(other_intvs):
        if intv.does_overlap(other_intv):
            return other_intv_idx
    return None


def interval_include_l(int1: list[int], intl: list[list[int]]) -> int:
    for int2i in range(len(intl)):
        if interval_include(int1, intl[int2i]):
            return int2i
    return -1


def interval_exclusive(
    int1: AmpliconInterval, intl: list[AmpliconInterval]
) -> tuple[set[int], list[AmpliconInterval]]:
    overlap_ints = set([])
    intl_ = [int1]
    for int2i in range(len(intl)):
        for inti_ in range(len(intl_))[::-1]:
            int_ = intl_[inti_]
            if int_.does_overlap(intl[int2i]):
                overlap_ints.add(int2i)
                del intl_[inti_]
                if int_.start < intl[int2i].start:
                    intl_.append(
                        AmpliconInterval(
                            int_.chr, int_.start, intl[int2i].start - 1
                        )
                    )
                if int_.end > intl[int2i].end:
                    intl_.append(
                        AmpliconInterval(
                            int_.chr, intl[int2i].end + 1, int_.end, -1
                        )
                    )
    return overlap_ints, intl_


def alignment2bp(
    rn: str,
    alignments: list[datatypes.ChimericAlignment],
    min_bp_match_cutoff: int,
    min_mapq: float,
    new_intv: Interval,
    orig_intv: Interval,
    gap_mapq: float = 10,
    max_nm: float | None = None,
):
    bp_list: list[Breakpoint] = []
    has_bp_assigned = [False] * (len(alignments) - 1)
    alignments.sort(key=lambda x: x.query_bounds.start)

    # Breakpoint from local alignment i and i + 1
    for ri in range(len(alignments) - 1):
        chimera1, chimera2 = alignments[ri], alignments[ri + 1]
        query_gap = chimera2.query_bounds.start - chimera1.query_bounds.end
        if query_gap + min_bp_match_cutoff < 0:
            continue
        if chimera1.mapq < min_mapq or chimera2.mapq < min_mapq:
            continue
        if max_nm is not None and (
            chimera1.edit_dist >= max_nm or chimera2.edit_dist >= max_nm
        ):
            continue

        ref_intv1, ref_intv2 = chimera1.ref_interval, chimera2.ref_interval
        if (
            ref_intv1.does_overlap(new_intv)
            and ref_intv2.does_overlap(orig_intv)
        ) or (
            ref_intv2.does_overlap(new_intv)
            and ref_intv1.does_overlap(orig_intv)
        ):
            bp_list.append(
                interval2bp(
                    chimera1,
                    chimera2,
                    BPAlignments(rn, ri, ri + 1),
                    query_gap,
                )
            )
            has_bp_assigned[ri] = True

    # Breakpoint from local alignment i - 1 and i + 1
    for ri in range(1, len(alignments) - 1):
        prev_chimera, curr_chimera, next_chimera = (
            alignments[ri - 1],
            alignments[ri],
            alignments[ri + 1],
        )
        if has_bp_assigned[ri - 1] or has_bp_assigned[ri]:
            continue
        if (
            curr_chimera.mapq >= gap_mapq
            or prev_chimera.mapq < min_mapq
            or next_chimera.mapq < min_mapq
        ):
            continue
        if max_nm is not None and (
            prev_chimera.edit_dist >= max_nm or next_chimera.edit_dist >= max_nm
        ):
            continue
        prev_intv, next_intv = (
            prev_chimera.ref_interval,
            next_chimera.ref_interval,
        )
        if (
            prev_intv.does_overlap(new_intv)
            and next_intv.does_overlap(orig_intv)
        ) or (
            next_intv.does_overlap(new_intv)
            and prev_intv.does_overlap(orig_intv)
        ):
            query_gap = (
                next_chimera.query_bounds.start - prev_chimera.query_bounds.end
            )

            bp_list.append(
                interval2bp(
                    prev_chimera,
                    next_chimera,
                    BPAlignments(rn, ri - 1, ri + 1),
                    query_gap,
                )
            )
    return bp_list


def alignment2bp_l(
    rn: str,
    alignments: list[datatypes.ChimericAlignment],
    min_bp_match_cutoff: int,
    min_mapq: int,
    gap_: int,
    intrvls: Sequence[Interval],
    gap_mapq: int = 10,
    max_nm: float | None = None,
) -> list[Breakpoint]:
    bp_list: list[Breakpoint] = []
    has_bp_assigned = [False] * (len(alignments) - 1)
    alignments.sort(key=lambda x: x.query_bounds.start)

    """
	Breakpoint from local alignment i and i + 1
	"""

    for i in range(len(alignments) - 1):
        """
		Add unmatched breakpoint to new_bp_list
		"""
        chimera1, chimera2 = alignments[i], alignments[i + 1]
        query_gap = chimera2.query_bounds.start - chimera1.query_bounds.end
        if query_gap + min_bp_match_cutoff < 0:
            continue
        if chimera1.mapq < min_mapq or chimera2.mapq < min_mapq:
            continue
        if max_nm is not None and (
            chimera1.edit_dist >= max_nm or chimera2.edit_dist >= max_nm
        ):
            continue
        if (io1 := interval_overlap_l(chimera1.ref_interval, intrvls)) is None:
            continue
        if (io2 := interval_overlap_l(chimera2.ref_interval, intrvls)) is None:
            continue
        if io1 != io2:
            continue
        ref_intv1, ref_intv2 = chimera1.ref_interval, chimera2.ref_interval

        ref_gap = abs(ref_intv2.start - ref_intv1.end)
        if ref_intv1.strand != ref_intv2.strand or (
            abs(ref_gap - query_gap) > max(gap_, abs(query_gap * 0.2))
        ):
            bp_list.append(
                interval2bp(
                    chimera1, chimera2, BPAlignments(rn, i, i + 1), query_gap
                )
            )
            has_bp_assigned[i] = True
    """
	Breakpoint from local alignment i - 1 and i + 1
	"""
    for i in range(1, len(alignments) - 1):
        """
		Add unmatched breakpoint to new_bp_list
		"""
        if has_bp_assigned[i - 1] or has_bp_assigned[i]:
            continue

        prev, curr, next = alignments[i - 1], alignments[i], alignments[i + 1]

        if (
            curr.mapq >= gap_mapq
            or prev.mapq < min_mapq
            or next.mapq < min_mapq
        ):
            continue
        if max_nm is not None and (
            prev.edit_dist >= max_nm or next.edit_dist >= max_nm
        ):
            continue

        prev_intv, next_intv = prev.ref_interval, next.ref_interval
        if (io1 := interval_overlap_l(prev_intv, intrvls)) is None:
            continue
        if (io2 := interval_overlap_l(next_intv, intrvls)) is None:
            continue
        if io1 != io2:
            continue

        # TODO: verify that we can just do abs(grr) instead of flipping order based on strand
        query_gap = next.query_bounds.start - prev.query_bounds.end
        ref_gap = abs(next_intv.start - prev_intv.end)
        if prev_intv.strand == next_intv.strand or (
            abs(ref_gap - query_gap) > max(gap_, abs(query_gap * 0.2))
        ):
            bp_list.append(
                interval2bp(
                    prev, next, BPAlignments(rn, i - 1, i + 1), query_gap
                )
            )
            has_bp_assigned[i] = True

    return bp_list


def filter_small_breakpoint_clusters(
    chr_to_cni_to_reads: defaultdict[str, defaultdict[int, set[str]]],
    min_cluster_cutoff: int,
) -> defaultdict[str, defaultdict[int, set[str]]]:
    """
    Fltering potential breakpoint clusters based on min # of supporting reads.
    """
    chr_to_cni_to_reads_filtered: defaultdict[
        str, defaultdict[int, set[str]]
    ] = defaultdict(lambda: defaultdict(set))

    for chr_, cni_to_reads in chr_to_cni_to_reads.items():
        for cni, reads in cni_to_reads.items():
            if len(reads) >= min_cluster_cutoff:
                chr_to_cni_to_reads_filtered[chr_][cni] = reads
    return chr_to_cni_to_reads_filtered


def cluster_bp_list(
    bp_list: list[Breakpoint], min_cluster_size: int, bp_distance_cutoff: int
) -> list[list[Breakpoint]]:
    """
    Clustering the breakpoints in bp_list
    """
    bp_dict: dict[ChrPairOrientation, list[int]] = defaultdict(list)
    for bpi, bp in enumerate(bp_list):
        bp_dict[
            ChrPairOrientation(
                bp.node1.chr, bp.node2.chr, bp.node1.strand, bp.node2.strand
            )
        ].append(bpi)

    bp_clusters: list[list[Breakpoint]] = []
    for chr_pair, bp_idxs in bp_dict.items():
        if len(bp_idxs) >= min_cluster_size:
            bp_clusters_: list[list[Breakpoint]] = []
            for bp_idx in bp_idxs:
                bp = bp_list[bp_idx]
                bpcim = -1
                for bpci in range(len(bp_clusters_)):
                    for lbp in bp_clusters_[bpci]:
                        if (
                            abs(bp.node1.pos - lbp.node1.pos)
                            < bp_distance_cutoff
                            and abs(bp.node2.pos - lbp.node2.pos)
                            < bp_distance_cutoff
                        ):
                            bpcim = bpci
                            break
                    if bpcim >= 0:
                        break
                if bpcim >= 0:
                    bp_clusters_[bpcim].append(bp)
                else:
                    bp_clusters_.append([bp])
            bp_clusters += bp_clusters_
        else:
            bp_clusters.append([bp_list[bpi] for bpi in bp_dict[chr_pair]])
    return bp_clusters


def interval2bp(
    chimera1: ChimericAlignment,
    chimera2: ChimericAlignment,
    read_info: BPAlignments,
    gap: int = 0,
) -> Breakpoint:
    """
    Convert split/chimeric alignment to breakpoint
    """
    intv1, intv2 = chimera1.ref_interval, chimera2.ref_interval
    chr_idx1, chr_idx2 = (
        CHR_TAG_TO_IDX[intv1.chr],
        CHR_TAG_TO_IDX[intv2.chr],
    )
    if (chr_idx2 < chr_idx1) or (
        chr_idx2 == chr_idx1 and intv2.start < intv1.end
    ):
        return Breakpoint(
            Node(intv1.chr, intv1.end, intv1.strand),
            Node(intv2.chr, intv2.start, intv2.strand.inverse),
            read_info,
            gap,
            was_reversed=True,
            mapq1=chimera1.mapq,
            mapq2=chimera2.mapq,
        )
    return Breakpoint(
        Node(intv2.chr, intv2.start, intv2.strand.inverse),
        Node(intv1.chr, intv1.end, intv1.strand),
        BPAlignments(
            read_info.name, read_info.alignment2, read_info.alignment1
        ),
        gap,
        was_reversed=False,
        mapq1=chimera1.mapq,
        mapq2=chimera2.mapq,
    )


def bpc2bp(bp_cluster: list[Breakpoint], bp_distance_cutoff: float):
    """
    Call exact breakpoint from a breakpoint cluster
    """
    # Copy the first BP to avoid mutating the original BP
    first_bp = bp_cluster[0]

    pos1 = 0 if first_bp.strand1 == Strand.FORWARD else 1_000_000_000
    pos2 = 0 if first_bp.strand2 == Strand.FORWARD else 1_000_000_000
    central_bp: Breakpoint = Breakpoint(
        Node(first_bp.node1.chr, pos1, first_bp.strand1),
        Node(first_bp.node2.chr, pos2, first_bp.strand2),
        # We only use the node information for the central BP
        BPAlignments(
            first_bp.alignment_info.name,
            first_bp.alignment_info.alignment1,
            first_bp.alignment_info.alignment2,
        ),
        first_bp.gap,
        first_bp.was_reversed,
        first_bp.mapq1,
        first_bp.mapq2,
    )

    bpr: list[BPAlignments] = []

    # Calculate basic dist. stats (mean/std) for breakpoints
    bp_stats = BreakpointStats(bp_distance_cutoff)
    for bp_ in bp_cluster:
        bp_stats.observe(bp_)

    bp_starts = []
    bp_ends = []

    for bp_ in bp_cluster:
        if bp_.pos1 > bp_stats.start.mean + 3 * bp_stats.start_window:
            continue
        if bp_.pos1 < bp_stats.start.mean - 3 * bp_stats.start_window:
            continue
        if bp_.pos2 > bp_stats.end.mean + 3 * bp_stats.end_window:
            continue
        if bp_.pos2 < bp_stats.end.mean - 3 * bp_stats.end_window:
            continue
        bp_starts.append(bp_.pos1)
        bp_ends.append(bp_.pos2)

    pos1 = central_bp.pos1
    if len(bp_starts) > 0:
        bp_start_ctr = Counter(bp_starts)
        if (
            len(bp_start_ctr.most_common(2)) == 1
            or bp_start_ctr.most_common(2)[0][1]
            > bp_start_ctr.most_common(2)[1][1]
        ):
            pos1 = bp_start_ctr.most_common(2)[0][0]
        elif len(bp_starts) % 2 == 1:
            pos1 = int(np.median(bp_starts))
        elif bp_cluster[
            -1
        ].strand1.is_forward:  # TODO: using last BP seems strange and below, but maintaining prev behavior
            pos1 = int(np.ceil(np.median(bp_starts)))
        else:
            pos1 = int(np.floor(np.median(bp_starts)))

    pos2 = central_bp.pos2
    if len(bp_ends) > 0:
        bp_end_ctr = Counter(bp_ends)
        if (
            len(bp_end_ctr.most_common(2)) == 1
            or bp_end_ctr.most_common(2)[0][1] > bp_end_ctr.most_common(2)[1][1]
        ):
            pos2 = bp_end_ctr.most_common(2)[0][0]
        elif len(bp_ends) % 2 == 1:
            pos2 = int(np.median(bp_ends))
        elif bp_cluster[-1].strand2.is_forward:
            pos2 = int(np.ceil(np.median(bp_ends)))
        else:
            pos2 = int(np.floor(np.median(bp_ends)))
    central_bp.node1 = Node(central_bp.node1.chr, pos1, central_bp.node1.strand)
    central_bp.node2 = Node(central_bp.node2.chr, pos2, central_bp.node2.strand)

    final_bp_cluster: list[Breakpoint] = []

    final_bp_stats = BreakpointStats(bp_distance_cutoff)
    for bp in bp_cluster:
        if bp_match(bp, central_bp, bp.gap * 1.2, bp_distance_cutoff):
            bpr.append(bp.alignment_info)
            final_bp_stats.observe(bp)
        else:
            final_bp_cluster.append(bp)
    if len(bpr) == 0:
        return central_bp, bpr, final_bp_stats, []
    return central_bp, bpr, final_bp_stats, final_bp_cluster


def bp_match(
    bp1: Breakpoint, bp2: Breakpoint, rgap: float, bp_dist_cutoff: float
):
    """
    Check if two breakpoints match
    A breakpoint (chr1, e1, chr2, s2) must either satisfy
        1. chr1 > chr2
        2. chr1 == chr2 and e1 >= s2
    """
    if bp1.chr1 != bp2.chr1 or bp1.chr2 != bp2.chr2:
        return False
    if bp1.strand1 != bp2.strand1 or bp1.strand2 != bp2.strand2:
        return False
    if rgap <= 0:
        return (
            abs(bp1.pos1 - bp2.pos1) < bp_dist_cutoff
            and abs(bp1.pos2 - bp2.pos2) < bp_dist_cutoff
        )

    rgap_ = rgap
    consume_rgap = [0, 0]
    if bp1.strand1.is_forward and bp1.pos1 <= bp2.pos1 - bp_dist_cutoff:
        rgap_ -= bp2.pos1 - bp_dist_cutoff - bp1.pos1 + 1
        consume_rgap[0] = 1
    if bp1.strand1.is_reverse and bp1.pos1 >= bp2.pos1 + bp_dist_cutoff:
        rgap_ -= bp1.pos1 - bp2.pos1 - bp_dist_cutoff + 1
        consume_rgap[0] = 1
    if bp1.strand2.is_forward and bp1.pos2 <= bp2.pos2 - bp_dist_cutoff:
        rgap_ -= bp2.pos2 - bp_dist_cutoff - bp1.pos2 + 1
        consume_rgap[1] = 1
    if bp1.strand2.is_reverse and bp1.pos2 >= bp2.pos2 + bp_dist_cutoff:
        rgap_ -= bp1.pos2 - bp2.pos2 - bp_dist_cutoff + 1
        consume_rgap[1] = 1
    return (
        (consume_rgap[0] == 1 and rgap_ >= 0)
        or (abs(bp1.pos1 - bp2.pos1)) < bp_dist_cutoff
    ) and (
        (consume_rgap[1] == 1 and rgap_ >= 0)
        or (abs(bp1.pos2 - bp2.pos2)) < bp_dist_cutoff
    )


def sort_chrom_names(chromlist: List[str]) -> List[str]:
    # TODO: use CHR_TAG_TO_IDX instead of this method?
    def sort_key(x: str):
        chr_val = x[3] if x.startswith("chr") else x
        return int(chr_val) if chr_val.isnumeric() else ord(chr_val)

    return sorted(chromlist, key=sort_key)


def get_intervals_from_seed_file(
    seed_file: io.TextIOWrapper,
) -> List[AmpliconInterval]:
    intervals = []
    for line in seed_file:
        seed = line.strip().split()
        intervals.append(AmpliconInterval(seed[0], int(seed[1]), int(seed[2])))
    return intervals


def get_interval_from_bed(file_row: Tuple[str, str, str]) -> List:
    return [file_row[0], int(file_row[1]), int(file_row[2]), -1]


def get_interval_from_cns(file_row: Tuple[str, str, str]) -> List:
    return [file_row[0], int(file_row[1]), int(file_row[2]), -1]


def check_valid_discordant_rc_partition(
    rc_list: list[int], partition: list[int], max_multiplicity: int = 5
) -> tuple[int, float] | None:
    """Verify if a given partition of discordant edges meets the CoRAL criteria.

    We impose a maximum number of traversals through a discordant edge by each
    walk, otherwise cycles cannot be detected. This is empirically determined
    per discordant edge by clustering the various reads (+ their associated
    counts) that support that edge. The `max_multiplicity` is used to separate
    supporting reads into different clusters when one read has a significantly
    higher count than another (likely indicating multiple ecDNA species).

    For additional detail, refer to the following from our paper:
    Supplemental Equation S4.8 - https://genome.cshlp.org/content/suppl/2024/09/27/gr.279131.124.DC1/Supplemental_Methods.pdf
    Supplemental Figure 13 - https://genome.cshlp.org/content/suppl/2024/09/27/gr.279131.124.DC1/Supplemental_Figures_.pdf

    Args:
        rc_list: Read counts for each discordant edge
        partition: Partition of the discordant edge indices
        max_multiplicity: Maximum allowed multiplicity of the discordant edge

    Returns:
        If the partition is invalid, returns None. Otherwise, returns a tuple
        containing:
            - The empirical value of R for the discordant edge
            - The score of the partition
    """
    if partition[0] == partition[1]:
        return (partition[0], 0.0)
    rc_list_p = rc_list[partition[0] : partition[1] + 1]
    if rc_list_p[-1] < rc_list_p[0] * 2.0:
        return (partition[1], 0.0)
    base_ri = 0
    while base_ri < len(rc_list_p) and rc_list_p[base_ri] < rc_list_p[0] * 2.0:
        base_ri += 1
    base_avg_rc = cast(float, np.average(rc_list_p[:base_ri]))
    if rc_list_p[-1] / base_avg_rc >= max_multiplicity + 0.5:
        return None
    score = -10.0
    best_ri = base_ri
    sum_deviation = 1.0
    for base_ri_ in range(base_ri, 0, -1):
        base_avg_rc = cast(float, np.average(rc_list_p[:base_ri_]))
        base_size = len(rc_list_p[:base_ri_])
        sizes = {}
        # Cluster breakpoints with higher multiplicities
        li = base_ri_
        multiplicity = 2
        if rc_list_p[base_ri_] / base_avg_rc < multiplicity - 0.5:
            continue
        while rc_list_p[base_ri_] / base_avg_rc >= multiplicity + 0.5:
            multiplicity += 1
        sum_gap = np.log2(rc_list_p[base_ri_]) - np.log2(
            rc_list_p[base_ri_ - 1]
        )
        # Note: sometimes int(round()) will convert 1.5 to 1
        # The following procedure works well
        for i in range(base_ri_, len(rc_list_p)):
            if rc_list_p[i] / base_avg_rc >= multiplicity + 0.5:
                sum_gap += np.log2(rc_list_p[i]) - np.log2(rc_list_p[i - 1])
                sizes[multiplicity] = [li, i - 1]
                li = i
                while rc_list_p[i] / base_avg_rc >= multiplicity + 0.5:
                    multiplicity += 1
        sizes[multiplicity] = [li, len(rc_list_p) - 1]
        if multiplicity > max_multiplicity:
            continue
        size_flag = True
        for m in range(2, multiplicity + 1):
            if m in sizes and sizes[m][1] - sizes[m][0] >= base_size:
                size_flag = False
                break
        if not size_flag:
            continue
        sum_deviation_ = sum(
            [
                np.abs(
                    m
                    - np.average(
                        rc_list_p[sizes[m][0] : sizes[m][1] + 1] / base_avg_rc  # type: ignore[operator]
                    )
                )
                for m in range(2, multiplicity + 1)
                if m in sizes
            ],
            0,
        )
        if sum_gap - sum_deviation_ > score:
            score = sum_gap - sum_deviation_
            sum_deviation = sum_deviation_
            best_ri = base_ri_
    if sum_deviation < 1.0:
        return (best_ri + partition[0] - 1, score)
    return None


def enumerate_partitions(
    k: int, start: int, end: int
) -> Generator[list[list[int]]]:
    """Generate all partitions of the interval [start, end] into k parts."""
    if k == 0:
        yield [[start, end]]
    else:
        for i in range(1, end - start - k + 2):
            for res in enumerate_partitions(k - 1, start + i, end):
                yield [[start, start + i - 1], *res]


def output_breakpoint_graph_lr(g: BreakpointGraph, ogfile: str) -> None:
    """Write a breakpoint graph to file in AA graph format with only long read information"""
    with open(ogfile, "w") as fp:
        fp.write(
            "SequenceEdge: StartPosition, EndPosition, PredictedCN, "
            "AverageCoverage, Size, NumberOfLongReads\n",
        )
        for se in g.sequence_edges:
            fp.write(
                f"sequence\t{se.chr}:{se.start}-\t{se.chr}:{se.end}+\t"
                f"{se.cn}\t{se.lr_nc * 1.0 / se.gap}\t{se.gap}\t{se.lr_count}\n"
            )
        fp.write(
            "BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads\n"
        )
        for srce in g.source_edges:
            fp.write(
                f"source\t0-->{srce.node.chr}:{srce.node.pos}{srce.node.strand}\t"
                f"{srce.cn}\t-1\n"
            )
        for ce in g.concordant_edges:
            fp.write(
                f"concordant\t{ce.node1.chr}:{ce.node1.pos}{ce.node1.strand}->"
                f"{ce.node2.chr}:{ce.node2.pos}{ce.node2.strand}\t"
                f"{ce.cn}\t{ce.lr_count}\n"
            )
        for de in g.discordant_edges:
            fp.write(
                f"discordant\t{de.node1.chr}:{de.node1.pos}{de.node1.strand}->"
                f"{de.node2.chr}:{de.node2.pos}{de.node2.strand}\t"
                f"{de.cn}\t{de.lr_count}\n"
            )
        # Only output the longest path constraints (post subpath filtering)
        fp.write("PathConstraint: Path, Support\n")
        for pc in g.longest_path_constraints:
            basic_pc = g.path_constraints[pc.pc_idx]
            fp.write(
                f"path_constraint\t"
                f"{core_utils.path_to_str(basic_pc.path, pc.edge_counts)}\t"
                f"{pc.support}\n"
            )
        fp.write("AmpliconIntervals: chr, start, end\n")
        for ai in g.amplicon_intervals:
            fp.write(f"interval\t{ai.chr}\t{ai.start}\t{ai.end}\n")


def output_breakpoint_info_lr(
    g: BreakpointGraph, filename: str, bp_stats: BreakpointStats
) -> None:
    """Write the list of breakpoints to file"""
    with open(filename, "w") as fp:
        fp.write(
            "chr1\tpos1\tchr2\tpos2\torientation\tlr_support\t"
            "lr_info=[avg1, avg2, std1, std2, mapq1, mapq2]\n",
        )
        for di in range(len(g.discordant_edges)):
            de = g.discordant_edges[di]
            fp.write(
                f"{de.node2.chr}\t{de.node2.pos}\t{de.node1.chr}\t{de.node1.pos}\t"
                f"{de.node2.strand}{de.node1.strand}\t{de.lr_count}\t{bp_stats[di]}\n"
            )


# TODO: further modularize the utilities in this file
def fetch_breakpoint_reads(
    bam_file: pysam.AlignmentFile,
) -> tuple[dict[str, list[ChimericAlignment]], datatypes.BasicStatTracker]:
    read_name_to_length: dict[str, int] = {}
    chimeric_strings: dict[str, list[str]] = defaultdict(list)
    edit_dist_stats = datatypes.BasicStatTracker()

    read: bam_types.BAMRead
    for read in bam_file.fetch():  # type: ignore[assignment]
        rn: str = read.query_name  # type: ignore[assignment]
        if read.flag < 256 and rn not in read_name_to_length:
            read_name_to_length[rn] = read.query_length
        try:
            sa_list = read.get_tag("SA:Z:")[:-1].split(";")  # type: ignore
            for sa in sa_list:
                if sa not in chimeric_strings[rn]:
                    chimeric_strings[rn].append(sa)
        except:
            # Highest mapq value from minimap2, = "most reliable" alignments
            if read.mapping_quality == 60:
                e = read.get_cigar_stats()[0][-1] / read.query_length
                edit_dist_stats.observe(e)
    logger.info(f"Fetched {len(chimeric_strings)} chimeric reads.")
    reads_wo_primary_alignment = []
    chimeric_alignments = {}
    for r in chimeric_strings:
        if r not in read_name_to_length:
            logger.warning(
                f"Found chimeric read name {r} without primary alignment; Read length: N/A."
            )
            logger.warning(f"All CIGAR strings: {chimeric_strings[r]}.")
            reads_wo_primary_alignment.append(r)
            continue
        chimeric_alignments[r] = cigar_parsing.alignment_from_satags(
            chimeric_strings[r], r
        )
    for r in reads_wo_primary_alignment:
        del chimeric_alignments[r]
    logger.info(
        f"Computed alignment intervals on all {len(chimeric_strings)} chimeric reads.",
    )
    return chimeric_alignments, edit_dist_stats


def get_cns_idx_intv_to_reads(
    d1_segs: dict[str, dict[int, set[str]]],
    cns_intervals_by_chr: dict[str, list[CNSInterval]],
):
    pass


def get_indel_alignments_from_read(
    chr_tag: core_types.ChrTag,
    read: bam_types.SAMSegment,
    min_del_len: int,
    nm_threshold: Optional[float] = None,
) -> list[datatypes.LargeIndelAlignment]:
    if read.mapping_quality < 20:
        return []
    indel_alignments = []
    blocks = read.get_blocks()

    passes_edit_filter = True
    if nm_threshold is not None:
        agg_del = 0
        for block, next_block in zip(blocks, blocks[1:]):
            gap = abs(next_block[0] - block[1])
            if gap > min_del_len:
                agg_del += gap
        read_nm = read.get_cigar_stats()[0][-1]  # Edit dist
        read_nm = (read_nm - agg_del) / read.query_length  # type: ignore[assignment]
        passes_edit_filter = read_nm < nm_threshold

    if passes_edit_filter:
        for block, next_block in zip(blocks, blocks[1:]):
            gap = abs(next_block[0] - block[1])
            if gap > min_del_len:
                indel_alignments.append(
                    datatypes.LargeIndelAlignment(
                        chr_tag,
                        next_block[0],
                        block[1],
                        blocks[0][0],
                        blocks[-1][1],
                        read.mapping_quality,
                    )
                )
    return indel_alignments


def does_bp_edge_join_sequence_edges(
    bp_edge: datatypes.ConcordantEdge | datatypes.DiscordantEdge,
    seq_edge_1: datatypes.SequenceEdge,
    seq_edge_2: datatypes.SequenceEdge,
) -> bool:
    return (
        bp_edge.node1.chr == seq_edge_1.chr
        and bp_edge.node1.pos == seq_edge_1.end
        and bp_edge.node2.pos == seq_edge_2.start
    )
