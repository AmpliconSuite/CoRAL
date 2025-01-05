"""Infer breakpoint graph(s) from long read alignments + associated CN calls."""

from __future__ import annotations

import functools
import io
import logging
import pathlib
import pickle
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, List

import intervaltree  # type: ignore[import-untyped]
import numpy as np
import pysam
import typer

from coral import cigar_parsing, types
from coral.breakpoint import breakpoint_utilities
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.breakpoint_types import CNSSegData
from coral.breakpoint.breakpoint_utilities import (
    alignment2bp,
    alignment2bp_l,
    bpc2bp,
    cluster_bp_list,
    interval_adjacent,
    interval_exclusive,
    interval_overlap,
    interval_overlap_l,
)
from coral.constants import CHR_TAG_TO_IDX
from coral.datatypes import (
    AmpliconInterval,
    Breakpoint,
    ChimericAlignment,
    CNSInterval,
    CNSIntervalTree,
    Interval,
    ReadInterval,
    WalkData,
)
from coral.models import path_constraints
from coral.types import Chr, CNSIdx, ReadName

edge_type_to_index = {"s": 0, "c": 1, "d": 2}


logger = logging.getLogger(__name__)


@dataclass
class LongReadBamToBreakpointMetadata:
    lr_bamfh: pysam.AlignmentFile  # Long read bam file

    lr_graph: list[BreakpointGraph] = field(
        default_factory=list
    )  # Breakpoint graph

    # Tunable hyperparameters
    max_seq_len: int = 2000000  # Maximum allowed length for a breakpoint edge
    cn_gain: float = 5.0
    min_bp_cov_factor: float = 1.0
    min_bp_match_cutoff_: int = 100
    interval_delta: int = 100000

    min_cluster_cutoff: int = (
        3  # Hard cutoff for considering a long read breakpoint cluster
    )
    max_breakpoint_distance_cutoff: int = 2000  # Used for breakpoint clustering - if the distance of two breakpoint positions are greater than this cutoff, then start a new cluster
    min_del_len: int = 600  # The minimum length of all +- (deletion) breakpoints returned by AA
    # small_del_cutoff = 10000 # +- breakpoints (small deletions) with the two ends less than this cutoff are treated specially

    # Map read name -> chimeric alignments (two or more records for one read)
    chimeric_alignments: dict[str, list[ChimericAlignment]] = field(
        default_factory=dict
    )
    # Map chr tag -> CN segment idx -> reads
    chr_cns_to_chimeras: dict[str, dict[int, list[ReadInterval]]] = field(
        default_factory=dict
    )

    # Map read name -> alignments with one record per read but large indels showing in CIGAR string
    large_indel_alignments: dict[str, list[pysam.AlignedSegment]] = field(
        default_factory=dict
    )
    # For edit distance filter of breakpoints
    nm_stats: list[float] = field(
        default_factory=functools.partial(lambda: [0.0] * 3)
    )  # Default of [0.0, 0.0, 0.0]
    nm_filter: bool = False

    # AA amplicon intervals
    amplicon_intervals: list[AmpliconInterval] = field(default_factory=list)
    amplicon_interval_connections: dict[
        tuple[int, int], set[AmpliconInterval]
    ] = field(default_factory=dict)

    cns_intervals: list[CNSInterval] = field(
        default_factory=list
    )  # CN segments
    cns_intervals_by_chr: dict[str, list[CNSInterval]] = field(
        default_factory=dict
    )
    # log2 CN for each CN segment
    log2_cn: list[float] = field(default_factory=list)
    # Interval tree structure for each chromosome
    cns_tree: dict[str, CNSIntervalTree] = field(default_factory=dict)
    # Normal long read coverage - used for CN assignment
    normal_cov: float = 0.0

    ccid2id: dict[int, int] = field(default_factory=dict)
    # List of breakpoints (discordant edges)
    new_bp_list: list[Breakpoint] = field(default_factory=list)
    # Statistics of breakpoints (discordant edges)
    new_bp_stats: list[list[int]] = field(default_factory=list)
    new_bp_ccids: list[int] = field(default_factory=list)
    source_edges: list[list[Any]] = field(default_factory=list)
    source_edge_ccids: list[int] = field(default_factory=list)
    path_constraints: dict[int, list[list[Any]]] = field(default_factory=dict)
    longest_path_constraints: dict[int, list[list[Any]]] = field(
        default_factory=dict
    )
    walks_by_amplicon: dict[int, WalkData[types.AmpliconWalk]] = field(
        default_factory=dict
    )
    walk_weights_by_amplicon: dict[int, WalkData[float]] = field(
        default_factory=dict
    )
    path_constraints_satisfied: dict[int, WalkData[list[int]]] = field(
        default_factory=dict
    )

    def set_raw_cns_data(self, cns_data: CNSSegData) -> None:
        self.cns_tree = cns_data.tree
        self.cns_intervals = cns_data.intervals
        self.cns_intervals_by_chr = cns_data.intervals_by_chr
        self.log2_cn = cns_data.log2_cn

    def read_cns(self, cns_file: io.TextIOWrapper) -> None:
        """Read in (cnvkit) *.cns file and estimate the normal long read coverage"""
        self.set_raw_cns_data(CNSSegData.from_file(cns_file))

        logger.debug(f"Total num LR copy number segments:{len(self.log2_cn)}.")
        log2_cn_order = np.argsort(self.log2_cn)
        cns_intervals_median = []
        log2_cn_median = []
        im = int(len(log2_cn_order) / 2.4)
        ip = im + 1
        total_int_len = 0
        cns_intervals_median.append(self.cns_intervals[log2_cn_order[ip]])
        cns_intervals_median.append(self.cns_intervals[log2_cn_order[im]])
        log2_cn_median.append(self.log2_cn[log2_cn_order[ip]])
        log2_cn_median.append(self.log2_cn[log2_cn_order[im]])
        total_int_len += (
            self.cns_intervals[log2_cn_order[ip]][2]
            - self.cns_intervals[log2_cn_order[ip]][1]
            + 1
        )
        total_int_len += (
            self.cns_intervals[log2_cn_order[im]][2]
            - self.cns_intervals[log2_cn_order[im]][1]
            + 1
        )
        i = 1
        while total_int_len < 10000000:
            cns_intervals_median.append(
                self.cns_intervals[log2_cn_order[ip + i]]
            )
            cns_intervals_median.append(
                self.cns_intervals[log2_cn_order[im - i]]
            )
            log2_cn_median.append(self.log2_cn[log2_cn_order[ip]])
            log2_cn_median.append(self.log2_cn[log2_cn_order[im]])
            total_int_len += (
                self.cns_intervals[log2_cn_order[ip + i]][2]
                - self.cns_intervals[log2_cn_order[ip + i]][1]
                + 1
            )
            total_int_len += (
                self.cns_intervals[log2_cn_order[im - i]][2]
                - self.cns_intervals[log2_cn_order[im - i]][1]
                + 1
            )
            i += 1
        logger.debug(
            f"Use {len(cns_intervals_median)} LR copy number segments."
        )
        logger.debug(
            f"Total length of LR copy number segments: {total_int_len}."
        )
        logger.debug(f"Average LR copy number: {np.average(log2_cn_median)}.")
        nnc = 0
        for i in range(len(cns_intervals_median)):
            nnc += sum(
                [
                    sum(nc)
                    for nc in self.lr_bamfh.count_coverage(
                        cns_intervals_median[i][0],
                        cns_intervals_median[i][1],
                        cns_intervals_median[i][2] + 1,
                        quality_threshold=0,
                        read_callback="nofilter",
                    )
                ]
            )
        self.normal_cov = nnc * 1.0 / total_int_len
        logger.info(
            f"LR normal cov ={self.normal_cov}, {nnc=}, {total_int_len=}."
        )
        self.min_cluster_cutoff = max(
            self.min_cluster_cutoff,
            int(self.min_bp_cov_factor * self.normal_cov),
        )
        logger.debug(f"Reset min_cluster_cutoff to {self.min_cluster_cutoff}.")

    def pos2cni(self, chr_tag: str, pos: int) -> set[int]:
        return self.cns_tree[chr_tag][pos]

    def hash_alignment_to_seg(
        self, chimeras_by_read: dict[str, list[ChimericAlignment]]
    ):
        """Speed up amplified interval search by hashing chimeric alignments
        from each long read to CN segments."""
        chr_cns_to_chimeras: dict[str, dict[int, list[ReadInterval]]] = (
            defaultdict(lambda: defaultdict(list))
        )
        for alignments in chimeras_by_read.values():
            for alignment in alignments:
                read_interval = alignment.ref_interval
                if (chr_tag := read_interval.chr_tag) in self.cns_tree:
                    cn_seg_idxs = self.cns_tree[chr_tag].get_cn_segment_indices(
                        read_interval
                    )
                    alignment.cns = cn_seg_idxs
                    for cni in cn_seg_idxs:
                        chr_cns_to_chimeras[chr_tag][cni].append(read_interval)
        logger.info("Completed hashing chimeric reads to CN segments.")
        self.chr_cns_to_chimeras = chr_cns_to_chimeras
        self.chimeric_alignments = chimeras_by_read
        return chimeras_by_read, chr_cns_to_chimeras

    def widen_seed_intervals(self) -> None:
        """Widen seed intervals to fully encompass CN segments that the interval falls within."""
        for interval in self.amplicon_intervals:
            chr_tag = interval.chr_tag
            lcni, rcni = self.cns_tree[chr_tag].get_cns_ends(interval)
            interval.start = self.cns_intervals_by_chr[chr_tag][lcni].start
            if self.pos2cni(
                chr_tag,
                interval.start - self.interval_delta,
            ):
                interval.start = (
                    self.cns_intervals_by_chr[chr_tag][lcni].start
                    - self.interval_delta
                )
            interval.end = self.cns_intervals_by_chr[chr_tag][rcni].end
            if self.pos2cni(
                chr_tag,
                interval.end + self.interval_delta,
            ):
                interval.end = (
                    self.cns_intervals_by_chr[chr_tag][rcni].end
                    + self.interval_delta
                )
            logger.debug(f"\tUpdated amplicon interval {interval}")

    def merge_amplicon_intervals(self) -> None:
        lastai = 0
        intervals_to_merge: list[tuple[int, int]] = []
        # Iterate through amplicon intervals in pairs, step size 1
        for ai, (curr, next) in enumerate(
            zip(self.amplicon_intervals, self.amplicon_intervals[1:])
        ):
            if not (curr.is_adjacent(next) or curr.does_overlap(next)):
                if ai > lastai:
                    logger.debug(
                        "Merging intervals from %d to %d." % (lastai, ai),
                    )
                    intervals_to_merge.append((lastai, ai))
                lastai = ai + 1
        num_amplicons = len(self.amplicon_intervals)
        if num_amplicons and lastai < num_amplicons - 1:
            logger.debug(
                f"Merging intervals from {lastai} to {num_amplicons-1}."
            )
            intervals_to_merge.append((lastai, num_amplicons - 1))
        for ai1, ai2 in intervals_to_merge[::-1]:
            start_intv, end_intv = (
                self.amplicon_intervals[ai1],
                self.amplicon_intervals[ai2],
            )
            # Reset interval
            start_intv.end = end_intv.end
            logger.debug(f"Reset amplicon interval {ai} to {start_intv}.")
            # Modify ccid
            for ai in range(ai1 + 1, ai2 + 1):
                intermediate_intv = self.amplicon_intervals[ai]
                if intermediate_intv.amplicon_id != start_intv.amplicon_id:
                    ccid = intermediate_intv.amplicon_id
                    logger.debug(
                        "Reset amplicon intervals with ccid %d." % ccid,
                    )
                    for ai_, interval in enumerate(self.amplicon_intervals):
                        if not interval.amplicon_id == ccid:
                            continue
                        interval.amplicon_id = start_intv.amplicon_id
                        logger.debug(
                            "Reset ccid of amplicon interval %d to %d."
                            % (
                                ai_,
                                start_intv.amplicon_id,
                            ),
                        )
                        logger.debug(f"Updated amplicon interval: {interval}")
            # Modify interval connections
            connection_map = {}
            for connection in self.amplicon_interval_connections:
                connection_map[connection] = connection
            for ai in range(ai1 + 1, ai2 + 1):
                ai_unsorted_ = sorted_ai_indices[ai1]
                ai_unsorted = sorted_ai_indices[ai]
                for connection in connection_map:
                    if ai_unsorted == connection_map[connection][0]:
                        connection_map[connection] = (
                            ai_unsorted_,
                            connection_map[connection][1],
                        )
                    if ai_unsorted == connection_map[connection][1]:
                        connection_map[connection] = (
                            connection_map[connection][0],
                            ai_unsorted_,
                        )
                    if (
                        connection_map[connection][1]
                        < connection_map[connection][0]
                    ):
                        connection_map[connection] = (
                            connection_map[connection][1],
                            connection_map[connection][0],
                        )
            for connection in connection_map:
                if connection != connection_map[connection]:
                    logger.debug(
                        f"Reset connection between amplicon intervals {connection} to {connection_map[connection]}."
                    )
                    if (
                        connection_map[connection]
                        not in self.amplicon_interval_connections
                    ):
                        self.amplicon_interval_connections[
                            connection_map[connection]
                        ] = self.amplicon_interval_connections[connection]
                    else:
                        self.amplicon_interval_connections[
                            connection_map[connection]
                        ] |= self.amplicon_interval_connections[connection]
                    del self.amplicon_interval_connections[connection]
                    if (
                        connection_map[connection][0]
                        == connection_map[connection][1]
                    ):
                        del self.amplicon_interval_connections[
                            connection_map[connection]
                        ]
            # Delete intervals

            for ai in range(ai1 + 1, ai2 + 1)[::-1]:
                ai_unsorted = sorted_ai_indices[ai]
                logger.debug(
                    f"Delete amplicon interval {ai_unsorted} - {self.amplicon_intervals[ai_unsorted]}."
                )
                del amplicon_intervals_sorted[ai]
                del sorted_ai_indices[ai]

    def find_amplicon_intervals(self):
        # Reset seed intervals
        logger.debug("Updating seed amplicon intervals based on CN segments.")
        self.widen_seed_intervals()
        ccid = 0
        for ai, interval in enumerate(self.amplicon_intervals):
            if interval.amplicon_id == -1:
                logger.debug(f"Begin processing amplicon interval {ai}")
                logger.debug(f"\tAmplicon interval {interval}")
                self.find_interval_i(ai, ccid)
                ccid += 1
        logger.debug(
            f"Identified {len(self.amplicon_intervals)} amplicon intervals in total.",
        )
        self.amplicon_intervals.sort()

        # amplicon_intervals_sorted = [
        #     self.amplicon_intervals[i] for i in sorted_ai_indices
        # ]
        for interval in self.amplicon_intervals:
            logger.debug(f"\tAmplicon {interval=}")

        # Merge amplicon intervals
        logger.debug("Begin merging adjacent intervals.")
        self.merge_amplicon_intervals()

        self.amplicon_intervals = [
            amplicon_intervals_sorted[ai]
            for ai in range(len(amplicon_intervals_sorted))
        ]
        ind_map = {
            sorted_ai_indices[i]: i for i in range(len(sorted_ai_indices))
        }
        connection_map = {
            connection: (
                min(ind_map[connection[0]], ind_map[connection[1]]),
                max(ind_map[connection[0]], ind_map[connection[1]]),
            )
            for connection in self.amplicon_interval_connections
        }
        self.amplicon_interval_connections = {
            connection_map[connection]: self.amplicon_interval_connections[
                connection
            ]
            for connection in self.amplicon_interval_connections
        }
        # Reset ccids
        ai_explored = np.zeros(len(self.amplicon_intervals))
        for ai in range(len(self.amplicon_intervals)):
            ai_ccid = interval.amplicon_id
            if ai_explored[ai] == 0:
                L = [ai]  # BFS queue
                while len(L) > 0:
                    ai_ = L.pop(0)
                    ai_explored[ai_] = 1
                    if interval.amplicon_id != ai_ccid:
                        interval.amplicon_id = ai_ccid
                    for ai1, ai2 in self.amplicon_interval_connections.keys():
                        if ai1 == ai_ and ai_explored[ai2] == 0:
                            L.append(ai2)
                        elif ai2 == ai_ and ai_explored[ai1] == 0:
                            L.append(ai1)

        logger.debug(
            "There are %d amplicon intervals after merging."
            % len(self.amplicon_intervals),
        )
        for ai in range(len(self.amplicon_intervals)):
            logger.debug("\tAmplicon interval %s after merging." % interval)

    def addbp(self, bp_: Breakpoint, bpr_, bp_stats_, ccid):
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            if (
                bp[0] == bp_[0]
                and bp[3] == bp_[3]
                and bp[2] == bp_[2]
                and bp[5] == bp_[5]
                and abs(bp[1] - bp_[1]) < 200
                and abs(bp[4] - bp_[4]) < 200
            ):
                self.new_bp_list[bpi][-1] |= set(bpr_)
                return bpi
        bpi = len(self.new_bp_list)
        self.new_bp_list.append(bp_ + [bpr_])
        self.new_bp_ccids.append(ccid)
        self.new_bp_stats.append(bp_stats_)
        return bpi

    def refine_breakpoints(
        self,
        interval: AmpliconInterval,
        cni_intv: Interval,
        new_bpi_refined: list[int],
    ):
        nint_segs = []
        nint_segs_ = []
        cni_chr = cni_intv.chr_tag
        ns = self.cns_intervals_by_chr[cni_chr][cni_intv.start].start
        ne = self.cns_intervals_by_chr[cni_chr][cni_intv.end].end
        cns_intv = Interval(cni_chr, ns, ne)

        for bpi in new_bpi_refined:
            bp = self.new_bp_list[bpi]
            try:
                if interval.contains(bp.chr1, bp.start) and cns_intv.contains(
                    bp.chr2, bp.end
                ):
                    cns_tree = self.cns_tree[bp.chr2]
                    nint_segs.append(
                        [cns_tree.get_single_cns_idx(bp.end), bp.end, bpi]
                    )
                elif interval.contains(bp.chr1, bp.start) and cns_intv.contains(
                    bp.chr1, bp.start
                ):
                    cns_tree = self.cns_tree[bp.chr1]
                    nint_segs.append(
                        [cns_tree.get_single_cns_idx(bp.start), bp.start, bpi]
                    )
                else:
                    logger.warning(
                        "\t\tExact breakpoint outside amplicon interval."
                    )
                    logger.warning("\t\tBreakpoint %s." % bp)
                    logger.warning("\t\tCurrent interval %s." % interval)
                    logger.warning(f"\t\tNew interval {cns_intv}.")
                    o1 = interval_overlap(
                        [bp[0], bp[1], bp[1]],
                        [cni_chr, ns, ne],
                    )
                    o2 = interval_overlap(
                        [bp[3], bp[4], bp[4]],
                        [cni_chr, ns, ne],
                    )
                    contains_start = cni_intv.contains(bp.chr1, bp.start)
                    contains_end = cni_intv.contains(bp.chr2, bp.end)
                    if o1 and o2:
                        nint_segs.append(
                            [
                                list(self.pos2cni(bp[0], bp[1]))[0].data,
                                bp[1],
                                bpi,
                            ],
                        )
                        nint_segs.append(
                            [
                                list(self.pos2cni(bp[3], bp[4]))[0].data,
                                bp[4],
                                bpi,
                            ],
                        )
                    elif o1:
                        nint_segs.append(
                            [
                                list(self.pos2cni(bp[0], bp[1]))[0].data,
                                bp[1],
                                bpi,
                            ],
                        )
                        nint_segs_.append(
                            [
                                bp[3],
                                list(self.pos2cni(bp[3], bp[4]))[0].data,
                                bp[4],
                                bpi,
                            ],
                        )
                    elif o2:
                        nint_segs_.append(
                            [
                                bp[0],
                                list(self.pos2cni(bp[0], bp[1]))[0].data,
                                bp[1],
                                bpi,
                            ],
                        )
                        nint_segs.append(
                            [
                                list(self.pos2cni(bp[3], bp[4]))[0].data,
                                bp[4],
                                bpi,
                            ],
                        )
                    else:
                        nint_segs_.append(
                            [
                                bp[0],
                                list(self.pos2cni(bp[0], bp[1]))[0].data,
                                bp[1],
                                bpi,
                            ],
                        )
                        nint_segs_.append(
                            [
                                bp[3],
                                list(self.pos2cni(bp[3], bp[4]))[0].data,
                                bp[4],
                                bpi,
                            ],
                        )
            except:
                pass

    def find_interval_i(self, ai: int, ccid: int):
        """Given an amplification interval I indexed by ai, search for
            amplification intervals connected with I iteratively (with BFS)
            by a breakpoint edge
        Assign I a connected component id ccid if not already assigned
        """
        logger.debug(f"\tStart BFS on amplicon interval {ai}.")
        interval_queue: list[int] = [ai]  # BFS queue
        while len(interval_queue) > 0:
            logger.debug(f"\t\tBFS queue: {interval_queue}")
            ai_: int = interval_queue.pop(0)
            interval = self.amplicon_intervals[ai_]
            logger.debug(f"\t\tNext amplicon interval {ai_}: {interval}.")
            chr = interval.chr_tag
            if interval.amplicon_id == -1:
                interval.amplicon_id = ccid
            logger.debug(f"\t\tReset connected component ID to {ccid}")

            # Identify all amplification intervals connected to interval indexed
            # by ai_ with a breakpoint edge
            try:
                si, ei = self.cns_tree[interval.chr_tag].get_cns_ends(interval)
            except KeyError:
                continue

            # Chr -> CNS idxs sharing a chimerical alignment -> associated reads
            chr_to_cns_to_reads: dict[Chr, dict[CNSIdx, set[ReadName]]] = (
                defaultdict(lambda: defaultdict(set))
            )
            for cns_idx in range(si, ei + 1):
                if cns_idx in self.chr_cns_to_chimeras[chr]:
                    for read_intv in self.chr_cns_to_chimeras[chr][cns_idx]:
                        read_name = read_intv.name
                        for alignment in self.chimeric_alignments[read_name]:
                            read_intv = alignment.ref_interval
                            if read_intv.chr_tag != chr:
                                continue

                            for cns_idx1 in alignment.cns:
                                if cns_idx1 <= si or cns_idx1 >= ei:
                                    chr_to_cns_to_reads[chr][cns_idx1].add(
                                        read_name
                                    )

            chr_to_cns_to_reads = (
                breakpoint_utilities.filter_small_breakpoint_clusters(
                    chr_to_cns_to_reads, self.min_cluster_cutoff
                )
            )

            new_intervals_refined = []
            new_intervals_connections = []
            for chr_ in chr_to_cns_to_reads:
                print(f"processing new d1 seg {chr_}")
                logger.debug(f"\t\tFound new intervals on chr {chr_=}")
                # Initial list of new amplicon intervals, based on CN seg idxs
                cni_intv_to_reads: dict[Interval, set[str]] = defaultdict(set)
                sorted_cns_idxs = sorted(chr_to_cns_to_reads[chr_])
                matching_reads = set()
                lasti = 0
                # Iterate through CN seg (index) pairs, step size 1
                for curr_cni, next_cni in zip(
                    sorted_cns_idxs, sorted_cns_idxs[1:]
                ):
                    nil = self.cns_intervals_by_chr[chr_][next_cni].start
                    lir = self.cns_intervals_by_chr[chr_][curr_cni].end
                    if next_cni - curr_cni > 2 or nil - lir > self.max_seq_len:
                        matching_reads |= chr_to_cns_to_reads[chr_][curr_cni]
                        cni_intv = Interval(
                            chr_,
                            sorted_cns_idxs[lasti],
                            curr_cni,
                        )
                        cni_intv_to_reads[cni_intv] = matching_reads
                        lasti = next_cni
                        matching_reads = set([])
                    else:
                        matching_reads |= chr_to_cns_to_reads[chr_][curr_cni]
                cni_intv = Interval(
                    chr_, sorted_cns_idxs[lasti], sorted_cns_idxs[-1]
                )
                matching_reads |= chr_to_cns_to_reads[chr_][sorted_cns_idxs[-1]]
                cni_intv_to_reads[cni_intv] = matching_reads

                # Refine initial intervals
                for cni_intv, matching_reads in cni_intv_to_reads.items():
                    chr_tag, start_cns, end_cns = (
                        cni_intv.chr_tag,
                        cni_intv.start,
                        cni_intv.end,
                    )
                    ns = self.cns_intervals_by_chr[chr_tag][start_cns].start
                    ne = self.cns_intervals_by_chr[chr_tag][end_cns].end
                    logger.debug(f"\t\tRefining new interval {[chr_, ns, ne]}.")
                    new_bp_list = []

                    for read_name in matching_reads:
                        new_bp_list += alignment2bp(
                            read_name,
                            self.chimeric_alignments[read_name],
                            self.min_bp_match_cutoff_,
                            20,
                            Interval(chr_tag, ns, ne),
                            interval,
                            max_nm=(
                                self.nm_stats[0] + 3 * self.nm_stats[1]
                                if self.nm_filter
                                else None
                            ),
                        )
                    logger.debug(
                        f"\t\tFound {len(new_bp_list)} reads connecting the two intervals."
                    )

                    new_bp_clusters = cluster_bp_list(
                        new_bp_list,
                        self.min_cluster_cutoff,
                        self.max_breakpoint_distance_cutoff,
                    )
                    logger.debug(
                        f"\t\tThese reads formed {len(new_bp_clusters)} clusters."
                    )
                    new_bp_refined = []
                    for c in new_bp_clusters:
                        logger.debug(f"\t\t\tNew cluster of size {len(c)}.")
                        bp_cluster_r = c
                        if len(c) >= self.min_cluster_cutoff:
                            num_subcluster = 0
                            while len(bp_cluster_r) >= self.min_cluster_cutoff:
                                logger.debug(
                                    f"\t\t\t\tSubcluster {num_subcluster}"
                                )
                                num_subcluster += 1  # TODO: this isn't in old code? just a pointless variable?
                                bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(
                                    bp_cluster_r,
                                    self.min_bp_match_cutoff_,
                                )
                                logger.debug(f"\t\t\t\t\tbp = {bp}")
                                logger.debug(
                                    f"\t\t\t\t\tNum long read support = {len(set(bpr))}"
                                )
                                logger.debug(
                                    f"\t\t\t\t\tbp_stats = {(bp_stats_)}"
                                )
                                if (
                                    len(set(bpr)) >= self.min_cluster_cutoff
                                ) or (
                                    len(set(bpr))
                                    >= max(
                                        self.normal_cov
                                        * self.min_bp_cov_factor,
                                        3.0,
                                    )
                                ):
                                    bpi = self.addbp(
                                        bp, set(bpr), bp_stats_, ccid
                                    )
                                    if bpi not in new_bp_refined:
                                        new_bp_refined.append(bpi)
                                    logger.debug("\t\t\t\tKept the cluster.")
                                else:
                                    logger.debug(
                                        "\t\t\t\tDiscarded the cluster."
                                    )
                            # logger.debug("\t\t\tDiscarded the cluster.") # TODO: confirm this is misleading?

                    nint_segs = []
                    nint_segs_ = []
                    if len(new_bp_refined) > 0:
                        self.refine_breakpoints(
                            self.amplicon_intervals[ai_],
                            cni_intv,
                            new_bp_refined,
                        )
                        nint_segs = sorted(
                            nint_segs, key=lambda item: (item[0], item[1])
                        )
                        nint_segs_ = sorted(
                            nint_segs_,
                            key=lambda item: (
                                CHR_TAG_TO_IDX[item[0]],
                                item[1],
                                item[2],
                            ),
                        )
                        lasti = 0
                        for cns_idx in range(len(nint_segs) - 1):
                            nil = self.cns_intervals_by_chr[chr_][
                                nint_segs[cns_idx + 1][0]
                            ][1]
                            ncn = self.cns_intervals_by_chr[chr_][
                                nint_segs[cns_idx + 1][0]
                            ][3]
                            lir = self.cns_intervals_by_chr[chr_][
                                nint_segs[cns_idx][0]
                            ][2]
                            lcn = self.cns_intervals_by_chr[chr_][
                                nint_segs[cns_idx][0]
                            ][3]
                            amp_flag = (ncn >= self.cn_gain) or (
                                lcn >= self.cn_gain
                            )
                            if (
                                (
                                    nint_segs[cns_idx + 1][0]
                                    - nint_segs[cns_idx][0]
                                    > 2
                                )
                                or (nil - lir > self.max_seq_len / 2)
                                or (
                                    nint_segs[cns_idx + 1][1]
                                    - nint_segs[cns_idx][1]
                                    > self.max_seq_len
                                )
                                or (
                                    not amp_flag
                                    and nil - lir > 2 * self.interval_delta
                                )
                                or (
                                    not amp_flag
                                    and nint_segs[cns_idx + 1][1]
                                    - nint_segs[cns_idx][1]
                                    > 3 * self.interval_delta
                                )
                            ):
                                amp_flag_l = (
                                    self.cns_intervals_by_chr[chr_][
                                        nint_segs[lasti][0]
                                    ][3]
                                    >= self.cn_gain
                                )
                                amp_flag_r = (
                                    self.cns_intervals_by_chr[chr_][
                                        nint_segs[cns_idx][0]
                                    ][3]
                                    >= self.cn_gain
                                )
                                if not amp_flag_l:
                                    l = max(
                                        nint_segs[lasti][1]
                                        - self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][0][1],
                                    )
                                else:
                                    l = max(
                                        self.cns_intervals_by_chr[chr_][
                                            nint_segs[lasti][0]
                                        ][1]
                                        - self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][0][1],
                                    )
                                if not amp_flag_r:
                                    read_intv = min(
                                        nint_segs[cns_idx][1]
                                        + self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][-1][2],
                                    )
                                else:
                                    read_intv = min(
                                        lir + self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][-1][2],
                                    )
                                if (
                                    self.cns_intervals_by_chr[chr_][
                                        nint_segs[lasti][0]
                                    ][3]
                                    and nint_segs[lasti][1]
                                    - int(self.max_seq_len / 2)
                                    > l
                                ):
                                    l = nint_segs[lasti][1] - int(
                                        self.max_seq_len / 2
                                    )
                                read_intv = min(
                                    nint_segs[cns_idx][1]
                                    + int(self.max_seq_len / 2),
                                    read_intv,
                                )
                                if len(list(self.pos2cni(chr_, l))) == 0:
                                    l = self.cns_intervals_by_chr[chr_][
                                        nint_segs[lasti][0]
                                    ][1]
                                if (
                                    len(list(self.pos2cni(chr_, read_intv)))
                                    == 0
                                ):
                                    read_intv = lir
                                new_intervals_refined.append(
                                    [chr_, l, read_intv, -1]
                                )
                                new_intervals_connections.append([])
                                for i_ in range(lasti, cns_idx + 1):
                                    new_intervals_connections[-1].append(
                                        nint_segs[i_][2]
                                    )
                                lasti = cns_idx + 1
                                logger.debug(
                                    "\t\tFixed new interval: %s."
                                    % new_intervals_refined[-1],
                                )
                                logger.debug(
                                    "\t\tList of breakpoints connected to the new interval:",
                                )
                                for bpi in new_intervals_connections[-1]:
                                    logger.debug(
                                        "\t\t\t%s" % self.new_bp_list[bpi][:6],
                                    )
                        if len(nint_segs) > 0:
                            amp_flag_l = (
                                self.cns_intervals_by_chr[chr_][
                                    nint_segs[lasti][0]
                                ][3]
                                >= self.cn_gain
                            )
                            amp_flag_r = (
                                self.cns_intervals_by_chr[chr_][
                                    nint_segs[-1][0]
                                ][3]
                                >= self.cn_gain
                            )
                            if not amp_flag_l:
                                l = max(
                                    nint_segs[lasti][1] - self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][0][1],
                                )
                            else:
                                l = max(
                                    self.cns_intervals_by_chr[chr_][
                                        nint_segs[lasti][0]
                                    ][1]
                                    - self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][0][1],
                                )
                            if not amp_flag_r:
                                read_intv = min(
                                    nint_segs[-1][1] + self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][-1][2],
                                )
                            else:
                                read_intv = min(
                                    self.cns_intervals_by_chr[chr_][
                                        nint_segs[-1][0]
                                    ][2]
                                    + self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][-1][2],
                                )
                            if (
                                nint_segs[lasti][1] - int(self.max_seq_len / 2)
                                > l
                            ):
                                l = (
                                    nint_segs[lasti][1]
                                    - int(self.max_seq_len / 2)
                                    > l
                                )
                            read_intv = min(
                                nint_segs[-1][1] + int(self.max_seq_len / 2),
                                read_intv,
                            )
                            if len(list(self.pos2cni(chr_, l))) == 0:
                                l = self.cns_intervals_by_chr[chr_][
                                    nint_segs[lasti][0]
                                ][1]
                            if len(list(self.pos2cni(chr_, read_intv))) == 0:
                                read_intv = self.cns_intervals_by_chr[chr_][
                                    nint_segs[-1][0]
                                ][2]
                            new_intervals_refined.append(
                                [chr_, l, read_intv, -1]
                            )
                            new_intervals_connections.append([])
                            for i_ in range(lasti, len(nint_segs)):
                                new_intervals_connections[-1].append(
                                    nint_segs[i_][2]
                                )
                            logger.debug(
                                "\t\tFixed new interval: %s."
                                % new_intervals_refined[-1],
                            )
                            logger.debug(
                                "\t\tList of breakpoints connected to the new interval:",
                            )
                            for bpi in new_intervals_connections[-1]:
                                logger.debug(
                                    "\t\t\t%s" % self.new_bp_list[bpi][:6],
                                )
                        lasti = 0
                        for cns_idx in range(len(nint_segs_) - 1):
                            # two intervals in nint_segs_ might be on different chrs
                            nil = self.cns_intervals_by_chr[
                                nint_segs_[cns_idx + 1][0]
                            ][nint_segs_[cns_idx + 1][1]][1]
                            ncn = self.cns_intervals_by_chr[
                                nint_segs_[cns_idx + 1][0]
                            ][nint_segs_[cns_idx + 1][1]][3]
                            lir = self.cns_intervals_by_chr[
                                nint_segs_[cns_idx][0]
                            ][nint_segs_[cns_idx][1]][2]
                            lcn = self.cns_intervals_by_chr[
                                nint_segs_[cns_idx][0]
                            ][nint_segs_[cns_idx][1]][3]
                            amp_flag = (ncn >= self.cn_gain) or (
                                lcn >= self.cn_gain
                            )
                            if (
                                (
                                    nint_segs_[cns_idx + 1][0]
                                    != nint_segs_[cns_idx][0]
                                )
                                or (
                                    nint_segs_[cns_idx + 1][1]
                                    - nint_segs_[cns_idx][1]
                                    > 2
                                )
                                or (nil - lir > self.max_seq_len / 2)
                                or (
                                    nint_segs_[cns_idx + 1][2]
                                    - nint_segs_[cns_idx][2]
                                    > self.max_seq_len
                                )
                                or (
                                    not amp_flag
                                    and nil - lir > 2 * self.interval_delta
                                )
                                or (
                                    not amp_flag
                                    and nint_segs_[cns_idx + 1][2]
                                    - nint_segs_[cns_idx][2]
                                    > 3 * self.interval_delta
                                )
                            ):
                                amp_flag_l = (
                                    self.cns_intervals_by_chr[
                                        nint_segs_[lasti][0]
                                    ][nint_segs_[lasti][1]][3]
                                    >= self.cn_gain
                                )
                                amp_flag_r = (
                                    self.cns_intervals_by_chr[
                                        nint_segs_[cns_idx][0]
                                    ][nint_segs_[cns_idx][1]][3]
                                    >= self.cn_gain
                                )
                                if not amp_flag_l:
                                    l = max(
                                        nint_segs_[lasti][2]
                                        - self.interval_delta,
                                        self.cns_intervals_by_chr[
                                            nint_segs_[lasti][0]
                                        ][0][1],
                                    )
                                else:
                                    l = max(
                                        self.cns_intervals_by_chr[
                                            nint_segs_[lasti][0]
                                        ][nint_segs_[lasti][1]][1]
                                        - self.interval_delta,
                                        self.cns_intervals_by_chr[
                                            nint_segs_[lasti][0]
                                        ][0][1],
                                    )
                                if not amp_flag_r:
                                    read_intv = min(
                                        nint_segs_[cns_idx][2]
                                        + self.interval_delta,
                                        self.cns_intervals_by_chr[
                                            nint_segs_[cns_idx][0]
                                        ][-1][2],
                                    )
                                else:
                                    read_intv = min(
                                        lir + self.interval_delta,
                                        self.cns_intervals_by_chr[
                                            nint_segs_[cns_idx][0]
                                        ][-1][2],
                                    )
                                l = max(
                                    nint_segs_[lasti][2]
                                    - int(self.max_seq_len / 2),
                                    l,
                                )
                                read_intv = min(
                                    nint_segs_[cns_idx][2]
                                    + int(self.max_seq_len / 2),
                                    read_intv,
                                )
                                if (
                                    len(
                                        list(
                                            self.pos2cni(
                                                nint_segs_[lasti][0], l
                                            )
                                        )
                                    )
                                    == 0
                                ):
                                    l = self.cns_intervals_by_chr[
                                        nint_segs_[lasti][0]
                                    ][nint_segs_[lasti][1]][1]
                                if (
                                    len(
                                        list(
                                            self.pos2cni(
                                                nint_segs_[cns_idx][0],
                                                read_intv,
                                            )
                                        )
                                    )
                                    == 0
                                ):
                                    read_intv = lir
                                new_intervals_refined.append(
                                    [nint_segs_[lasti][0], l, read_intv, -1]
                                )
                                new_intervals_connections.append([])
                                lasti = cns_idx + 1
                                logger.debug(
                                    "\t\tFixed new interval: %s."
                                    % new_intervals_refined[-1],
                                )
                                logger.debug(
                                    "\t\tSkip breakpoints connected to the new interval.",
                                )
                        if len(nint_segs_) > 0:
                            amp_flag_l = (
                                self.cns_intervals_by_chr[nint_segs_[lasti][0]][
                                    nint_segs_[lasti][1]
                                ][3]
                                >= self.cn_gain
                            )
                            amp_flag_r = (
                                self.cns_intervals_by_chr[nint_segs_[-1][0]][
                                    nint_segs_[-1][1]
                                ][3]
                                >= self.cn_gain
                            )
                            if not amp_flag_l:
                                l = max(
                                    nint_segs_[lasti][2] - self.interval_delta,
                                    self.cns_intervals_by_chr[
                                        nint_segs_[lasti][0]
                                    ][0][1],
                                )
                            else:
                                l = max(
                                    self.cns_intervals_by_chr[
                                        nint_segs_[lasti][0]
                                    ][nint_segs_[lasti][1]][1]
                                    - self.interval_delta,
                                    self.cns_intervals_by_chr[
                                        nint_segs_[lasti][0]
                                    ][0][1],
                                )
                            if not amp_flag_r:
                                read_intv = min(
                                    nint_segs_[-1][2] + self.interval_delta,
                                    self.cns_intervals_by_chr[
                                        nint_segs_[-1][0]
                                    ][-1][2],
                                )
                            else:
                                read_intv = min(
                                    self.cns_intervals_by_chr[
                                        nint_segs_[-1][0]
                                    ][nint_segs_[-1][1]][2]
                                    + self.interval_delta,
                                    self.cns_intervals_by_chr[
                                        nint_segs_[-1][0]
                                    ][-1][2],
                                )
                            l = max(
                                nint_segs_[lasti][2]
                                - int(self.max_seq_len / 2),
                                l,
                            )
                            read_intv = min(
                                nint_segs_[-1][2] + int(self.max_seq_len / 2),
                                read_intv,
                            )
                            if (
                                len(list(self.pos2cni(nint_segs_[lasti][0], l)))
                                == 0
                            ):
                                l = self.cns_intervals_by_chr[
                                    nint_segs_[lasti][0]
                                ][nint_segs_[lasti][1]][1]
                            if (
                                len(
                                    list(
                                        self.pos2cni(
                                            nint_segs_[lasti][0], read_intv
                                        )
                                    )
                                )
                                == 0
                            ):
                                read_intv = self.cns_intervals_by_chr[
                                    nint_segs_[lasti][0]
                                ][nint_segs_[-1][1]][2]
                            new_intervals_refined.append(
                                [nint_segs_[lasti][0], l, read_intv, -1]
                            )
                            new_intervals_connections.append([])
                            logger.debug(
                                "\t\tFixed new interval: %s."
                                % new_intervals_refined[-1],
                            )
                            logger.debug(
                                "\t\tSkip breakpoints connected to the new interval:",
                            )

            # BFS
            logger.debug("\t\tProcessing new intervals.")
            for ni in range(len(new_intervals_refined)):
                ei, intl = interval_exclusive(
                    new_intervals_refined[ni], self.amplicon_intervals
                )
                if len(intl) == 0:
                    ei_str = ""
                    for ei_ in ei:
                        ei_str += "%s " % self.amplicon_intervals[ei_]
                    ei_str = ei_str.rstrip()
                    logger.debug(
                        "\t\tNew interval %s overlaps with existing interval %s."
                        % (new_intervals_refined[ni], ei_str),
                    )
                    for bpi in new_intervals_connections[ni]:
                        bp = self.new_bp_list[bpi][:6]
                        for ei_ in ei:
                            connection = (min(ai_, ei_), max(ai_, ei_))
                            if (
                                ei_ != ai_
                                and interval_overlap(
                                    [bp[0], bp[1], bp[1]],
                                    self.amplicon_intervals[ei_],
                                )
                                or interval_overlap(
                                    [bp[3], bp[4], bp[4]],
                                    self.amplicon_intervals[ei_],
                                )
                            ):
                                try:
                                    self.amplicon_interval_connections[
                                        connection
                                    ].add(bpi)
                                except:
                                    self.amplicon_interval_connections[
                                        connection
                                    ] = set([bpi])
                    for ei_ in ei:
                        if ei_ != ai_ and self.amplicon_intervals[ei_][3] < 0:
                            interval_queue.append(ei_)
                else:
                    for int_ in intl:
                        nai = len(self.amplicon_intervals)
                        self.amplicon_intervals.append(int_)
                        logger.debug(
                            "\t\tAdded new interval %s to the amplicon interval list."
                            % int_,
                        )
                        logger.debug(
                            "\t\tNew interval index: %d." % nai,
                        )
                        self.amplicon_interval_connections[(ai_, nai)] = set([])
                        if len(ei) == 0:
                            for bpi in new_intervals_connections[ni]:
                                self.amplicon_interval_connections[
                                    (ai_, nai)
                                ].add(bpi)
                        else:
                            for bpi in new_intervals_connections[ni]:
                                bp = self.new_bp_list[bpi][:6]
                                for ei_ in ei:
                                    connection = (min(ai_, ei_), max(ai_, ei_))
                                    if interval_overlap(
                                        [bp[0], bp[1], bp[1]],
                                        self.amplicon_intervals[ei_],
                                    ) or interval_overlap(
                                        [bp[3], bp[4], bp[4]],
                                        self.amplicon_intervals[ei_],
                                    ):
                                        try:
                                            self.amplicon_interval_connections[
                                                connection
                                            ].add(bpi)
                                        except:
                                            self.amplicon_interval_connections[
                                                connection
                                            ] = set(
                                                [bpi],
                                            )
                                    else:
                                        self.amplicon_interval_connections[
                                            (ai_, nai)
                                        ].add(bpi)
                        interval_queue.append(nai)

    def find_breakpoints(self):
        """Search for breakpoints from chimeric alignments within identified amplified intervals
        Then cluster the breakpoints from chimeric alignments
        """
        new_bp_list_ = []
        for r in self.chimeric_alignments:
            new_bp_list_ += alignment2bp_l(
                r,
                self.chimeric_alignments[r],
                self.min_bp_match_cutoff_,
                20,
                100,
                self.amplicon_intervals,
                max_nm=(self.nm_stats[0] + 3 * self.nm_stats[1])
                if self.nm_filter
                else None,
            )

        logger.debug(f"Found {len(new_bp_list_)} reads with new breakpoints.")

        new_bp_clusters = cluster_bp_list(
            new_bp_list_,
            self.min_cluster_cutoff,
            self.max_breakpoint_distance_cutoff,
        )
        logger.debug(f"These reads formed {(len(new_bp_clusters))} clusters.")

        for c in new_bp_clusters:
            logger.debug("New cluster of size %d." % (len(c)))
            if len(c) >= self.min_cluster_cutoff:
                num_subcluster = 0
                bp_cluster_r = c
                while len(bp_cluster_r) >= self.min_cluster_cutoff:
                    bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(
                        bp_cluster_r,
                        self.min_bp_match_cutoff_,
                    )
                    logger.debug(
                        "\tSubcluster %d" % (num_subcluster),
                    )
                    logger.debug(
                        "\tbp = %s" % (bp),
                    )
                    logger.debug(
                        "\tNum long read support = %d" % (len(set(bpr))),
                    )
                    logger.debug(
                        "\tbp_stats = %s" % (bp_stats_),
                    )
                    if (
                        num_subcluster == 0
                        and len(set(bpr)) >= self.min_cluster_cutoff
                    ) or (
                        len(set(bpr))
                        >= max(self.normal_cov * self.min_bp_cov_factor, 3.0)
                    ):
                        io1 = interval_overlap_l(
                            [bp[0], bp[1], bp[1]], self.amplicon_intervals
                        )
                        io2 = interval_overlap_l(
                            [bp[3], bp[4], bp[4]], self.amplicon_intervals
                        )
                        if io1 >= 0 and io2 >= 0:
                            assert (
                                self.amplicon_intervals[io1][3]
                                == self.amplicon_intervals[io2][3]
                            )
                            bpi = self.addbp(
                                bp,
                                set(bpr),
                                bp_stats_,
                                self.amplicon_intervals[io1][3],
                            )
                            try:
                                self.amplicon_interval_connections[
                                    (min(io1, io2), max(io1, io2))
                                ].add(bpi)
                            except:
                                self.amplicon_interval_connections[
                                    (min(io1, io2), max(io1, io2))
                                ] = set([bpi])
                    else:
                        logger.debug(
                            f"\tDiscarded the subcluster {num_subcluster}."
                        )
                    num_subcluster += 1
            else:
                logger.debug("\tDiscarded the cluster.")

    def find_smalldel_breakpoints(self):
        """Search for breakpoints from a single alignment record, within identified amplified intervals
        For each alignment record with large indels, first try to match the resulting breakpoint to the list of AA breakpoints
        Cluster the unmatched breakpoints
        """
        new_bp_list_ = []
        if self.nm_filter:
            for ai in self.amplicon_intervals:
                for read in self.lr_bamfh.fetch(ai[0], ai[1], ai[2] + 1):
                    rn = read.query_name
                    rq = read.mapping_quality
                    rnm = read.get_cigar_stats()[0][-1]
                    if rq < 20:
                        continue
                    blocks = read.get_blocks()
                    agg_del = 0
                    for bi in range(len(blocks) - 1):
                        if (
                            abs(blocks[bi + 1][0] - blocks[bi][1])
                            > self.min_del_len
                        ):
                            agg_del += abs(blocks[bi + 1][0] - blocks[bi][1])
                    rnm = (rnm - agg_del) / read.query_length
                    if rnm < self.nm_stats[0] + 3 * self.nm_stats[1]:
                        for bi in range(len(blocks) - 1):
                            if (
                                abs(blocks[bi + 1][0] - blocks[bi][1])
                                > self.min_del_len
                            ):
                                try:
                                    self.large_indel_alignments[rn].append(
                                        [
                                            ai[0],
                                            blocks[bi + 1][0],
                                            blocks[bi][1],
                                            blocks[0][0],
                                            blocks[-1][1],
                                            rq,
                                        ],
                                    )
                                except:
                                    self.large_indel_alignments[rn] = [
                                        [
                                            ai[0],
                                            blocks[bi + 1][0],
                                            blocks[bi][1],
                                            blocks[0][0],
                                            blocks[-1][1],
                                            rq,
                                        ],
                                    ]
        else:
            for ai in self.amplicon_intervals:
                for read in self.lr_bamfh.fetch(ai[0], ai[1], ai[2] + 1):
                    rn = read.query_name
                    rq = read.mapping_quality
                    if rq < 20:
                        continue
                    blocks = read.get_blocks()
                    for bi in range(len(blocks) - 1):
                        if (
                            abs(blocks[bi + 1][0] - blocks[bi][1])
                            > self.min_del_len
                        ):
                            try:
                                self.large_indel_alignments[rn].append(
                                    [
                                        ai[0],
                                        blocks[bi + 1][0],
                                        blocks[bi][1],
                                        blocks[0][0],
                                        blocks[-1][1],
                                        rq,
                                    ],
                                )
                            except:
                                self.large_indel_alignments[rn] = [
                                    [
                                        ai[0],
                                        blocks[bi + 1][0],
                                        blocks[bi][1],
                                        blocks[0][0],
                                        blocks[-1][1],
                                        rq,
                                    ],
                                ]
        logger.info(
            "Fetched %d reads with large indels in CIGAR."
            % (len(self.large_indel_alignments)),
        )

        for r in self.large_indel_alignments.keys():
            for rr_gap_i in range(len(self.large_indel_alignments[r])):
                rr_gap = self.large_indel_alignments[r][rr_gap_i][:3]
                rr_gap_ = rr_gap
                if rr_gap[2] > rr_gap[1]:
                    rr_gap_[2] = rr_gap[1]
                    rr_gap_[1] = rr_gap[2]
                new_bp_list_.append(
                    [
                        rr_gap_[0],
                        rr_gap_[1],
                        "-",
                        rr_gap_[0],
                        rr_gap_[2],
                        "+",
                        (r, rr_gap_i, rr_gap_i),
                        0,
                        0,
                        -1,
                        -1,
                    ],
                )
        logger.debug(
            "Found %d reads with new small del breakpoints."
            % (len(new_bp_list_)),
        )

        new_bp_clusters = cluster_bp_list(
            new_bp_list_,
            self.min_cluster_cutoff,
            self.max_breakpoint_distance_cutoff,
        )
        logger.debug(
            "These reads formed %d clusters." % (len(new_bp_clusters)),
        )
        for c in new_bp_clusters:
            logger.debug("New cluster of size %d." % (len(c)))
            if len(c) >= self.min_cluster_cutoff:
                num_subcluster = 0
                bp_cluster_r = c
                while len(bp_cluster_r) >= self.min_cluster_cutoff:
                    bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(
                        bp_cluster_r,
                        self.min_bp_match_cutoff_,
                    )
                    logger.debug(
                        "\tSubcluster %d" % (num_subcluster),
                    )
                    logger.debug(
                        "\tbp = %s" % (bp),
                    )
                    logger.debug(
                        "\tNum long read support = %d" % (len(set(bpr))),
                    )
                    logger.debug(
                        "\tbp_stats = %s" % (bp_stats_),
                    )
                    if (
                        num_subcluster == 0
                        and len(set(bpr)) >= self.min_cluster_cutoff
                    ) or (
                        len(set(bpr))
                        >= max(self.normal_cov * self.min_bp_cov_factor, 3.0)
                    ):
                        io1 = interval_overlap_l(
                            [bp[0], bp[1], bp[1]], self.amplicon_intervals
                        )
                        io2 = interval_overlap_l(
                            [bp[3], bp[4], bp[4]], self.amplicon_intervals
                        )
                        if io1 >= 0 and io2 >= 0:
                            assert (
                                self.amplicon_intervals[io1][3]
                                == self.amplicon_intervals[io2][3]
                            )
                            bpi = self.addbp(
                                bp,
                                set(bpr),
                                bp_stats_,
                                self.amplicon_intervals[io1][3],
                            )
                            try:
                                self.amplicon_interval_connections[
                                    (min(io1, io2), max(io1, io2))
                                ].add(bpi)
                            except:
                                self.amplicon_interval_connections[
                                    (min(io1, io2), max(io1, io2))
                                ] = set([bpi])
                    else:
                        logger.debug(
                            "\tDiscarded the subcluster %d." % (num_subcluster),
                        )
                    num_subcluster += 1
            else:
                logger.debug(
                    "\tDiscarded the cluster.",
                )

    def build_graph(self):
        """Organize the identified discordant edges into a list of breakpoint graphs, stored in lr_graph
        Each graph represent a connected component of amplified intervals, i.e., amplicon
        """
        # Split amplified intervals according to discordant edges
        split_int = dict()
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            for ai in range(len(self.amplicon_intervals)):
                seg = interval
                if bp[0] == seg[0] and seg[1] < bp[1] < seg[2]:
                    if bp[2] == "+":
                        try:
                            split_int[ai].append(
                                (bp[1], bp[1] + 1, bpi, 1, "+")
                            )
                        except:
                            split_int[ai] = [(bp[1], bp[1] + 1, bpi, 1, "+")]
                    if bp[2] == "-":
                        try:
                            split_int[ai].append(
                                (bp[1] - 1, bp[1], bpi, 1, "-")
                            )
                        except:
                            split_int[ai] = [(bp[1] - 1, bp[1], bpi, 1, "-")]
                if bp[3] == seg[0] and seg[1] < bp[4] < seg[2]:
                    if bp[5] == "+":
                        try:
                            split_int[ai].append(
                                (bp[4], bp[4] + 1, bpi, 4, "+")
                            )
                        except:
                            split_int[ai] = [(bp[4], bp[4] + 1, bpi, 4, "+")]
                    if bp[5] == "-":
                        try:
                            split_int[ai].append(
                                (bp[4] - 1, bp[4], bpi, 4, "-")
                            )
                        except:
                            split_int[ai] = [(bp[4] - 1, bp[4], bpi, 4, "-")]

        # Split amplified intervals according to source edges
        for srci in range(len(self.source_edges)):
            srce = self.source_edges[srci]
            for ai in range(len(self.amplicon_intervals)):
                seg = interval
                if srce[3] == seg[0] and seg[1] < srce[4] < seg[2]:
                    if srce[5] == "+":
                        try:
                            split_int[ai].append(
                                (
                                    srce[4],
                                    srce[4] + 1,
                                    len(self.new_bp_list) + srci,
                                    4,
                                    "+",
                                ),
                            )
                        except:
                            split_int[ai] = [
                                (
                                    srce[4],
                                    srce[4] + 1,
                                    len(self.new_bp_list) + srci,
                                    4,
                                    "+",
                                ),
                            ]
                    if srce[5] == "-":
                        try:
                            split_int[ai].append(
                                (
                                    srce[4] - 1,
                                    srce[4],
                                    len(self.new_bp_list) + srci,
                                    4,
                                    "-",
                                ),
                            )
                        except:
                            split_int[ai] = [
                                (
                                    srce[4] - 1,
                                    srce[4],
                                    len(self.new_bp_list) + srci,
                                    4,
                                    "-",
                                ),
                            ]

        # Construct graphs with sequence and concordant edges
        logger.debug(
            "Will split the following %d amplicon intervals into sequence edges and build breakpoint graphs."
            % (len(split_int)),
        )
        amplicon_id = 1
        for sseg in self.amplicon_intervals:
            if sseg[3] not in self.ccid2id:
                self.ccid2id[sseg[3]] = amplicon_id
                amplicon_id += 1
        for cci in range(len(self.ccid2id)):
            # Initialize breakpoint graph objects
            self.lr_graph.append(BreakpointGraph())

        for ai in split_int:
            logger.debug("Will split the amplicon interval at index %d." % (ai))
            logger.debug("\tSplit interval at %s." % (split_int[ai]))
        for ai, interval in enumerate(split_int):
            split_int[ai].sort(key=lambda item: item[0])
            sseg = interval
            amplicon_idx = self.ccid2id[sseg[3]] - 1
            for ssi in range(len(split_int[ai])):
                if ssi == 0:
                    self.lr_graph[amplicon_idx].add_node(
                        (sseg[0], sseg[1], "-")
                    )
                    self.lr_graph[amplicon_idx].add_node(
                        (sseg[0], split_int[ai][ssi][0], "+")
                    )
                    self.lr_graph[amplicon_idx].add_node(
                        (sseg[0], split_int[ai][ssi][1], "-")
                    )
                    self.lr_graph[amplicon_idx].add_sequence_edge(
                        sseg[0],
                        sseg[1],
                        split_int[ai][ssi][0],
                    )
                    self.lr_graph[amplicon_idx].add_concordant_edge(
                        sseg[0],
                        split_int[ai][ssi][0],
                        "+",
                        sseg[0],
                        split_int[ai][ssi][1],
                        "-",
                    )
                elif split_int[ai][ssi][0] > split_int[ai][ssi - 1][0]:
                    self.lr_graph[amplicon_idx].add_node(
                        (sseg[0], split_int[ai][ssi - 1][1], "-")
                    )
                    self.lr_graph[amplicon_idx].add_node(
                        (sseg[0], split_int[ai][ssi][0], "+")
                    )
                    self.lr_graph[amplicon_idx].add_node(
                        (sseg[0], split_int[ai][ssi][1], "-")
                    )
                    self.lr_graph[amplicon_idx].add_sequence_edge(
                        sseg[0],
                        split_int[ai][ssi - 1][1],
                        split_int[ai][ssi][0],
                    )
                    self.lr_graph[amplicon_idx].add_concordant_edge(
                        sseg[0],
                        split_int[ai][ssi][0],
                        "+",
                        sseg[0],
                        split_int[ai][ssi][1],
                        "-",
                    )
            self.lr_graph[amplicon_idx].add_node(
                (sseg[0], split_int[ai][-1][1], "-")
            )
            self.lr_graph[amplicon_idx].add_node((sseg[0], sseg[2], "+"))
            self.lr_graph[amplicon_idx].add_sequence_edge(
                sseg[0], split_int[ai][-1][1], sseg[2]
            )
        for ai in range(len(self.amplicon_intervals)):
            if ai not in split_int:
                sseg = interval
                amplicon_idx = self.ccid2id[sseg[3]] - 1
                self.lr_graph[amplicon_idx].add_node((sseg[0], sseg[1], "-"))
                self.lr_graph[amplicon_idx].add_node((sseg[0], sseg[2], "+"))
                self.lr_graph[amplicon_idx].add_sequence_edge(
                    sseg[0], sseg[1], sseg[2]
                )
        for amplicon_idx in range(len(self.lr_graph)):
            self.lr_graph[amplicon_idx].sort_edges()

        # Add nodes corresponding to interval ends
        for ai in range(len(self.amplicon_intervals)):
            sseg = interval
            amplicon_idx = self.ccid2id[sseg[3]] - 1
            self.lr_graph[amplicon_idx].amplicon_intervals.append(
                [sseg[0], sseg[1], sseg[2]]
            )
            self.lr_graph[amplicon_idx].add_endnode((sseg[0], sseg[1], "-"))
            self.lr_graph[amplicon_idx].add_endnode((sseg[0], sseg[2], "+"))
        logger.debug(
            "The following nodes correspond to interval ends.",
        )
        for amplicon_idx in range(len(self.lr_graph)):
            for node in self.lr_graph[amplicon_idx].endnodes.keys():
                logger.debug(
                    "\tAmplicon %d, node %s." % (amplicon_idx + 1, str(node)),
                )

        # Construct graphs with discordant and source edges
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            bp_ccid = self.new_bp_ccids[bpi]
            io1 = interval_overlap_l(
                [bp[0], bp[1], bp[1]], self.amplicon_intervals
            )
            io2 = interval_overlap_l(
                [bp[3], bp[4], bp[4]], self.amplicon_intervals
            )
            assert (
                self.amplicon_intervals[io1][3]
                == self.amplicon_intervals[io2][3]
            )
            amplicon_idx = self.ccid2id[self.amplicon_intervals[io1][3]] - 1
            if self.amplicon_intervals[io1][3] != bp_ccid:
                logger.debug(
                    "Reset the ccid for breakpoint %s at index %d from %d to %d."
                    % (bp[:6], bpi, bp_ccid, self.amplicon_intervals[io1][3]),
                )
                self.new_bp_ccids[bpi] = self.amplicon_intervals[io1][3]
            self.lr_graph[amplicon_idx].add_discordant_edge(
                bp[0],
                bp[1],
                bp[2],
                bp[3],
                bp[4],
                bp[5],
                lr_count=len(bp[-1]),
                reads=bp[-1],
            )
        for srci in range(len(self.source_edges)):
            srce = self.source_edges[srci]
            src_ccid = self.source_edge_ccids[srci]
            amplicon_idx = self.ccid2id[src_ccid] - 1
            self.lr_graph[amplicon_idx].add_source_edge(
                srce[3], srce[4], srce[5]
            )

        # Print summary statistics for each amplicon
        for amplicon_idx in range(len(self.lr_graph)):
            logger.debug(
                f"Num sequence edges in amplicon {amplicon_idx + 1} = {len(self.lr_graph[amplicon_idx].sequence_edges)}."
            )
            logger.debug(
                f"Num concordant edges in amplicon {amplicon_idx + 1} = {len(self.lr_graph[amplicon_idx].concordant_edges)}."
            )
            logger.debug(
                f"Num discordant edges in amplicon {amplicon_idx + 1} = {len(self.lr_graph[amplicon_idx].discordant_edges)}."
            )
            logger.debug(
                f"Num source edges in amplicon {amplicon_idx + 1} = {len(self.lr_graph[amplicon_idx].source_edges)}."
            )

    def assign_cov(self):
        """Extract the long read coverage from bam file, if missing, for each sequence edge"""
        for amplicon_idx in range(len(self.lr_graph)):
            for seqi in range(len(self.lr_graph[amplicon_idx].sequence_edges)):
                seg = self.lr_graph[amplicon_idx].sequence_edges[seqi]
                if seg[5] == -1:
                    logger.debug(f"Finding LR cov for sequence edge {seg[:3]}.")
                    """
					For long read, use the total number of nucleotides
					"""
                    rl_list = [
                        read
                        for read in self.lr_bamfh.fetch(
                            seg[0], seg[1], seg[2] + 1
                        )
                        if read.infer_read_length()
                    ]
                    self.lr_graph[amplicon_idx].sequence_edges[seqi][5] = len(
                        rl_list
                    )
                    self.lr_graph[amplicon_idx].sequence_edges[seqi][6] = sum(
                        [
                            sum(nc)
                            for nc in self.lr_bamfh.count_coverage(
                                seg[0],
                                seg[1],
                                seg[2] + 1,
                                quality_threshold=0,
                                read_callback="nofilter",
                            )
                        ],
                    )
                    logger.debug(
                        f"LR cov assigned for sequence edge {seg[:3]}."
                    )
        """
		Extract the long read coverage from bam file, if missing, for each concordant edge 
		"""
        for amplicon_idx in range(len(self.lr_graph)):
            for eci in range(len(self.lr_graph[amplicon_idx].concordant_edges)):
                ec = self.lr_graph[amplicon_idx].concordant_edges[eci]
                logger.debug(f"Finding cov for concordant edge {ec[:6]}.")
                rls = set(
                    [
                        read.query_name
                        for read in self.lr_bamfh.fetch(
                            contig=ec[0], start=ec[1], stop=ec[1] + 1
                        )
                    ],
                )
                rrs = set(
                    [
                        read.query_name
                        for read in self.lr_bamfh.fetch(
                            contig=ec[3], start=ec[4], stop=ec[4] + 1
                        )
                    ],
                )
                rls1 = set(
                    [
                        read.query_name
                        for read in self.lr_bamfh.fetch(
                            contig=ec[0],
                            start=ec[1] - self.min_bp_match_cutoff_ - 1,
                            stop=ec[1] - self.min_bp_match_cutoff_,
                        )
                    ],
                )
                rrs1 = set(
                    [
                        read.query_name
                        for read in self.lr_bamfh.fetch(
                            contig=ec[3],
                            start=ec[4] + self.min_bp_match_cutoff_,
                            stop=ec[4] + self.min_bp_match_cutoff_ + 1,
                        )
                    ],
                )
                rbps = set([])
                for bpi in self.lr_graph[amplicon_idx].nodes[
                    (ec[0], ec[1], ec[2])
                ][2]:
                    for r in self.lr_graph[amplicon_idx].discordant_edges[bpi][
                        10
                    ]:
                        rbps.add(r[0])
                for bpi in self.lr_graph[amplicon_idx].nodes[
                    (ec[3], ec[4], ec[5])
                ][2]:
                    for r in self.lr_graph[amplicon_idx].discordant_edges[bpi][
                        10
                    ]:
                        rbps.add(r[0])
                self.lr_graph[amplicon_idx].concordant_edges[eci][9] = rls | rrs
                self.lr_graph[amplicon_idx].concordant_edges[eci][8] = len(
                    (rls & rrs & rls1 & rrs1) - rbps,
                )
                logger.debug(f"LR cov assigned for concordant edge {ec[:6]}.")

    def compute_path_constraints(self):
        """Convert reads mapped within the amplicons into subpath constraints"""
        for amplicon_idx in range(len(self.lr_graph)):
            self.path_constraints[amplicon_idx] = [[], [], []]
            self.longest_path_constraints[amplicon_idx] = [[], [], []]
            bp_reads = dict()
            concordant_reads = dict()
            ld = len(self.lr_graph[amplicon_idx].discordant_edges)
            for di in range(ld):
                bp = self.lr_graph[amplicon_idx].discordant_edges[di]
                for r_ in bp[10]:
                    if r_[1] == r_[2]:
                        if r_[0] in bp_reads:
                            bp_reads[r_[0]][1].append([r_[1], r_[2], di])
                        else:
                            bp_reads[r_[0]] = [[], [[r_[1], r_[2], di]]]
                    elif r_[0] in bp_reads:
                        bp_reads[r_[0]][0].append([r_[1], r_[2], di])
                    else:
                        bp_reads[r_[0]] = [[[r_[1], r_[2], di]], []]
            logger.debug(
                f"There are {len(bp_reads)} reads covering >=1 breakpoint in amplicon {amplicon_idx + 1}."
            )
            for rn in bp_reads:
                bp_reads_rn = bp_reads[rn][0]
                bp_reads_rn_sdel = bp_reads[rn][1]
                paths = []
                if len(bp_reads_rn) == 1 and len(bp_reads_rn_sdel) == 0:
                    rints = [
                        aint[:4] for aint in self.chimeric_alignments[rn][1]
                    ]
                    ai1 = bp_reads_rn[0][0]
                    ai2 = bp_reads_rn[0][1]
                    bpi = bp_reads_rn[0][2]
                    logger.debug(f"Read {rn} covers a single breakpoint.")
                    logger.debug(
                        f"Alignment intervals on reference = {rints}; mapq = {self.chimeric_alignments[rn][2]}; bp = ({ai1}, {ai2}, {bpi})"
                    )
                    path = path_constraints.chimeric_alignment_to_path_i(
                        self.lr_graph[amplicon_idx],
                        rints,
                        ai1,
                        ai2,
                        bpi,
                    )
                    paths.append(path)
                    logger.debug(f"Resulting subpath = {path}")
                elif len(bp_reads_rn) > 1 and len(bp_reads_rn_sdel) == 0:
                    bp_reads_rn = sorted(
                        bp_reads_rn, key=lambda item: min(item[0], item[1])
                    )
                    bp_reads_rn_split = [[0]]
                    last_ai = max(bp_reads_rn[0][0], bp_reads_rn[0][1])
                    for i in range(1, len(bp_reads_rn)):
                        if min(bp_reads_rn[i][0], bp_reads_rn[i][1]) == last_ai:
                            bp_reads_rn_split[-1].append(i)
                        else:
                            bp_reads_rn_split.append([i])
                        last_ai = max(bp_reads_rn[i][0], bp_reads_rn[i][1])
                    logger.debug(f"Read {rn} covers multiple breakpoints.")
                    logger.debug(
                        f"Blocks of local alignments: {bp_reads_rn_split}"
                    )
                    qints = self.chimeric_alignments[rn][0]
                    skip = 0
                    for qi in range(len(qints) - 1):
                        if (
                            qints[qi + 1][0] - qints[qi][1]
                            < -self.min_bp_match_cutoff_
                        ):
                            skip = 1
                            break
                    if skip == 1:
                        logger.debug(
                            "Discarded the read due to overlapping local alignments."
                        )
                        logger.debug(
                            f"Alignment intervals on reference = {self.chimeric_alignments[rn][1]}; mapq = {self.chimeric_alignments[rn][2]}."
                        )
                        logger.debug(
                            f"Alignment intervals on the read = {qints}."
                        )
                        continue
                    for ai_block in bp_reads_rn_split:
                        rints = [
                            aint[:4] for aint in self.chimeric_alignments[rn][1]
                        ]
                        ai_list = [bp_reads_rn[bi][:2] for bi in ai_block]
                        bp_list = [bp_reads_rn[bi][2] for bi in ai_block]
                        if len(set(bp_list)) < len(bp_list):
                            logger.debug(
                                "\tDiscarded the block due to repeated breakpoints.",
                            )
                            logger.debug(
                                "\tBlocks of local alignments: %s" % ai_block,
                            )
                            continue
                        path = path_constraints.chimeric_alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            rints,
                            ai_list,
                            bp_list,
                        )
                        paths.append(path)
                        logger.debug(
                            "\tAlignment intervals on reference = %s; mapq = %s; bps = %s"
                            % (
                                rints,
                                self.chimeric_alignments[rn][2],
                                bp_reads_rn,
                            ),
                        )
                        logger.debug(
                            "\tResulting subpath = %s" % path,
                        )
                elif len(bp_reads_rn) == 0 and len(bp_reads_rn_sdel) == 1:
                    rints = self.large_indel_alignments[rn][0]
                    rq = rints[-1]
                    logger.debug(
                        "Read %s covers a single small del breakpoint." % rn,
                    )
                    if rints[3] < rints[4]:
                        if rints[2] < rints[1]:
                            rints = [
                                [rints[0], rints[3], rints[2], "+"],
                                [rints[0], rints[1], rints[4], "+"],
                            ]
                        else:
                            logger.debug(
                                "\tDiscarded the read due to inconsistent alignment information.",
                            )
                            continue
                    elif rints[2] > rints[1]:
                        rints = [
                            [rints[0], rints[3], rints[2], "-"],
                            [rints[0], rints[1], rints[4], "-"],
                        ]
                    else:
                        logger.debug(
                            "\tDiscarded the read due to inconsistent alignment information.",
                        )
                        continue
                    bpi = bp_reads_rn_sdel[0][2]
                    path = []
                    if rints[0][3] == "+":
                        logger.debug(
                            "\tAlignment intervals on reference = %s; mapq = %s; bp = (1, 0, %d)"
                            % (rints, rq, bpi),
                        )
                        path = path_constraints.chimeric_alignment_to_path_i(
                            self.lr_graph[amplicon_idx],
                            rints,
                            1,
                            0,
                            bpi,
                        )
                        paths.append(path)
                    else:
                        logger.debug(
                            "\tAlignment intervals on reference = %s; mapq = %s; bp = (0, 1, %d)"
                            % (rints, rq, bpi),
                        )
                        path = path_constraints.chimeric_alignment_to_path_i(
                            self.lr_graph[amplicon_idx],
                            rints,
                            0,
                            1,
                            bpi,
                        )
                        paths.append(path)
                    logger.debug(
                        "\tResulting subpath = %s" % path,
                    )
                elif len(bp_reads_rn) == 0 and len(bp_reads_rn_sdel) > 1:
                    rints = self.large_indel_alignments[rn]
                    rq = rints[0][-1]
                    rints_ = set(
                        [
                            (
                                rint[0],
                                min(rint[3], rint[4]),
                                max(rint[3], rint[4]),
                            )
                            for rint in rints
                        ],
                    )
                    logger.debug(
                        "Read %s covers multiple small del breakpoints." % rn,
                    )
                    if len(rints_) > 1 or len(rints) <= 1:
                        logger.debug(
                            "\tDiscarded the read due to inconsistent alignment information.",
                        )
                        continue
                    rints_ = [
                        [
                            rint[0],
                            min(rint[3], rint[4]),
                            max(rint[3], rint[4]),
                            "+",
                        ]
                        for rint in rints
                    ]
                    rints = sorted(
                        rints, key=lambda item: min(item[1], item[2])
                    )
                    for ri in range(len(rints)):
                        rint = rints[ri]
                        rints_.append(
                            [
                                rint[0],
                                min(rint[3], rint[4]),
                                max(rint[3], rint[4]),
                                "+",
                            ]
                        )
                        rints_[ri][2] = min(rint[1], rint[2])
                        rints_[ri + 1][1] = max(rint[1], rint[2])
                    bp_reads_rn_sdel_split = [[]]
                    bp_reads_rn_sdel = sorted(
                        bp_reads_rn_sdel, key=lambda item: item[0]
                    )
                    last_ai = 0
                    for i in range(len(bp_reads_rn_sdel)):
                        if i == 0 or bp_reads_rn_sdel[i][0] == last_ai + 1:
                            bp_reads_rn_sdel_split[-1].append(i)
                        else:
                            bp_reads_rn_sdel_split.append([i])
                        last_ai = bp_reads_rn_sdel[i][0]
                    logger.debug(
                        "\tBlocks of local alignments: %s"
                        % bp_reads_rn_sdel_split,
                    )
                    for ai_block in bp_reads_rn_sdel_split:
                        ai_list = [
                            [
                                bp_reads_rn_sdel[bi][0],
                                bp_reads_rn_sdel[bi][0] + 1,
                            ]
                            for bi in ai_block
                        ]
                        bp_list = [bp_reads_rn_sdel[bi][2] for bi in ai_block]
                        if len(set(bp_list)) < len(bp_list):
                            logger.debug(
                                "\tDiscarded the block due to repeated breakpoints.",
                            )
                            logger.debug(
                                "\tBlocks of local alignments: %s" % ai_block,
                            )
                            continue
                        path = path_constraints.chimeric_alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            rints_,
                            ai_list,
                            bp_list,
                        )
                        paths.append(path)
                        logger.debug(
                            f"\tAlignment intervals on reference = {rints_}; mapq = {rq}; bps = {bp_reads_rn_sdel}"
                        )
                        logger.debug(f"Resulting subpath = {path}")
                else:
                    rints = [
                        aint[:4] for aint in self.chimeric_alignments[rn][1]
                    ]
                    rints_ = self.large_indel_alignments[rn]
                    rint_split = []
                    skip = 0
                    logger.debug(
                        f"Read {rn} covers breakpoints and small del breakpoints."
                    )
                    logger.debug(
                        "\tAlignment intervals on reference = %s; mapq = %s"
                        % (rints, self.chimeric_alignments[rn][2]),
                    )
                    logger.debug(
                        "\tSmall del alignment intervals on reference = %s"
                        % rints_,
                    )
                    for rint_ in rints_:
                        fount_split_rint = 0
                        for ri in range(len(rints)):
                            rint = rints[ri]
                            if (
                                rint_[0] == rint[0]
                                and min(rint_[1], rint_[2])
                                > min(rint[1], rint[2])
                                and max(rint_[1], rint_[2])
                                < max(rint[1], rint[2])
                            ):
                                fount_split_rint = 1
                                rint_split.append(ri)
                                break
                        if fount_split_rint == 0:
                            skip = 1
                            break
                    if skip == 1:
                        logger.debug(
                            "\tDiscarded the read due to inconsistent alignment information.",
                        )
                        continue
                    for rsi in range(len(rint_split)):
                        ri = rint_split[rsi]
                        rints.insert(ri, rints[ri][:])
                        if rints[ri][3] == "+":
                            rints[ri][2] = min(rints_[rsi][1], rints_[rsi][2])
                            rints[ri + 1][1] = max(
                                rints_[rsi][1], rints_[rsi][2]
                            )
                        else:
                            rints[ri][2] = max(rints_[rsi][1], rints_[rsi][2])
                            rints[ri + 1][1] = min(
                                rints_[rsi][1], rints_[rsi][2]
                            )
                        for i in range(len(bp_reads_rn)):
                            if (
                                bp_reads_rn[i][0] >= ri
                                and bp_reads_rn[i][1] >= ri
                            ):
                                bp_reads_rn[i][0] += 1
                                bp_reads_rn[i][1] += 1
                        for i in range(len(bp_reads_rn_sdel)):
                            if bp_reads_rn_sdel[i][0] == rsi:
                                if rints[ri][3] == "+":
                                    bp_reads_rn.append(
                                        [ri + 1, ri, bp_reads_rn_sdel[i][2]]
                                    )
                                else:
                                    bp_reads_rn.append(
                                        [ri, ri + 1, bp_reads_rn_sdel[i][2]]
                                    )
                    bp_reads_rn = sorted(
                        bp_reads_rn, key=lambda item: min(item[0], item[1])
                    )
                    bp_reads_rn_split = [[0]]
                    last_ai = max(bp_reads_rn[0][0], bp_reads_rn[0][1])
                    for i in range(1, len(bp_reads_rn)):
                        if min(bp_reads_rn[i][0], bp_reads_rn[i][1]) == last_ai:
                            bp_reads_rn_split[-1].append(i)
                        else:
                            bp_reads_rn_split.append([i])
                        last_ai = max(bp_reads_rn[i][0], bp_reads_rn[i][1])
                    logger.debug(
                        "\tBlocks of local alignments: %s" % bp_reads_rn_split,
                    )
                    qints = self.chimeric_alignments[rn][0]
                    skip = 0
                    for qi in range(len(qints) - 1):
                        if (
                            qints[qi + 1][0] - qints[qi][1]
                            < -self.min_bp_match_cutoff_
                        ):
                            skip = 1
                            break
                    if skip == 1:
                        logger.debug(
                            "\tDiscarded the read due to overlapping local alignments.",
                        )
                        logger.debug(
                            "\tAlignment intervals on reference = %s; mapq = %s."
                            % (
                                self.chimeric_alignments[rn][1],
                                self.chimeric_alignments[rn][2],
                            ),
                        )
                        logger.debug(
                            "\tAlignment intervals on the read = %s." % qints,
                        )
                        continue
                    for ai_block in bp_reads_rn_split:
                        ai_list = [bp_reads_rn[bi][:2] for bi in ai_block]
                        bp_list = [bp_reads_rn[bi][2] for bi in ai_block]
                        if len(set(bp_list)) < len(bp_list):
                            logger.debug(
                                "\tDiscarded the block due to repeated breakpoints.",
                            )
                            logger.debug(
                                "\tBlocks of local alignments: %s" % ai_block,
                            )
                            continue
                        path = path_constraints.chimeric_alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            rints,
                            ai_list,
                            bp_list,
                        )
                        paths.append(path)
                        logger.debug(
                            "\tAlignment intervals on reference = %s; mapq (unsplit) = %s; bps = %s"
                            % (
                                rints,
                                self.chimeric_alignments[rn][2],
                                bp_reads_rn,
                            ),
                        )
                        logger.debug(
                            "\tResulting subpath = %s" % path,
                        )
                for pi in range(len(paths)):
                    path = paths[pi]
                    if len(path) > 5 and path_constraints.valid_path(
                        self.lr_graph[amplicon_idx], path
                    ):
                        if path in self.path_constraints[amplicon_idx][0]:
                            pci = self.path_constraints[amplicon_idx][0].index(
                                path
                            )
                            self.path_constraints[amplicon_idx][1][pci] += 1
                        elif (
                            path[::-1] in self.path_constraints[amplicon_idx][0]
                        ):
                            pci = self.path_constraints[amplicon_idx][0].index(
                                path[::-1]
                            )
                            self.path_constraints[amplicon_idx][1][pci] += 1
                        else:
                            self.path_constraints[amplicon_idx][0].append(path)
                            self.path_constraints[amplicon_idx][1].append(1)
                            self.path_constraints[amplicon_idx][2].append(
                                amplicon_idx
                            )
            logger.debug(
                "There are %d distinct subpaths due to reads involving breakpoints in amplicon %d."
                % (
                    len(self.path_constraints[amplicon_idx][0]),
                    amplicon_idx + 1,
                )
            )

            # Extract reads in concordant_edges_reads
            lc = len(self.lr_graph[amplicon_idx].concordant_edges)
            for ci in range(lc):
                for rn in self.lr_graph[amplicon_idx].concordant_edges[ci][9]:
                    if (
                        rn not in self.large_indel_alignments
                        and rn not in self.chimeric_alignments
                    ):
                        concordant_reads[rn] = amplicon_idx
            logger.debug(
                "There are %d concordant reads within amplicon intervals in amplicon %d."
                % (len(concordant_reads), amplicon_idx + 1)
            )
            for aint in self.amplicon_intervals:
                if amplicon_idx != self.ccid2id[aint[3]] - 1:
                    continue
                for read in self.lr_bamfh.fetch(aint[0], aint[1], aint[2] + 1):
                    rn = read.query_name
                    q = read.mapq
                    if q >= 20 and rn in concordant_reads:
                        path = path_constraints.alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            [
                                read.reference_name,
                                read.reference_start,
                                read.reference_end,
                            ],
                        )
                        if len(path) > 5 and path_constraints.valid_path(
                            self.lr_graph[amplicon_idx], path
                        ):
                            if path in self.path_constraints[amplicon_idx][0]:
                                pci = self.path_constraints[amplicon_idx][
                                    0
                                ].index(path)
                                self.path_constraints[amplicon_idx][1][pci] += 1
                            elif (
                                path[::-1]
                                in self.path_constraints[amplicon_idx][0]
                            ):
                                pci = self.path_constraints[amplicon_idx][
                                    0
                                ].index(path[::-1])
                                self.path_constraints[amplicon_idx][1][pci] += 1
                            else:
                                self.path_constraints[amplicon_idx][0].append(
                                    path
                                )
                                self.path_constraints[amplicon_idx][1].append(1)
                                self.path_constraints[amplicon_idx][2].append(
                                    concordant_reads[rn]
                                )
            logger.debug(
                f"There are {len(self.path_constraints[amplicon_idx][0])} distinct subpaths in amplicon {amplicon_idx + 1}."
            )

    def closebam(self):
        """Close the short read and long read bam file"""
        self.lr_bamfh.close()


def reconstruct_graph(
    lr_bam_filename: pathlib.Path,
    cnv_seed_file: typer.FileText,
    cn_seg_file: typer.FileText,
    output_dir: str,
    output_bp: bool,
    skip_cycle_decomp: bool,
    output_all_path_contraints: bool,
    min_bp_support: float,
    cycle_decomp_alpha: float,
    cycle_decomp_time_limit: int,
    cycle_decomp_threads: int,
    postprocess_greedy_sol: bool,
    log_file: str,
):
    seed_intervals = breakpoint_utilities.get_intervals_from_seed_file(
        cnv_seed_file
    )  # type: ignore[arg-type] # file is passed as IO stream
    logger.debug(f"Parsed {len(seed_intervals)} seed amplicon intervals.")
    for interval in seed_intervals:
        logger.debug(f"Seed interval: {interval}")
    bam_file = pysam.AlignmentFile(str(lr_bam_filename), "rb")
    b2bn = LongReadBamToBreakpointMetadata(
        lr_bamfh=bam_file,
        amplicon_intervals=seed_intervals,
    )
    # if filter_bp_by_edit_distance:
    # b2bn.nm_filter = True
    b2bn.min_bp_cov_factor = min_bp_support
    logger.info("Opened LR bam files.")

    b2bn.read_cns(cn_seg_file)
    logger.info("Completed parsing CN segment files.")
    chimeric_alignments, edit_dist_stats = (
        breakpoint_utilities.fetch_breakpoint_reads(bam_file)
    )
    logger.info("Completed fetching reads containing breakpoints.")
    chimeric_alignments_by_read, chr_cns_to_chimeric_alignments = (
        b2bn.hash_alignment_to_seg(chimeric_alignments)
    )
    b2bn.find_amplicon_intervals()
    logger.info("Completed finding amplicon intervals.")
    b2bn.find_smalldel_breakpoints()
    logger.info("Completed finding small del breakpoints.")
    b2bn.find_breakpoints()
    logger.info("Completed finding all discordant breakpoints.")
    if output_bp:
        b2bn.build_graph()
        logger.info("Breakpoint graph built for all amplicons.")
        for gi in range(len(b2bn.lr_graph)):
            bp_stats_i = []
            for discordant_edge in b2bn.lr_graph[gi].discordant_edges:
                for bpi in range(len(b2bn.new_bp_list)):
                    bp_ = b2bn.new_bp_list[bpi]
                    if (
                        discordant_edge[0] == bp_[0]
                        and discordant_edge[1] == bp_[1]
                        and discordant_edge[2] == bp_[2]
                        and discordant_edge[3] == bp_[3]
                        and discordant_edge[4] == bp_[4]
                        and discordant_edge[5] == bp_[5]
                    ):
                        bp_stats_i.append(b2bn.new_bp_stats[bpi])
                        break
            file_prefix = f"{output_dir}/amplicon{gi+1}"
            breakpoint_utilities.output_breakpoint_info_lr(
                b2bn.lr_graph[gi], file_prefix + "_breakpoints.txt", bp_stats_i
            )
            with open(f"{file_prefix}_bp.graph", "wb") as file:
                pickle.dump(b2bn.lr_graph[gi], file)
        # Unable to dump full metadata class until we write custom Cython `__reduce__` methods for pysam objects
        # with open(f"{output_dir}/full_bb.graph", "wb") as file:
        #     pickle.dump(b2bn, file)
        logger.info(
            f"Wrote breakpoint information, for all amplicons, to {output_dir}/amplicon*_breakpoints.txt."
        )
    else:
        b2bn.build_graph()
        logger.info("Breakpoint graph built for all amplicons.")
        b2bn.assign_cov()
        logger.info(
            "Fetched read coverage for all sequence and concordant edges."
        )
        for gi in range(len(b2bn.lr_graph)):
            b2bn.lr_graph[gi].compute_cn_lr(b2bn.normal_cov)
        logger.info("Computed CN for all edges.")
        for gi in range(len(b2bn.lr_graph)):
            breakpoint_utilities.output_breakpoint_graph_lr(
                b2bn.lr_graph[gi], f"{output_dir}/amplicon{gi+1}_graph.txt"
            )
            file_prefix = f"{output_dir}/amplicon{gi+1}"
            with open(
                f"{file_prefix}_bp.graph", "wb"
            ) as file:  # TODO: merge this logic into above fn
                pickle.dump(b2bn.lr_graph[gi], file)
        # Unable to dump full metadata class until we write custom Cython `__reduce__` methods for pysam objects
        # with open(f"{output_dir}/full_bb.graph", "wb") as file:
        #     pickle.dump(b2bn, file)
        logger.info(
            f"Wrote breakpoint graph for all complicons to {output_dir}/amplicon*_graph.txt."
        )

    return b2bn
