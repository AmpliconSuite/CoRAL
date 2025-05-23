"""Infer breakpoint graph(s) from long read alignments + associated CN calls."""

from __future__ import annotations

import copy
import functools
import io
import logging
import pathlib
import pickle
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np
import pysam
import typer

from coral import bam_types, core_types, core_utils
from coral.breakpoint import breakpoint_utilities, path_utilities
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.breakpoint_types import CNSSegData
from coral.breakpoint.breakpoint_utilities import (
    alignment2bp,
    alignment2bp_l,
    bpc2bp,
    cluster_bp_list,
    interval_exclusive,
    interval_overlap_l,
)
from coral.constants import CHR_TAG_TO_IDX
from coral.datatypes import (
    AmpliconInterval,
    BPAlignments,
    BPIndexedAlignmentContainer,
    BPIndexedAlignments,
    BPToChrCNI,
    BPToCNI,
    Breakpoint,
    BreakpointStats,
    ChimericAlignment,
    CNInterval,
    CNSIntervalTree,
    FinalizedPathConstraint,
    Interval,
    LargeIndelAlignment,
    Node,
    OptimizationWalk,
    OutputPCOptions,
    PathConstraint,
    ReferenceInterval,
    SourceEdge,
    SplitInterval,
    Strand,
    WalkData,
)
from coral.models import path_constraints

edge_type_to_index = {"s": 0, "c": 1, "d": 2}


logger = logging.getLogger(__name__)


@dataclass
class LongReadBamToBreakpointMetadata:
    lr_bamfh: pysam.AlignmentFile  # Long read bam file
    bam: bam_types.BAMWrapper

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
    # small_del_cutoff = 10000 # +- breakpoints (small deletions) with the two
    # ends less than this cutoff are treated specially

    read_length: dict[str, int] = field(
        default_factory=dict
    )  # Map read name -> read length

    # Map read name -> chimeric alignments (two or more records for one read)
    chimeric_alignments: dict[str, list[ChimericAlignment]] = field(
        default_factory=lambda: defaultdict(list)
    )
    # Map chr tag -> CN segment idx -> reads
    chr_to_cni_to_chimeric_ref_intvs: defaultdict[
        str, defaultdict[int, list[ReferenceInterval]]
    ] = field(default_factory=lambda: defaultdict(lambda: defaultdict(list)))

    # Map read name -> alignments with one record per read but large indels
    # showing in CIGAR string
    large_indel_alignments: dict[str, list[LargeIndelAlignment]] = field(
        default_factory=lambda: defaultdict(list)
    )
    # For edit distance filter of breakpoints
    nm_stats: list[float] = field(
        default_factory=functools.partial(lambda: [0.0] * 3)
    )  # Default of [0.0, 0.0, 0.0]
    nm_filter: bool = False

    # AA amplicon intervals
    amplicon_intervals: list[AmpliconInterval] = field(default_factory=list)
    amplicon_interval_connections: dict[tuple[int, int], set[int]] = field(
        default_factory=lambda: defaultdict(set)
    )

    cns_intervals: list[CNInterval] = field(
        default_factory=list
    )  # CN segments
    cns_intervals_by_chr: dict[str, list[CNInterval]] = field(
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
    new_bp_stats: list[BreakpointStats] = field(default_factory=list)
    new_bp_ccids: list[int] = field(default_factory=list)
    source_edges: list[SourceEdge] = field(default_factory=list)
    source_edge_ccids: list[int] = field(default_factory=list)
    path_constraints: dict[int, list[PathConstraint]] = field(
        default_factory=lambda: defaultdict(list)
    )
    longest_path_constraints: dict[int, list[FinalizedPathConstraint]] = field(
        default_factory=lambda: defaultdict(list)
    )
    walks_by_amplicon: dict[int, WalkData[OptimizationWalk]] = field(
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

        logger.info(f"Total num LR CN segments:{len(self.log2_cn)}.")
        log2_cn_order = np.argsort(self.log2_cn)
        cns_intervals_median: list[CNInterval] = []
        log2_cn_median = []
        im = int(len(log2_cn_order) / 2.4)
        ip = im + 1
        total_int_len = 0
        cns_intervals_median.append(self.cns_intervals[log2_cn_order[ip]])
        cns_intervals_median.append(self.cns_intervals[log2_cn_order[im]])
        log2_cn_median.append(self.log2_cn[log2_cn_order[ip]])
        log2_cn_median.append(self.log2_cn[log2_cn_order[im]])
        total_int_len += len(self.cns_intervals[log2_cn_order[ip]])
        total_int_len += len(self.cns_intervals[log2_cn_order[im]])
        i = 1
        while total_int_len < 10_000_000:
            cns_intervals_median.append(
                self.cns_intervals[log2_cn_order[ip + i]]
            )
            cns_intervals_median.append(
                self.cns_intervals[log2_cn_order[im - i]]
            )
            log2_cn_median.append(self.log2_cn[log2_cn_order[ip]])
            log2_cn_median.append(self.log2_cn[log2_cn_order[im]])
            total_int_len += len(self.cns_intervals[log2_cn_order[ip + i]])
            total_int_len += len(self.cns_intervals[log2_cn_order[im - i]])
            i += 1
        logger.debug(
            f"Use {len(cns_intervals_median)} LR copy number segments."
        )
        logger.debug(
            f"Total length of LR copy number segments: {total_int_len}."
        )
        logger.debug(f"Average LR copy number: {np.average(log2_cn_median)}.")
        nnc: float = 0
        for i in range(len(cns_intervals_median)):
            median_intv = cns_intervals_median[i]
            nnc += self.bam.count_raw_coverage(median_intv)
        self.normal_cov = nnc * 1.0 / total_int_len
        logger.info(
            f"LR normal cov ={self.normal_cov}, {nnc=}, {total_int_len=}."
        )
        self.min_cluster_cutoff = max(
            self.min_cluster_cutoff,
            int(self.min_bp_cov_factor * self.normal_cov),
        )
        logger.debug(f"Reset min_cluster_cutoff to {self.min_cluster_cutoff}.")

    def pos2cni(self, chr: str, pos: int) -> core_types.CNSIdx:
        return self.cns_tree[chr][pos]

    def hash_alignment_to_seg(
        self, chimeras_by_read: dict[str, list[ChimericAlignment]]
    ) -> tuple[
        dict[str, list[ChimericAlignment]],
        defaultdict[str, defaultdict[int, list[ReferenceInterval]]],
    ]:
        """Speed up amplified interval search by hashing chimeric alignments
        from each long read to CN segments."""
        chr_cns_to_chimeras: defaultdict[
            str, defaultdict[int, list[ReferenceInterval]]
        ] = defaultdict(lambda: defaultdict(list))

        for alignments in chimeras_by_read.values():
            for alignment in alignments:
                read_interval = alignment.ref_interval
                if (chr_tag := read_interval.chr) not in self.cns_tree:
                    continue
                cn_seg_idxs = self.cns_tree[chr_tag].get_cn_segment_indices(
                    read_interval
                )
                alignment.cns_idxs = cn_seg_idxs
                for cni in cn_seg_idxs:
                    chr_cns_to_chimeras[chr_tag][cni].append(read_interval)
        logger.info("Completed hashing chimeric reads to CN segments.")
        self.chr_to_cni_to_chimeric_ref_intvs = chr_cns_to_chimeras
        self.chimeric_alignments = chimeras_by_read
        return chimeras_by_read, chr_cns_to_chimeras

    def widen_seed_intervals(self) -> None:
        """Widen seed intervals to fully encompass CN segments that the
        interval falls within."""
        for interval in self.amplicon_intervals:
            chr_tag = interval.chr
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
            logger.info(f"\tWidened amplicon interval {interval}")

    def merge_amplicon_intervals(self) -> None:
        lastai = 0
        intervals_to_merge: list[tuple[int, int]] = []
        # Iterate through amplicon intervals in pairs, step size 1
        for ai, (curr, next) in enumerate(
            zip(self.amplicon_intervals, self.amplicon_intervals[1:])
        ):
            if not (curr.is_adjacent(next) or curr.does_overlap(next)):
                if ai > lastai:
                    logger.info(
                        f"Merging intervals from {lastai} to {ai} "
                        f"{self.amplicon_intervals[lastai]} and "
                        f"{self.amplicon_intervals[ai]}."
                    )
                    intervals_to_merge.append((lastai, ai))
                lastai = ai + 1
        num_amplicons = len(self.amplicon_intervals)
        if num_amplicons and lastai < num_amplicons - 1:
            logger.info(
                f"Merging intervals from {lastai} to {num_amplicons-1} "
                f"{self.amplicon_intervals[lastai]} and "
                f"{self.amplicon_intervals[num_amplicons-1]}."
            )
            intervals_to_merge.append((lastai, num_amplicons - 1))

        orig_intervals = copy.deepcopy(self.amplicon_intervals)
        original_interval_idxs: list[int] = sorted(
            range(len(self.amplicon_intervals)),
            key=lambda i: self.amplicon_intervals[i],
        )

        for ai1, ai2 in reversed(intervals_to_merge):
            start_intv, end_intv = (
                self.amplicon_intervals[ai1],
                self.amplicon_intervals[ai2],
            )
            # Reset interval
            prev_intv = copy.deepcopy(start_intv)
            start_intv.end = end_intv.end
            logger.info(
                f"Extend amplicon interval {ai1} {prev_intv} to {start_intv}."
            )

            # These will all be collapsed into the starting interval
            intvs_to_remove = range(ai1 + 1, ai2 + 1)

            # Modify ccid
            for ai in intvs_to_remove:
                intermediate_intv = self.amplicon_intervals[ai]
                if intermediate_intv.amplicon_id != start_intv.amplicon_id:
                    ccid = intermediate_intv.amplicon_id
                    logger.info(f"Reset amplicon intervals with ccid {ccid}.")
                    for ai_, interval in enumerate(self.amplicon_intervals):
                        if interval.amplicon_id != ccid:
                            continue
                        logger.info(
                            f"Reset ccid of amplicon interval {ai_} ("
                            f"{self.amplicon_intervals[ai_]}) to "
                            f"{start_intv.amplicon_id}."
                        )
                        interval.amplicon_id = start_intv.amplicon_id
                        logger.info(f"Updated amplicon interval: {interval}")

            # Create duplicate connections, to modify while merging
            connection_map = {
                connection: connection
                for connection in self.amplicon_interval_connections
            }

            merged_idx = original_interval_idxs[ai1]
            for ai in intvs_to_remove:
                for prev_conn in connection_map:
                    prev_src, prev_dst = connection_map[prev_conn]

                    # Update connections from prev interval to merged interval,
                    # in ascending order
                    if ai == prev_src:
                        connection_map[prev_conn] = (
                            min(merged_idx, prev_dst),
                            max(merged_idx, prev_dst),
                        )
                    if ai == prev_dst:
                        prev_src = connection_map[prev_conn][0]
                        connection_map[prev_conn] = (
                            min(prev_src, merged_idx),
                            max(prev_src, merged_idx),
                        )

            for prev_conn, new_conn in connection_map.items():
                if prev_conn != new_conn:
                    logger.info(
                        f"Reset connection between amplicon intervals "
                        f"{prev_conn} to {new_conn}, i.e., "
                        f"{orig_intervals[prev_conn[0]]} and "
                        f"{orig_intervals[prev_conn[1]]}, to "
                        f"{original_interval_idxs.index(new_conn[0])} and "
                        f"{original_interval_idxs.index(new_conn[1])}."
                    )
                    self.amplicon_interval_connections[new_conn] |= (
                        self.amplicon_interval_connections[prev_conn]
                    )
                    del self.amplicon_interval_connections[prev_conn]

                    if new_conn[0] == new_conn[1]:
                        del self.amplicon_interval_connections[new_conn]

            # Delete intervals
            for ai in reversed(intvs_to_remove):
                intv = self.amplicon_intervals[ai]
                logger.info(f"Delete amplicon interval {ai} - {intv}.")
                del self.amplicon_intervals[ai]
                del original_interval_idxs[ai]

        for ai, interval in enumerate(self.amplicon_intervals):
            logger.info(
                f"\t{ai} (Prev. {original_interval_idxs[ai]}) - "
                f"Amplicon interval: {interval}"
            )
        # Need to remap connections based on reduced set of intervals
        orig_to_new_idx = {
            old_idx: new_idx
            for new_idx, old_idx in enumerate(original_interval_idxs)
        }
        logger.info(f"{self.amplicon_interval_connections=}")
        self.amplicon_interval_connections = defaultdict(
            set,
            {
                (orig_to_new_idx[conn[0]], orig_to_new_idx[conn[1]]): bp_idxs
                for conn, bp_idxs in self.amplicon_interval_connections.items()
            },
        )

        logger.info("Resetting amplicon ccids.")
        self.reset_amplicon_ccids()

    def reset_amplicon_ccids(self) -> None:
        ai_explored = np.zeros(len(self.amplicon_intervals))
        for ai, interval in enumerate(self.amplicon_intervals):
            ai_ccid = interval.amplicon_id
            if ai_explored[ai] == 0:
                ai_index_queue: list[int] = [ai]  # BFS queue
                while len(ai_index_queue) > 0:
                    ai_ = ai_index_queue.pop(0)
                    ai_explored[ai_] = 1
                    curr_intv = self.amplicon_intervals[ai_]
                    if curr_intv.amplicon_id != ai_ccid:
                        logger.info(
                            f"Reset amplicon interval of {curr_intv} to "
                            f"{ai_ccid}."
                        )
                        curr_intv.amplicon_id = ai_ccid
                    for ai1, ai2 in self.amplicon_interval_connections:
                        if ai1 == ai_ and ai_explored[ai2] == 0:
                            ai_index_queue.append(ai2)
                        elif ai2 == ai_ and ai_explored[ai1] == 0:
                            ai_index_queue.append(ai1)

    def sort_amplicon_intervals_and_connections(self) -> None:
        """Sort amplicon intervals in ascending order and remap connections
        accordingly."""
        sorted_intv_idxs = sorted(
            range(len(self.amplicon_intervals)),
            key=lambda x: self.amplicon_intervals[x],
        )
        self.amplicon_intervals = [
            self.amplicon_intervals[idx] for idx in sorted_intv_idxs
        ]

        # Remap connections using sorted indices in ascending interval idx order
        old_to_sorted_idx = {
            old_idx: sorted_idx
            for sorted_idx, old_idx in enumerate(sorted_intv_idxs)
        }
        self.amplicon_interval_connections = defaultdict(
            set,
            {
                (
                    min(
                        old_to_sorted_idx[old_conn[0]],
                        old_to_sorted_idx[old_conn[1]],
                    ),
                    max(
                        old_to_sorted_idx[old_conn[0]],
                        old_to_sorted_idx[old_conn[1]],
                    ),
                ): bp_idxs
                for old_conn, bp_idxs in self.amplicon_interval_connections.items()
            },
        )

    def find_amplicon_intervals(self) -> None:
        # Reset seed intervals
        logger.debug("Widening seed amplicon intervals based on CN segments.")
        self.widen_seed_intervals()
        ccid = 0
        for ai, interval in enumerate(self.amplicon_intervals):
            if interval.amplicon_id == -1:
                logger.debug(f"Begin processing amplicon interval {ai}")
                logger.debug(f"\tAmplicon interval {interval}")
                self.find_interval_i(ai, ccid)
                ccid += 1
        logger.info(
            f"Identified {len(self.amplicon_intervals)} amplicon intervals in "
            "total.",
        )
        for i, interval in enumerate(self.amplicon_intervals):
            logger.info(f"\t{i} - Amplicon interval: {interval}")

        for i, (ai, a2) in enumerate(self.amplicon_interval_connections):
            logger.info(f"\t{i} - Connection: {ai} -> {a2}")

        self.sort_amplicon_intervals_and_connections()

        # Merge amplicon intervals
        logger.info("Begin merging adjacent intervals.")
        self.merge_amplicon_intervals()

        logger.info(
            f"There are {len(self.amplicon_intervals)} amplicon intervals after"
            " merging."
        )
        for i, interval in enumerate(self.amplicon_intervals):
            logger.info(f"\t{i} - Amplicon interval: {interval}")

    def addbp(
        self,
        bp_: Breakpoint,
        bpr_: set[BPAlignments],
        bp_stats_: BreakpointStats,
        ccid: int,
    ) -> int:
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            if bp.is_close(bp_):
                logger.info(f"Updated existing breakpoint at {bpi}: {bp}.")
                logger.debug(f"Added {len(bpr_)} reads: {bpr_}.")
                self.new_bp_list[bpi].all_alignments |= set(bpr_)
                return bpi

        bpi = len(self.new_bp_list)
        logger.info(f"Added new breakpoint at {bpi}: {bp_}.")
        self.new_bp_list.append(bp_)
        self.new_bp_ccids.append(ccid)
        self.new_bp_stats.append(bp_stats_)
        return bpi

    def add_breakpoints_and_connections(
        self, new_bp_clusters: list[list[Breakpoint]], ccid: int | None = None
    ) -> list[int]:
        new_refined_bpi = []
        for initial_cluster in new_bp_clusters:
            if len(initial_cluster) < self.min_cluster_cutoff:
                continue
            bp_cluster_r = initial_cluster
            num_subcluster = 0
            while len(bp_cluster_r) >= self.min_cluster_cutoff:
                logger.info(
                    f"\t\t\tSubcluster {num_subcluster} of size "
                    f"{len(bp_cluster_r)}"
                )
                bp, bpr, bp_stats, bp_cluster_r = bpc2bp(
                    bp_cluster_r,
                    self.min_bp_match_cutoff_,
                )
                bp.all_alignments |= set(bpr)
                read_support = len(set(bpr))

                logger.info(f"\t\t\t\t{bp=}")
                logger.info(f"\t\t\t\tNum long read support = {read_support}")
                logger.info(f"\t\t\t\tbp_stats = {(bp_stats)}")

                if (read_support < self.min_cluster_cutoff) and (
                    read_support
                    < max(
                        self.normal_cov * self.min_bp_cov_factor,
                        3.0,
                    )
                ):
                    logger.info("\t\t\tDiscarded the subcluster.")
                    num_subcluster += 1
                    continue

                # Create breakpoint during initial interval refinement
                if ccid is not None:
                    logger.info("\t\t\tKept the subcluster.")
                    bpi = self.addbp(bp, set(bpr), bp_stats, ccid)
                    if bpi not in new_refined_bpi:
                        new_refined_bpi.append(bpi)
                elif num_subcluster == 0:
                    bp_start_intv = Interval(bp.chr1, bp.pos1, bp.pos1)
                    bp_end_intv = Interval(bp.chr2, bp.pos2, bp.pos2)
                    io1 = interval_overlap_l(
                        bp_start_intv, self.amplicon_intervals
                    )
                    io2 = interval_overlap_l(
                        bp_end_intv, self.amplicon_intervals
                    )
                    if io1 is None or io2 is None:
                        continue
                    amp_intv1 = self.amplicon_intervals[io1]
                    amp_intv2 = self.amplicon_intervals[io2]
                    if amp_intv1.amplicon_id != amp_intv2.amplicon_id:
                        raise ValueError(
                            f"Breakpoint endpoints matched to intervals "
                            f"{io1} and {io2} from different amplicons "
                            f"{amp_intv1.amplicon_id} and "
                            f"{amp_intv2.amplicon_id}."
                        )
                    bpi = self.addbp(
                        bp,
                        set(bpr),
                        bp_stats,
                        amp_intv1.amplicon_id,
                    )
                    self.amplicon_interval_connections[
                        (min(io1, io2), max(io1, io2))
                    ].add(bpi)
                    logger.info(
                        f"Added breakpoint connection {bpi} between "
                        f"amplicon intervals {io1} and {io2}."
                    )
                num_subcluster += 1

        return new_refined_bpi

    def get_refined_amplicon_interval_same_chr(
        self, chr_tag: core_types.ChrTag, left: BPToCNI, right: BPToCNI
    ) -> AmpliconInterval:
        """Get refined amplicon interval for two intervals on the same chromosome."""
        cns_intvs = self.cns_intervals_by_chr[chr_tag]
        cns_tree = self.cns_tree[chr_tag]

        left_cns, right_cns = cns_intvs[left.cni], cns_intvs[right.cni]
        first_cns, last_cns = cns_intvs[0], cns_intvs[-1]

        is_left_amplified = left_cns.cn >= self.cn_gain
        is_right_amplified = right_cns.cn >= self.cn_gain

        if not is_left_amplified:
            l = max(left.pos - self.interval_delta, first_cns.start)
        else:
            l = max(left_cns.start - self.interval_delta, first_cns.start)
        if not is_right_amplified:
            r = min(right.pos + self.interval_delta, last_cns.end)
        else:
            r = min(right_cns.end + self.interval_delta, last_cns.end)
        if left_cns.cn and left.pos - int(self.max_seq_len / 2) > l:
            l = left.pos - int(self.max_seq_len / 2)
        r = min(right.pos + int(self.max_seq_len / 2), r)
        if cns_tree.get_single_cns_idx(l) is None:
            l = left_cns.start
        if cns_tree.get_single_cns_idx(r) is None:
            r = right_cns.end

        return AmpliconInterval(chr_tag, l, r, -1)

    def get_refined_amplicon_intervals_same_chr(
        self, chr_: core_types.ChrTag, nint_segs: Sequence[BPToCNI]
    ) -> tuple[list[AmpliconInterval], list[list[core_types.BPIdx]]]:
        new_intervals_refined: list[AmpliconInterval] = []
        new_intervals_connections: list[list[core_types.BPIdx]] = []
        prev_i = 0
        for i in range(len(nint_segs) - 1):
            cns_intvs = self.cns_intervals_by_chr[chr_]

            curr_cni = nint_segs[i][0]
            curr_cns = cns_intvs[curr_cni]
            curr_pos = nint_segs[i][1]

            next_cni = nint_segs[i + 1][0]
            next_cns = cns_intvs[next_cni]
            next_pos = nint_segs[i + 1][1]

            cns_gap = next_cns.start - curr_cns.end

            is_amplified = (next_cns.cn >= self.cn_gain) or (
                curr_cns.cn >= self.cn_gain
            )
            if (
                (next_cni - curr_cni > 2)
                or (cns_gap > self.max_seq_len / 2)
                or (next_pos - curr_pos > self.max_seq_len)
                or (not is_amplified and cns_gap > 2 * self.interval_delta)
                or (
                    not is_amplified
                    and next_pos - curr_pos > 3 * self.interval_delta
                )
            ):
                refined_intv = self.get_refined_amplicon_interval_same_chr(
                    chr_, nint_segs[prev_i], nint_segs[i]
                )
                new_intervals_refined.append(refined_intv)

                new_intervals_connections.append(
                    [nint_segs[i_].bp_idx for i_ in range(prev_i, i + 1)]
                )
                logger.debug(f"\t\t\tFixed new interval: {refined_intv}.")
                logger.debug(
                    "\t\t\tList of breakpoints connected to the new interval:",
                )
                for bpi in new_intervals_connections[-1]:
                    logger.debug(f"\t\t\t\t{self.new_bp_list[bpi]}")

                prev_i = i + 1

        if len(nint_segs) > 0:
            refined_intv = self.get_refined_amplicon_interval_same_chr(
                chr_, nint_segs[prev_i], nint_segs[-1]
            )

            new_intervals_refined.append(refined_intv)
            logger.info(f"\t\t\tFixed new interval: {refined_intv}.")

            new_intervals_connections.append(
                [nint_segs[i_].bp_idx for i_ in range(prev_i, len(nint_segs))]
            )
            logger.info(
                "\t\t\tList of breakpoints connected to the new interval:",
            )
            for bpi in new_intervals_connections[-1]:
                logger.info(f"\t\t\t\t{self.new_bp_list[bpi]}")

        return new_intervals_refined, new_intervals_connections

    def get_refined_amplicon_interval(
        self, left: BPToChrCNI, right: BPToChrCNI
    ) -> AmpliconInterval:
        """Get refined amplicon interval for two intervals on potentially
        different chromosomes."""
        # TODO: unify this with get_refined_amplicon_interval_same_chr
        left_cns = self.cns_intervals_by_chr[left.chr][left.cni]
        right_cns = self.cns_intervals_by_chr[right.chr][right.cni]

        # First cns uses left chr, last CNS uses right chr
        # since intervals are pre-sorted with chr as first sort key
        first_cns = self.cns_intervals_by_chr[left.chr][0]
        last_cns = self.cns_intervals_by_chr[right.chr][-1]

        is_left_amplified = left_cns.cn >= self.cn_gain
        is_right_amplified = right_cns.cn >= self.cn_gain

        if not is_left_amplified:
            l = max(left.pos - self.interval_delta, first_cns.start)
        else:
            l = max(left_cns.start - self.interval_delta, first_cns.start)
        if not is_right_amplified:
            r = min(right.pos + self.interval_delta, last_cns.end)
        else:
            r = min(right_cns.end + self.interval_delta, last_cns.end)
        if left_cns.cn and left.pos - int(self.max_seq_len / 2) > l:
            l = left.pos - int(self.max_seq_len / 2)
        r = min(right.pos + int(self.max_seq_len / 2), r)
        if self.cns_tree[left.chr].get_single_cns_idx(l) is None:
            l = left_cns.start
        if self.cns_tree[right.chr].get_single_cns_idx(r) is None:
            r = right_cns.end

        return AmpliconInterval(left.chr, l, r, -1)

    def get_refined_amplicon_intervals(
        self, nint_segs_: Sequence[BPToChrCNI]
    ) -> tuple[list[AmpliconInterval], list[list[core_types.BPIdx]]]:
        new_intervals_refined: list[AmpliconInterval] = []
        new_intervals_connections: list[list[core_types.BPIdx]] = []
        prev_i = 0
        for i in range(len(nint_segs_) - 1):
            curr_seg = nint_segs_[i]
            next_seg = nint_segs_[i + 1]

            # two intervals in nint_segs_ might be on different chrs
            next_cns = self.cns_intervals_by_chr[next_seg.chr][next_seg.cni]
            curr_cns = self.cns_intervals_by_chr[curr_seg.chr][curr_seg.cni]

            nil = next_cns.start
            ncn = next_cns.cn
            lir = curr_cns.end
            lcn = curr_cns.cn
            is_amplified = (ncn >= self.cn_gain) or (lcn >= self.cn_gain)
            if (
                (nint_segs_[i + 1][0] != nint_segs_[i][0])
                or (nint_segs_[i + 1][1] - nint_segs_[i][1] > 2)
                or (nil - lir > self.max_seq_len / 2)
                or (nint_segs_[i + 1][2] - nint_segs_[i][2] > self.max_seq_len)
                or (not is_amplified and nil - lir > 2 * self.interval_delta)
                or (
                    not is_amplified
                    and nint_segs_[i + 1][2] - nint_segs_[i][2]
                    > 3 * self.interval_delta
                )
            ):
                new_intv = self.get_refined_amplicon_interval(
                    nint_segs_[prev_i], nint_segs_[i]
                )
                new_intervals_refined.append(new_intv)
                logger.debug(f"\t\tFixed new interval: {new_intv}.")

                new_intervals_connections.append([])
                logger.debug(
                    "\t\tSkip breakpoints connected to the new interval."
                )
                prev_i = i + 1
        if len(nint_segs_) > 0:
            left, right = nint_segs_[prev_i], nint_segs_[-1]
            new_intv = self.get_refined_amplicon_interval(left, right)
            new_intervals_refined.append(new_intv)
            logger.debug(f"\t\tFixed new interval: {new_intv}.")
            new_intervals_connections.append([])
            logger.debug("\t\tSkip breakpoints connected to the new interval:")
        return new_intervals_refined, new_intervals_connections

    def refine_initial_interval(
        self,
        interval: Interval,
        ccid: int,
        chr_: core_types.ChrTag,
        cni_intv: Interval,
        reads: set[core_types.ReadName],
    ) -> tuple[list[AmpliconInterval], list[list[core_types.BPIdx]]]:
        chr_tag, start_cns, end_cns = (
            cni_intv.chr,
            cni_intv.start,
            cni_intv.end,
        )
        ns = self.cns_intervals_by_chr[chr_tag][start_cns].start
        ne = self.cns_intervals_by_chr[chr_tag][end_cns].end
        cns_intv = Interval(chr_, ns, ne)
        logger.info(
            f"\t\tRefining new interval {[chr_, ns, ne]}, # reads = {len(reads)}."
        )

        new_bp_list = []

        for read_name in reads:
            new_bp_list += alignment2bp(
                read_name,
                self.chimeric_alignments[read_name],
                self.min_bp_match_cutoff_,
                20,
                cns_intv,
                interval,
                max_nm=(
                    self.nm_stats[0] + 3 * self.nm_stats[1]
                    if self.nm_filter
                    else None
                ),
            )

        logger.info(
            f"\t\tFound {len(new_bp_list)} reads connecting the two intervals."
        )

        new_bp_clusters = cluster_bp_list(
            new_bp_list,
            self.min_cluster_cutoff,
            self.max_breakpoint_distance_cutoff,
        )
        logger.info(f"\t\tThese reads formed {len(new_bp_clusters)} clusters.")
        new_bp_refined = self.add_breakpoints_and_connections(
            new_bp_clusters, ccid=ccid
        )

        new_intervals_refined: list[AmpliconInterval] = []
        new_intervals_connections: list[list[core_types.BPIdx]] = []
        if len(new_bp_refined) > 0:
            nint_segs, nint_segs_ = self.get_refined_cni_intvs(
                interval, cni_intv, new_bp_refined
            )
            new_intvs, new_conns = self.get_refined_amplicon_intervals_same_chr(
                chr_, nint_segs
            )
            new_intervals_refined.extend(new_intvs)
            new_intervals_connections.extend(new_conns)
            new_intvs, new_conns = self.get_refined_amplicon_intervals(
                nint_segs_
            )
            new_intervals_refined.extend(new_intvs)
            new_intervals_connections.extend(new_conns)
            logger.debug(f"\t\t\tCreated {new_intervals_connections=}")
        return new_intervals_refined, new_intervals_connections

    def get_refined_cni_intvs(
        self,
        interval: Interval,
        cni_intv: Interval,
        new_bpi_refined: list[int],
    ) -> tuple[Sequence[BPToCNI], Sequence[BPToChrCNI]]:
        same_chr_segs: list[BPToCNI] = []
        diff_chr_segs: list[BPToChrCNI] = []

        cni_chr = cni_intv.chr
        ns = self.cns_intervals_by_chr[cni_chr][cni_intv.start].start
        ne = self.cns_intervals_by_chr[cni_chr][cni_intv.end].end
        cns_intv = Interval(cni_chr, ns, ne)

        for bpi in new_bpi_refined:
            bp = self.new_bp_list[bpi]
            try:
                if interval.contains(bp.chr1, bp.pos1) and cns_intv.contains(
                    bp.chr2, bp.pos2
                ):
                    cns_tree = self.cns_tree[bp.chr2]
                    if (cni := cns_tree.get_single_cns_idx(bp.pos2)) is None:
                        raise ValueError(
                            f"Unable to find CNS index for {bp.pos2}."
                        )
                    same_chr_segs.append(BPToCNI(cni, bp.pos2, bpi))
                elif interval.contains(bp.chr2, bp.pos2) and cns_intv.contains(
                    bp.chr1, bp.pos1
                ):
                    cns_tree = self.cns_tree[bp.chr1]
                    if (cni := cns_tree.get_single_cns_idx(bp.pos1)) is None:
                        raise ValueError(
                            f"Unable to find CNS index for {bp.pos1}."
                        )
                    same_chr_segs.append(BPToCNI(cni, bp.pos1, bpi))
                else:
                    logger.warning(
                        "\t\tExact breakpoint outside amplicon interval."
                    )
                    logger.warning(f"\t\tBreakpoint {bp}.")
                    logger.warning(f"\t\tCurrent interval {interval}.")
                    logger.warning(f"\t\tNew interval {cns_intv}.")

                    chr1_tree = self.cns_tree[bp.chr1]
                    chr2_tree = self.cns_tree[bp.chr2]
                    if (
                        start_cni := chr1_tree.get_single_cns_idx(bp.pos1)
                    ) is None:
                        raise ValueError(
                            f"Unable to find CNS index for {bp.pos1}."
                        )
                    if (
                        end_cni := chr2_tree.get_single_cns_idx(bp.pos2)
                    ) is None:
                        raise ValueError(
                            f"Unable to find CNS index for {bp.pos2}."
                        )

                    # TODO: fix, separate set of segments generated here
                    if cni_intv.contains(bp.chr1, bp.pos1):
                        same_chr_segs.append(BPToCNI(start_cni, bp.pos1, bpi))
                    else:
                        diff_chr_segs.append(
                            BPToChrCNI(bp.chr1, start_cni, bp.pos1, bpi)
                        )
                    if cni_intv.contains(bp.chr2, bp.pos2):
                        same_chr_segs.append(BPToCNI(end_cni, bp.pos2, bpi))
                    else:
                        diff_chr_segs.append(
                            BPToChrCNI(bp.chr2, end_cni, bp.pos2, bpi)
                        )
            except:
                logger.exception(f"Unable to parse {bp=}")

        return sorted(same_chr_segs), sorted(
            diff_chr_segs, key=lambda item: (CHR_TAG_TO_IDX[item[0]], item[1:])
        )

    def find_interval_i(self, ai: int, ccid: int) -> None:
        """Given an amplification interval I indexed by ai, search for
            amplification intervals connected with I iteratively (with BFS)
            by a breakpoint edge
        Assign I a connected component id ccid if not already assigned
        """

        logger.debug(f"Start BFS on amplicon interval {ai}.")
        interval_queue: list[int] = [ai]  # BFS queue
        while len(interval_queue) > 0:
            logger.info(f"\tBFS queue: {interval_queue}")
            ai_: int = interval_queue.pop(0)
            interval = self.amplicon_intervals[ai_]
            logger.info(f"\tNext amplicon interval {ai_}: {interval}.")
            ai_chr = interval.chr
            if interval.amplicon_id == -1:
                interval.amplicon_id = ccid
            logger.debug(f"\t\tReset connected component ID to {ccid}")

            # Identify all amplification intervals connected to interval indexed
            # by ai_ with a breakpoint edge
            try:
                si, ei = self.cns_tree[interval.chr].get_cns_ends(interval)
            except KeyError:
                logger.error(f"Unable to find CNS ends for {interval=}")
                continue

            # Chr -> CNS idxs sharing a chimerical alignment -> associated reads
            chr_to_cni_to_reads: defaultdict[
                core_types.ChrTag,
                defaultdict[core_types.CNSIdx, set[core_types.ReadName]],
            ] = defaultdict(lambda: defaultdict(set))

            cni_to_chimeric_ref_intvs = self.chr_to_cni_to_chimeric_ref_intvs[
                ai_chr
            ]
            for i in range(si, ei + 1):
                for ref_intv in cni_to_chimeric_ref_intvs[i]:
                    read_name = ref_intv.name
                    for alignment in self.chimeric_alignments[read_name]:
                        chimera_chr = alignment.ref_interval.chr
                        for cns_idx in alignment.cns_idxs:
                            if chimera_chr == ai_chr and si < cns_idx < ei:
                                continue
                            chr_to_cni_to_reads[chimera_chr][cns_idx].add(
                                read_name
                            )
            chr_to_cni_to_reads = (
                breakpoint_utilities.filter_small_breakpoint_clusters(
                    chr_to_cni_to_reads, self.min_cluster_cutoff
                )
            )

            new_intervals_refined: list[AmpliconInterval] = []
            new_intervals_connections: list[list[core_types.BPIdx]] = []

            for chr_, cni_to_reads in sorted(
                chr_to_cni_to_reads.items(),
                key=lambda item: CHR_TAG_TO_IDX[item[0]],
            ):
                # Verify if there are any CN segments on this chr
                if len(cni_to_reads) == 0:
                    continue
                logger.debug(f"\t\tFound new intervals on chr {chr_=}")
                sorted_cns_idxs = sorted(cni_to_reads)
                # Initial list of new amplicon intervals, based on CN seg idxs
                cni_intv_to_reads: defaultdict[
                    tuple[
                        core_types.ChrTag, core_types.CNSIdx, core_types.CNSIdx
                    ],
                    set[str],
                ] = defaultdict(set)
                matching_reads = set()
                cns_intvs = self.cns_intervals_by_chr[chr_]

                prev_cni = sorted_cns_idxs[0]
                # Iterate through CN seg (index) pairs, step size 1
                for curr_cni, next_cni in zip(
                    sorted_cns_idxs, sorted_cns_idxs[1:]
                ):
                    nil = cns_intvs[next_cni].start
                    lir = cns_intvs[curr_cni].end
                    if next_cni - curr_cni > 2 or nil - lir > self.max_seq_len:
                        matching_reads |= chr_to_cni_to_reads[chr_][curr_cni]
                        cni_intv_to_reads[(chr_, prev_cni, curr_cni)] = (
                            matching_reads
                        )
                        prev_cni = next_cni
                        matching_reads = set()
                    else:
                        matching_reads |= chr_to_cni_to_reads[chr_][curr_cni]
                matching_reads |= chr_to_cni_to_reads[chr_][sorted_cns_idxs[-1]]
                cni_intv_to_reads[(chr_, prev_cni, sorted_cns_idxs[-1])] = (
                    matching_reads
                )

                # Refine initial intervals
                for cni_intv, matching_reads in cni_intv_to_reads.items():
                    new_intvs, new_conns = self.refine_initial_interval(
                        interval,
                        ccid,
                        chr_,
                        # Convert tuple (used for dict key) to Interval
                        Interval(*cni_intv),
                        matching_reads,
                    )
                    new_intervals_refined.extend(new_intvs)
                    new_intervals_connections.extend(new_conns)

            # BFS
            logger.debug("\t\tProcessing new intervals.")
            for ni in range(len(new_intervals_refined)):
                eis, intl = interval_exclusive(
                    new_intervals_refined[ni], self.amplicon_intervals
                )
                if len(intl) == 0:
                    ei_str = ""
                    for ei_ in eis:
                        ei_str += f"{self.amplicon_intervals[ei_]} "
                    ei_str = ei_str.rstrip()
                    logger.debug(
                        f"\t\tNew interval {new_intervals_refined[ni]} "
                        f"overlaps with existing interval {ei_str}."
                    )
                    for bpi in new_intervals_connections[ni]:
                        bp = self.new_bp_list[bpi]
                        for ei_ in eis:
                            connection = (min(ai_, ei_), max(ai_, ei_))
                            e_intv = self.amplicon_intervals[ei_]
                            if (
                                ei_ != ai_
                                and e_intv.contains(bp.chr1, bp.pos1)
                                or e_intv.contains(bp.chr2, bp.pos2)
                            ):
                                self.amplicon_interval_connections[
                                    connection
                                ].add(bpi)
                    for ei_ in eis:
                        e_intv = self.amplicon_intervals[ei_]
                        if ei_ != ai_ and e_intv.amplicon_id < 0:
                            interval_queue.append(ei_)
                else:
                    for int_ in intl:
                        nai = len(self.amplicon_intervals)
                        self.amplicon_intervals.append(int_)
                        logger.debug(
                            f"\t\tAdded new interval {int_} to the amplicon "
                            "interval list."
                        )
                        logger.debug(f"\t\tNew interval index: {nai}")
                        self.amplicon_interval_connections[(ai_, nai)] = set()
                        if len(eis) == 0:
                            for bpi in new_intervals_connections[ni]:
                                self.amplicon_interval_connections[
                                    (ai_, nai)
                                ].add(bpi)
                        else:
                            for bpi in new_intervals_connections[ni]:
                                bp = self.new_bp_list[bpi]
                                for ei_ in eis:
                                    e_intv = self.amplicon_intervals[ei_]
                                    connection = (min(ai_, ei_), max(ai_, ei_))
                                    if e_intv.contains(
                                        bp.chr1, bp.pos1
                                    ) or e_intv.contains(bp.chr2, bp.pos2):
                                        self.amplicon_interval_connections[
                                            connection
                                        ].add(bpi)
                                    else:
                                        self.amplicon_interval_connections[
                                            (ai_, nai)
                                        ].add(bpi)
                        interval_queue.append(nai)

    def find_breakpoints(self) -> None:
        """Search for breakpoints from chimeric alignments within identified amplified intervals
        Then cluster the breakpoints from chimeric alignments
        """
        new_bp_list_: list[Breakpoint] = []
        for r, chimeras in self.chimeric_alignments.items():
            new_bp_list_ += alignment2bp_l(
                r,
                chimeras,
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

        self.add_breakpoints_and_connections(new_bp_clusters)

    def find_smalldel_breakpoints(self) -> None:
        """Search for breakpoints from a single alignment record, within
        identified amplified intervals
        For each alignment record with large indels, first try to match the
        resulting breakpoint to the list of AA breakpoints
        Cluster the unmatched breakpoints
        """
        new_bp_list_: list[Breakpoint] = []
        nm_threshold = (
            self.nm_stats[0] + 3 * self.nm_stats[1] if self.nm_filter else None
        )
        for ai in self.amplicon_intervals:
            for read in self.bam.fetch_interval(ai):
                indel_alignments = (
                    breakpoint_utilities.get_indel_alignments_from_read(
                        ai.chr, read, self.min_del_len, nm_threshold
                    )
                )
                if indel_alignments:
                    self.large_indel_alignments[read.query_name].extend(
                        indel_alignments
                    )
        logger.info(
            f"Fetched {len(self.large_indel_alignments)} reads with large indels in CIGAR."
        )

        for r in self.large_indel_alignments:
            for i, alignment in enumerate(self.large_indel_alignments[r]):
                new_bp_list_.append(
                    Breakpoint(
                        node1=Node(
                            alignment.chr_tag,
                            max(alignment.next_start, alignment.curr_end),
                            Strand.REVERSE,
                        ),
                        node2=Node(
                            alignment.chr_tag,
                            min(alignment.next_start, alignment.curr_end),
                            Strand.FORWARD,
                        ),
                        alignment_info=BPAlignments(r, i, i),
                        gap=0,
                        was_reversed=True,
                        mapq1=-1,
                        mapq2=-1,
                    )
                )
        logger.info(
            f"Found {len(new_bp_list_)} reads with new small del breakpoints."
        )

        new_bp_clusters = cluster_bp_list(
            new_bp_list_,
            self.min_cluster_cutoff,
            self.max_breakpoint_distance_cutoff,
        )
        logger.info(f"These reads formed {len(new_bp_clusters)} clusters.")
        self.add_breakpoints_and_connections(new_bp_clusters)

    def build_graphs(self) -> None:
        """Organize the identified discordant edges into a list of breakpoint
        graphs, stored in lr_graph.
        Each graph represent a connected component of amplified intervals,
        i.e., amplicons.
        """
        # Split amplified intervals according to discordant edges
        split_intvs_by_ai: dict[int, list[SplitInterval]] = defaultdict(list)
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            for ai, interval in enumerate(self.amplicon_intervals):
                seg = interval
                if bp.chr1 == seg.chr and seg.start < bp.pos1 < seg.end:
                    if bp.strand1 == Strand.FORWARD:
                        split_intvs_by_ai[ai].append(
                            SplitInterval(bp.pos1, bp.pos1 + 1, Strand.FORWARD)
                        )
                    else:
                        split_intvs_by_ai[ai].append(
                            SplitInterval(bp.pos1 - 1, bp.pos1, Strand.REVERSE)
                        )
                if bp.chr2 == seg.chr and seg.start < bp.pos2 < seg.end:
                    if bp.strand2 == Strand.FORWARD:
                        split_intvs_by_ai[ai].append(
                            SplitInterval(bp.pos2, bp.pos2 + 1, Strand.FORWARD)
                        )
                    else:
                        split_intvs_by_ai[ai].append(
                            SplitInterval(bp.pos2 - 1, bp.pos2, Strand.REVERSE)
                        )
        # Construct graphs with sequence and concordant edges
        logger.info(
            f"Will split the following {len(split_intvs_by_ai)} amplicon "
            "intervals into sequence edges and build breakpoint graphs."
        )
        amplicon_id = 1
        for sseg in self.amplicon_intervals:
            if sseg.amplicon_id not in self.ccid2id:
                self.ccid2id[sseg.amplicon_id] = amplicon_id
                amplicon_id += 1

        for cci in range(len(self.ccid2id)):
            # Initialize breakpoint graph objects
            self.lr_graph.append(BreakpointGraph())

        for ai, split_intvs in split_intvs_by_ai.items():
            logger.info(
                f"Will split the amplicon interval at index {ai} "
                f"{self.amplicon_intervals[ai]}."
            )
            for split_intv in split_intvs:
                logger.info(f"\tSplit interval at {split_intv}.")
        for ai, split_intvs in split_intvs_by_ai.items():
            split_intvs.sort(key=lambda item: item.start)
            sseg = self.amplicon_intervals[ai]
            amplicon_idx = self.ccid2id[sseg.amplicon_id] - 1
            bp_graph = self.lr_graph[amplicon_idx]
            for ssi in range(len(split_intvs_by_ai[ai])):
                curr_intv = split_intvs_by_ai[ai][ssi]
                if ssi == 0:
                    bp_graph.add_sequence_edge(
                        sseg.chr,
                        sseg.start,
                        curr_intv.start,
                    )
                    bp_graph.add_concordant_edge(
                        Node(sseg.chr, curr_intv.start, Strand.FORWARD),
                        Node(sseg.chr, curr_intv.end, Strand.REVERSE),
                    )
                elif (
                    curr_intv.start
                    > (prev_intv := split_intvs_by_ai[ai][ssi - 1]).start
                ):
                    bp_graph.add_sequence_edge(
                        sseg.chr,
                        prev_intv.end,
                        curr_intv.start,
                    )
                    bp_graph.add_concordant_edge(
                        Node(sseg.chr, curr_intv.start, Strand.FORWARD),
                        Node(sseg.chr, curr_intv.end, Strand.REVERSE),
                    )
            last_intv = split_intvs_by_ai[ai][-1]
            bp_graph.add_sequence_edge(sseg.chr, last_intv.end, sseg.end)
        for ai in range(len(self.amplicon_intervals)):
            if ai not in split_intvs_by_ai:
                sseg = self.amplicon_intervals[ai]
                amplicon_idx = self.ccid2id[sseg.amplicon_id] - 1
                self.lr_graph[amplicon_idx].add_sequence_edge(
                    sseg.chr, sseg.start, sseg.end
                )
        for amplicon_idx in range(len(self.lr_graph)):
            self.lr_graph[amplicon_idx].sort_edges()

        # Add nodes corresponding to interval ends
        for ai in range(len(self.amplicon_intervals)):
            sseg = self.amplicon_intervals[ai]
            amplicon_idx = self.ccid2id[sseg.amplicon_id] - 1
            self.lr_graph[amplicon_idx].amplicon_intervals.append(
                AmpliconInterval(sseg.chr, sseg.start, sseg.end, amplicon_idx)
            )
            self.lr_graph[amplicon_idx].add_endnode(
                Node(sseg.chr, sseg.start, Strand.REVERSE)
            )
            self.lr_graph[amplicon_idx].add_endnode(
                Node(sseg.chr, sseg.end, Strand.FORWARD)
            )
        logger.debug(
            "The following nodes correspond to interval ends.",
        )
        for amplicon_idx in range(len(self.lr_graph)):
            for node in self.lr_graph[amplicon_idx].endnode_adjacencies:
                logger.debug(f"\tAmplicon {amplicon_idx + 1}, node {node}.")

        # Construct graphs with discordant and source edges
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            bp_ccid = self.new_bp_ccids[bpi]
            io1 = interval_overlap_l(
                Interval(bp.chr1, bp.pos1, bp.pos1), self.amplicon_intervals
            )
            io2 = interval_overlap_l(
                Interval(bp.chr2, bp.pos2, bp.pos2), self.amplicon_intervals
            )
            assert io1 is not None and io2 is not None
            intv1 = self.amplicon_intervals[io1]
            intv2 = self.amplicon_intervals[io2]
            assert intv1.amplicon_id == intv2.amplicon_id

            amplicon_idx = self.ccid2id[intv1.amplicon_id] - 1
            if intv1.amplicon_id != bp_ccid:
                logger.debug(
                    f"Reset the ccid for breakpoint {bp} at index {bpi} from "
                    f"{bp_ccid} to {intv1.amplicon_id}."
                )
                self.new_bp_ccids[bpi] = intv1.amplicon_id
            self.lr_graph[amplicon_idx].add_discordant_edge(
                bp.node1,
                bp.node2,
                read_support=len(bp.all_alignments),
                alignments=bp.all_alignments,
            )

        # Print summary statistics for each amplicon
        for amplicon_idx, bp_graph in enumerate(self.lr_graph):
            idx = amplicon_idx + 1
            logger.debug(
                f"Num sequence edges in amplicon {idx} = {bp_graph.num_seq_edges}."
            )
            logger.debug(
                f"Num concordant edges in amplicon {idx} = {bp_graph.num_conc_edges}."
            )
            logger.debug(
                f"Num discordant edges in amplicon {idx} = {bp_graph.num_disc_edges}."
            )
            logger.debug(
                f"Num source edges in amplicon {idx} = {bp_graph.num_src_edges}."
            )

    def assign_cov(self) -> None:
        """Extract the long read coverage from bam file, if missing, for
        each sequence edge."""
        for amplicon_idx in range(len(self.lr_graph)):
            for seq_edge in self.lr_graph[amplicon_idx].sequence_edges:
                if seq_edge.lr_count == -1:
                    logger.debug(
                        f"Finding LR cov for sequence edge {seq_edge}."
                    )
                    """
					For long read, use the total number of nucleotides
					"""
                    rl_list = [
                        read
                        for read in self.lr_bamfh.fetch(
                            seq_edge.chr, seq_edge.start, seq_edge.end + 1
                        )
                        if read.infer_read_length()
                    ]
                    seq_edge.lr_count = len(rl_list)
                    seq_edge.lr_nc = self.bam.count_raw_coverage(
                        seq_edge.interval
                    )
                    logger.debug(
                        f"LR cov assigned for sequence edge {seq_edge}."
                    )
        """
		Extract the long read coverage from bam file, if missing, for each
        concordant edge
		"""
        for amplicon_idx in range(len(self.lr_graph)):
            bp_graph = self.lr_graph[amplicon_idx]
            for eci in range(len(bp_graph.concordant_edges)):
                ec = bp_graph.concordant_edges[eci]
                logger.debug(f"Finding cov for concordant edge {ec}.")
                rls = {
                    read.query_name for read in self.bam.fetch_node(ec.node1)
                }
                rrs = {
                    read.query_name for read in self.bam.fetch_node(ec.node2)
                }
                rls1 = {
                    read.query_name
                    for read in self.lr_bamfh.fetch(
                        contig=ec.node1.chr,
                        start=ec.node1.pos - self.min_bp_match_cutoff_ - 1,
                        stop=ec.node1.pos - self.min_bp_match_cutoff_,
                    )
                }
                rrs1 = {
                    read.query_name
                    for read in self.lr_bamfh.fetch(
                        contig=ec.node2.chr,
                        start=ec.node2.pos + self.min_bp_match_cutoff_,
                        stop=ec.node2.pos + self.min_bp_match_cutoff_ + 1,
                    )
                }
                rbps = set()
                for bpi in (
                    self.lr_graph[amplicon_idx]
                    .node_adjacencies[ec.node1]
                    .discordant
                ):
                    for r in bp_graph.discordant_edges[bpi].alignments:
                        rbps.add(r[0])
                for bpi in bp_graph.node_adjacencies[ec.node2].discordant:
                    for r in bp_graph.discordant_edges[bpi].alignments:
                        rbps.add(r[0])
                bp_graph.concordant_edges[eci].read_names = rls | rrs
                bp_graph.concordant_edges[eci].lr_count = len(
                    (rls & rrs & rls1 & rrs1) - rbps,
                )
                logger.debug(f"LR cov assigned for concordant edge {ec}.")

    def compute_bp_graph_path_constraints(
        self, bp_graph: BreakpointGraph
    ) -> list[PathConstraint]:
        bp_path_constraints: list[PathConstraint] = []
        read_to_alignments: defaultdict[str, BPIndexedAlignmentContainer] = (
            defaultdict(BPIndexedAlignmentContainer)
        )
        concordant_reads = {}

        for di, discordant_edge in enumerate(bp_graph.discordant_edges):
            for bp_alignments in discordant_edge.alignments:
                if bp_alignments.alignment1 == bp_alignments.alignment2:
                    read_to_alignments[bp_alignments.name].equal.append(
                        BPIndexedAlignments(
                            alignment1=bp_alignments.alignment1,
                            alignment2=bp_alignments.alignment2,
                            discordant_idx=di,
                        )
                    )
                else:
                    read_to_alignments[bp_alignments.name].unequal.append(
                        BPIndexedAlignments(
                            alignment1=bp_alignments.alignment1,
                            alignment2=bp_alignments.alignment2,
                            discordant_idx=di,
                        )
                    )
        logger.debug(
            f"There are {len(read_to_alignments)} reads covering >=1 breakpoint"
            f" in amplicon."
        )

        for rn, disc_alignments in read_to_alignments.items():
            paths = path_utilities.get_bp_graph_paths(
                bp_graph,
                rn,
                disc_alignments,
                self.chimeric_alignments.get(rn, []),
                self.large_indel_alignments.get(rn, []),
                self.min_bp_match_cutoff_,
            )

            for path in paths:
                if len(path) <= 5 or not path_constraints.valid_path(
                    bp_graph, path
                ):
                    continue
                existing_paths = [pc.path for pc in bp_path_constraints]
                if path in existing_paths:
                    pci = existing_paths.index(path)
                    bp_path_constraints[pci].support += 1
                elif path[::-1] in existing_paths:
                    pci = existing_paths.index(path[::-1])
                    bp_path_constraints[pci].support += 1
                else:
                    bp_path_constraints.append(
                        PathConstraint(
                            path=path,
                            support=1,
                            amplicon_id=bp_graph.amplicon_idx,
                        )
                    )
        logger.debug(
            f"There are {len(bp_path_constraints)} distinct "
            f"subpaths due to reads involving breakpoints in amplicon "
            f"{bp_graph.amplicon_idx + 1}."
        )

        # Extract reads in concordant_edges_reads
        lc = len(bp_graph.concordant_edges)
        for ci in range(lc):
            for rn in bp_graph.concordant_edges[ci].read_names:
                if (
                    rn not in self.large_indel_alignments
                    and rn not in self.chimeric_alignments
                ):
                    concordant_reads[rn] = bp_graph.amplicon_idx + 1
        logger.debug(
            f"There are {len(concordant_reads)} concordant reads within "
            f"amplicon intervals in amplicon {bp_graph.amplicon_idx + 1}."
        )
        for aint in bp_graph.amplicon_intervals:
            for read in self.lr_bamfh.fetch(aint.chr, aint.start, aint.end + 1):
                rn = read.query_name  # type: ignore[arg-type, assignment]
                if read.mapping_quality >= 20 and rn in concordant_reads:
                    path = path_constraints.alignment_to_path(
                        bp_graph,
                        Interval(
                            read.reference_name,  # type: ignore[arg-type]
                            read.reference_start,
                            read.reference_end,  # type: ignore[arg-type]
                        ),
                    )
                    if len(path) <= 5 or not path_constraints.valid_path(
                        bp_graph, path
                    ):
                        continue
                    existing_paths = [pc.path for pc in bp_path_constraints]
                    if path in existing_paths:
                        pci = existing_paths.index(path)
                        bp_path_constraints[pci].support += 1
                    elif path[::-1] in existing_paths:
                        pci = existing_paths.index(path[::-1])
                        bp_path_constraints[pci].support += 1
                    else:
                        logger.info(
                            f"Adding path {path} to path constraints, {aint=}"
                        )
                        bp_path_constraints.append(
                            PathConstraint(
                                path=path,
                                support=1,
                                amplicon_id=bp_graph.amplicon_idx,
                            )
                        )

        logger.debug(
            f"There are {len(bp_path_constraints)} distinct "
            f"subpaths in amplicon {bp_graph.amplicon_idx + 1}."
        )
        return bp_path_constraints

    def compute_path_constraints(self) -> None:
        """Convert reads mapped within the amplicons into subpath constraints"""
        for amplicon_idx, bp_graph in enumerate(self.lr_graph):
            bp_graph.amplicon_idx = amplicon_idx
            bp_path_constraints = self.compute_bp_graph_path_constraints(
                bp_graph
            )
            bp_graph.path_constraints = bp_path_constraints
            bp_graph.longest_path_constraints = (
                path_constraints.longest_path_dict(bp_path_constraints)
            )

    def closebam(self):
        """Close the short read and long read bam file"""
        self.lr_bamfh.close()


@core_utils.profile_fn
def reconstruct_graphs(
    lr_bam_filename: pathlib.Path,
    cnv_seed_file: typer.FileText,
    cn_seg_file: typer.FileText,
    output_prefix: str,
    output_path_constraints: OutputPCOptions,
    min_bp_support: float,
) -> LongReadBamToBreakpointMetadata:
    seed_intervals = breakpoint_utilities.get_intervals_from_seed_file(
        cnv_seed_file
    )  # type: ignore[arg-type] # file is passed as IO stream
    logger.info(f"Parsed {len(seed_intervals)} seed amplicon intervals.")
    for interval in seed_intervals:
        logger.info(f"Seed interval: {interval}")

    b2bn = LongReadBamToBreakpointMetadata(
        lr_bamfh=pysam.AlignmentFile(str(lr_bam_filename), "rb"),
        bam=bam_types.BAMWrapper(str(lr_bam_filename), "rb"),
        amplicon_intervals=seed_intervals,
    )
    # if filter_bp_by_edit_distance:
    # b2bn.nm_filter = True
    b2bn.min_bp_cov_factor = min_bp_support
    logger.info("Opened LR bam files.")

    b2bn.read_cns(cn_seg_file)
    logger.info("Completed parsing CN segment files.")

    # Write chimeric alignments to pickle file to avoid re-fetching from BAM
    # on successive reconstructions
    pickle_path = pathlib.Path(f"{output_prefix}_chimeric_alignments.pickle")
    if pickle_path.exists():
        try:
            with pickle_path.open("rb") as file:
                chimeric_alignments = pickle.load(file)
        except Exception as e:
            logger.error(
                f"Unable to load chimeric alignments: {e}, re-fetching"
            )
            chimeric_alignments, edit_dist_stats = (
                breakpoint_utilities.fetch_breakpoint_reads(b2bn.bam)
            )
            with pickle_path.open("wb") as file:
                pickle.dump(chimeric_alignments, file)
    else:
        start = time.time()
        chimeric_alignments, edit_dist_stats = (
            breakpoint_utilities.fetch_breakpoint_reads(b2bn.bam)
        )
        with pickle_path.open("wb") as file:
            pickle.dump(chimeric_alignments, file)
        logger.error(f"Time to fetch breakpoint reads: {time.time() - start}")
    logger.info("Completed fetching reads containing breakpoints.")

    chimeric_alignments_by_read, chr_cns_to_chimeric_alignments = (
        b2bn.hash_alignment_to_seg(chimeric_alignments)
    )
    b2bn.hash_alignment_to_seg(chimeric_alignments_by_read)

    b2bn.find_amplicon_intervals()
    logger.info("Completed finding amplicon intervals.")
    b2bn.find_smalldel_breakpoints()
    logger.info("Completed finding small del breakpoints.")
    b2bn.find_breakpoints()
    logger.info("Completed finding all discordant breakpoints.")
    b2bn.build_graphs()
    logger.info("Breakpoint graphs built for all amplicons.")
    
    b2bn.assign_cov()
    logger.info(
        "Fetched read coverage for all sequence and concordant edges."
    )
    for gi in range(len(b2bn.lr_graph)):
        b2bn.lr_graph[gi].compute_cn_lr(b2bn.normal_cov)
    logger.info("Computed CN for all edges.")
    b2bn.compute_path_constraints()
    logger.info("Computed all subpath constraints.")

    for gi in range(len(b2bn.lr_graph)):
        breakpoint_utilities.output_breakpoint_graph_lr(
            b2bn.lr_graph[gi], f"{output_prefix}_amplicon{gi+1}_graph.txt", output_path_constraints
        )
        #file_prefix = f"{output_prefix}_amplicon{gi+1}"
    logger.info(
        f"Wrote breakpoint graph for all complicons to {output_prefix}_amplicon*_graph.txt."
    )

    return b2bn
