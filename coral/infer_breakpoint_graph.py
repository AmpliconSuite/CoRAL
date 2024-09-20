from __future__ import annotations

import functools
import logging
import math
import sys
import time
from dataclasses import dataclass, field
from typing import Any, List, Optional

import intervaltree  # type: ignore[import-untyped]
import numpy as np
import pysam

from coral import breakpoint_graph, cigar_parsing, state_provider
from coral.breakpoint_graph import BreakpointGraph
from coral.constants import CHR_TAG_TO_IDX
from coral.types import AmpliconInterval, CnsInterval

edge_type_to_index = {"s": 0, "c": 1, "d": 2}


from typing import Dict, List, Optional

import pysam


@dataclass
class BamToBreakpointNanopore:
    lr_bamfh: Optional[pysam.AlignmentFile] = None  # Long read bam file

    lr_graph: list[BreakpointGraph] = field(default_factory=list)  # Breakpoint graph

    max_seq_len: int = 2000000
    cn_gain: float = 5.0
    min_bp_cov_factor: float = 1.0
    min_bp_match_cutoff_: int = 100
    interval_delta: int = 100000
    min_cluster_cutoff: int = 3  # Hard cutoff for considering a long read breakpoint cluster
    max_breakpoint_distance_cutoff: int = 2000  # Used for breakpoint clustering - if the distance of two breakpoint positions are greater than this cutoff, then start a new cluster
    min_del_len: int = 600  # The minimum length of all +- (deletion) breakpoints returned by AA
    # small_del_cutoff = 10000 # +- breakpoints (small deletions) with the two ends less than this cutoff are treated specially

    read_length: dict[str, int] = field(default_factory=dict)  # Map read name -> read length

    # Map read name -> chimeric alignments (two or more records for one read)
    chimeric_alignments: dict[str, list[pysam.AlignedSegment]] = field(default_factory=dict)
    chimeric_alignments_seg: dict[str, list[pysam.AlignedSegment]] = field(default_factory=dict)

    # Map read name -> alignments with one record per read but large indels showing in CIGAR string
    large_indel_alignments: dict[str, list[pysam.AlignedSegment]] = field(default_factory=dict)
    # For edit distance filter of breakpoints
    nm_stats: list[float] = field(default_factory=functools.partial(lambda: [0.0] * 3))  # Default of [0.0, 0.0, 0.0]
    nm_filter: bool = False

    amplicon_intervals: list[AmpliconInterval] = field(default_factory=list)  # AA amplicon intervals
    amplicon_interval_connections: dict[str, List[AmpliconInterval]] = field(default_factory=dict)

    cns_intervals: list[CnsInterval] = field(default_factory=list)  # CN segments
    cns_intervals_by_chr: dict[str, list[CnsInterval]] = field(default_factory=dict)
    log2_cn: list[float] = field(default_factory=list)  # log2 CN for each CN segment
    cns_tree: dict[str, intervaltree.IntervalTree] = field(
        default_factory=dict
    )  # Interval tree structure for each chromosome
    normal_cov: float = 0.0  # Normal long read coverage - used for CN assignment

    ccid2id: dict[int, int] = field(default_factory=dict)
    new_bp_list: list[list[int]] = field(default_factory=list)  # List of breakpoints (discordant edges)
    new_bp_stats: list[list[int]] = field(default_factory=list)  # Statistics of breakpoints (discordant edges)
    new_bp_ccids: list[int] = field(default_factory=list)
    source_edges: list[list[Any]] = field(default_factory=list)
    source_edge_ccids: list[int] = field(default_factory=list)
    path_constraints: dict[int, list[list[Any]]] = field(default_factory=dict)
    longest_path_constraints: dict[int, list[list[Any]]] = field(default_factory=dict)
    cycles: dict[int, list[list[list[int]]]] = field(default_factory=dict)  # cycles, paths
    cycle_weights: dict[int, list[list[float]]] = field(default_factory=dict)  # cycles, paths
    path_constraints_satisfied: dict[int, list[list[list[int]]]] = field(default_factory=dict)  # cycles, paths

    def __init__(self, lr_bamfile, seedfile):
        self.lr_bamfh = pysam.AlignmentFile(lr_bamfile, "rb")
        with open(seedfile) as fp:
            for line in fp:
                s = line.strip().split()
                self.amplicon_intervals.append((s[0], int(s[1]), int(s[2]), -1))
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Parsed %d seed amplicon intervals." % (len(self.amplicon_intervals)),
        )
        for ai in self.amplicon_intervals:
            logging.debug(
                "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Seed interval %s" % (ai[:3]),
            )

    def read_cns(self, cns):
        """Read in (cnvkit) *.cns file and estimate the normal long read coverage"""
        self.cns_intervals = []
        self.log2_cn = []
        idx = 0
        with open(cns) as fp:
            for line in fp:
                s = line.strip().split()
                if s[0] != "chromosome":
                    self.cns_intervals.append((s[0], int(s[1]), int(s[2]) - 1))
                    if s[0] not in self.cns_tree:
                        self.cns_tree[s[0]] = intervaltree.IntervalTree()
                        self.cns_intervals_by_chr[s[0]] = []
                        idx = 0
                    self.cns_tree[s[0]][int(s[1]) : int(s[2])] = idx
                    idx += 1
                    if cns.endswith(".cns"):
                        self.cns_intervals_by_chr[s[0]].append(
                            [s[0], int(s[1]), int(s[2]) - 1, 2 * (2 ** float(s[4]))],
                        )
                        self.log2_cn.append(float(s[4]))
                    elif cns.endswith(".bed"):
                        self.cns_intervals_by_chr[s[0]].append(
                            [s[0], int(s[1]), int(s[2]) - 1, float(s[3])],
                        )
                        self.log2_cn.append(np.log2(float(s[3]) / 2.0))
                    else:
                        sys.stderr.write(cns + "\n")
                        sys.stderr.write("Invalid cn_seg file format!\n")

        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Total num LR copy number segments: %d." % (len(self.log2_cn)),
        )
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
        total_int_len += self.cns_intervals[log2_cn_order[ip]][2] - self.cns_intervals[log2_cn_order[ip]][1] + 1
        total_int_len += self.cns_intervals[log2_cn_order[im]][2] - self.cns_intervals[log2_cn_order[im]][1] + 1
        i = 1
        while total_int_len < 10000000:
            cns_intervals_median.append(self.cns_intervals[log2_cn_order[ip + i]])
            cns_intervals_median.append(self.cns_intervals[log2_cn_order[im - i]])
            log2_cn_median.append(self.log2_cn[log2_cn_order[ip]])
            log2_cn_median.append(self.log2_cn[log2_cn_order[im]])
            total_int_len += (
                self.cns_intervals[log2_cn_order[ip + i]][2] - self.cns_intervals[log2_cn_order[ip + i]][1] + 1
            )
            total_int_len += (
                self.cns_intervals[log2_cn_order[im - i]][2] - self.cns_intervals[log2_cn_order[im - i]][1] + 1
            )
            i += 1
        logging.debug(state_provider.get_timed_task_log(f"Use {len(cns_intervals_median)} LR copy number segments."))
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Total length of LR copy number segments: %d." % (total_int_len),
        )
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Average LR copy number: %f." % (np.average(log2_cn_median)),
        )
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
                ],
            )
        self.normal_cov = nnc * 1.0 / total_int_len
        logging.info(
            "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "LR normal cov = %f." % (self.normal_cov),
        )
        self.min_cluster_cutoff = max(
            self.min_cluster_cutoff,
            self.min_bp_cov_factor * self.normal_cov,
        )
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Reset min_cluster_cutoff to %f." % (self.min_cluster_cutoff),
        )

    def fetch(self):
        for read in self.lr_bamfh.fetch():
            rn = read.query_name
            if read.flag < 256 and rn not in self.read_length:
                self.read_length[rn] = read.query_length
            try:
                sa_list = read.get_tag("SA:Z:")[:-1].split(";")
                for sa in sa_list:
                    try:
                        if sa not in self.chimeric_alignments[rn]:
                            self.chimeric_alignments[rn].append(sa)
                    except:
                        self.chimeric_alignments[rn] = [sa]
            except:
                if read.mapping_quality == 60:
                    e = read.get_cigar_stats()[0][-1] / read.query_length
                    self.nm_stats[0] += e
                    self.nm_stats[1] += e * e
                    self.nm_stats[2] += 1
        self.nm_stats[0] /= self.nm_stats[2]
        self.nm_stats[1] = math.sqrt(self.nm_stats[1] / self.nm_stats[2] - self.nm_stats[0] ** 2)
        logging.info(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Fetched %d chimeric reads." % (len(self.chimeric_alignments)),
        )
        reads_wo_primary_alignment = []
        for r in self.chimeric_alignments.keys():
            if r not in self.read_length:
                logging.warning(
                    f"#TIME"
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "Found chimeric read without primary alignment.",
                )
                logging.warning(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "\tRead name: %s; Read length: N/A." % r,
                )
                logging.warning(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "\tAll CIGAR strings: %s." % (self.chimeric_alignments[r]),
                )
                reads_wo_primary_alignment.append(r)
                continue
            rl = self.read_length[r]
            self.chimeric_alignments[r] = cigar_parsing.alignment_from_satags(
                self.chimeric_alignments[r],
                rl,
            )
        for r in reads_wo_primary_alignment:
            del self.chimeric_alignments[r]
        logging.info(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Computed alignment intervals on all chimeric reads.",
        )

    def pos2cni(self, chr, pos):
        return self.cns_tree[chr][pos]

    def hash_alignment_to_seg(self):
        """Speed up amplified interval search by hashing chimeric alignments from each long read to CN segments"""
        for read in self.chimeric_alignments.keys():
            for ri in range(len(self.chimeric_alignments[read][1])):
                rint = self.chimeric_alignments[read][1][ri]
                if rint[0] in self.cns_tree:
                    lcni_ = self.pos2cni(rint[0], min(rint[1], rint[2]))
                    rcni_ = self.pos2cni(rint[0], max(rint[1], rint[2]))
                    assert len(lcni_) <= 1 and len(rcni_) <= 1
                    lcni, rcni = -1, -1
                    if len(lcni_) == 1:
                        lcni = list(lcni_)[0].data
                    if len(rcni_) == 1:
                        rcni = list(rcni_)[0].data
                    cniset = set([lcni, rcni])
                    if len(cniset) > 1 and -1 in cniset:
                        cniset.remove(-1)
                    self.chimeric_alignments[read][1][ri].append(cniset)
                    if rint[0] not in self.chimeric_alignments_seg:
                        self.chimeric_alignments_seg[rint[0]] = dict()
                    for cni in cniset:
                        if cni != -1:
                            try:
                                self.chimeric_alignments_seg[rint[0]][cni].append(read)
                            except:
                                self.chimeric_alignments_seg[rint[0]][cni] = [read]
                else:
                    self.chimeric_alignments[read][1][ri].append(set([-1]))

    def find_amplicon_intervals(self):
        # Reset seed intervals
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Updating seed amplicon intervals according to CN segments.",
        )
        for ai in range(len(self.amplicon_intervals)):
            chr = self.amplicon_intervals[ai][0]
            lcni = list(self.pos2cni(chr, self.amplicon_intervals[ai][1]))[0].data
            rcni = list(self.pos2cni(chr, self.amplicon_intervals[ai][2]))[0].data
            self.amplicon_intervals[ai][1] = self.cns_intervals_by_chr[chr][lcni][1]
            if (
                len(
                    list(
                        self.pos2cni(
                            chr,
                            self.cns_intervals_by_chr[chr][lcni][1] - self.interval_delta,
                        ),
                    ),
                )
                > 0
            ):
                self.amplicon_intervals[ai][1] = self.cns_intervals_by_chr[chr][lcni][1] - self.interval_delta
            self.amplicon_intervals[ai][2] = self.cns_intervals_by_chr[chr][rcni][2]
            if (
                len(
                    list(
                        self.pos2cni(
                            chr,
                            self.cns_intervals_by_chr[chr][rcni][2] + self.interval_delta,
                        ),
                    ),
                )
                > 0
            ):
                self.amplicon_intervals[ai][2] = self.cns_intervals_by_chr[chr][rcni][2] + self.interval_delta
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "\tUpdated amplicon interval %s" % self.amplicon_intervals[ai],
            )

        ccid = 0
        for ai in range(len(self.amplicon_intervals)):
            if self.amplicon_intervals[ai][3] == -1:
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "Begin processing amplicon interval %d" % ai,
                )
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "\tAmplicon interval %s" % self.amplicon_intervals[ai],
                )
                self.find_interval_i(ai, ccid)
                ccid += 1
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Identified %d amplicon intervals in total." % len(self.amplicon_intervals),
        )
        sorted_ai_indices = sorted(
            range(len(self.amplicon_intervals)),
            key=lambda i: (
                CHR_TAG_TO_IDX[self.amplicon_intervals[i][0]],
                self.amplicon_intervals[i][1],
            ),
        )
        amplicon_intervals_sorted = [self.amplicon_intervals[i] for i in sorted_ai_indices]
        for ai in range(len(amplicon_intervals_sorted)):
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "\tAmplicon interval %s" % amplicon_intervals_sorted[ai],
            )

        # Merge amplicon intervals
        logging.debug(
            "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Begin merging adjacent intervals.",
        )
        lastai = 0
        intervals_to_merge = []
        for ai in range(len(amplicon_intervals_sorted) - 1):
            if not (
                interval_adjacent(amplicon_intervals_sorted[ai + 1], amplicon_intervals_sorted[ai])
                or interval_overlap(
                    amplicon_intervals_sorted[ai],
                    amplicon_intervals_sorted[ai + 1],
                )
            ):
                if ai > lastai:
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Merging intervals from %d to %d." % (lastai, ai),
                    )
                    intervals_to_merge.append([lastai, ai])
                lastai = ai + 1
        if len(amplicon_intervals_sorted) > 0 and lastai < len(amplicon_intervals_sorted) - 1:
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "Merging intervals from %d to %d." % (lastai, len(amplicon_intervals_sorted) - 1),
            )
            intervals_to_merge.append([lastai, len(amplicon_intervals_sorted) - 1])
        for int_range in intervals_to_merge[::-1]:
            # Reset interval
            amplicon_intervals_sorted[int_range[0]][2] = amplicon_intervals_sorted[int_range[1]][2]
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "Reset amplicon interval %d to %s." % (int_range[0], amplicon_intervals_sorted[int_range[0]]),
            )
            # Modify ccid
            for ai in range(int_range[0] + 1, int_range[1] + 1):
                if amplicon_intervals_sorted[ai][3] != amplicon_intervals_sorted[int_range[0]][3]:
                    ccid = amplicon_intervals_sorted[ai][3]
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Reset amplicon intervals with ccid %d." % ccid,
                    )
                    for ai_ in range(len(amplicon_intervals_sorted)):
                        if amplicon_intervals_sorted[ai_][3] == ccid:
                            amplicon_intervals_sorted[ai_][3] = amplicon_intervals_sorted[int_range[0]][3]
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "Reset ccid of amplicon interval %d to %d."
                                % (ai_, amplicon_intervals_sorted[int_range[0]][3]),
                            )
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "Updated amplicon interval: %s" % amplicon_intervals_sorted[ai_],
                            )
            # Modify interval connections
            connection_map = dict()
            for connection in self.amplicon_interval_connections.keys():
                connection_map[connection] = connection
            for ai in range(int_range[0] + 1, int_range[1] + 1):
                ai_unsorted_ = sorted_ai_indices[int_range[0]]
                ai_unsorted = sorted_ai_indices[ai]
                for connection in connection_map:
                    if ai_unsorted == connection_map[connection][0]:
                        connection_map[connection] = (ai_unsorted_, connection_map[connection][1])
                    if ai_unsorted == connection_map[connection][1]:
                        connection_map[connection] = (connection_map[connection][0], ai_unsorted_)
                    if connection_map[connection][1] < connection_map[connection][0]:
                        connection_map[connection] = (
                            connection_map[connection][1],
                            connection_map[connection][0],
                        )
            for connection in connection_map:
                if connection != connection_map[connection]:
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Reset connection between amplicon intervals %s to %s."
                        % (connection, connection_map[connection]),
                    )
                    if connection_map[connection] not in self.amplicon_interval_connections:
                        self.amplicon_interval_connections[connection_map[connection]] = (
                            self.amplicon_interval_connections[connection]
                        )
                    else:
                        self.amplicon_interval_connections[connection_map[connection]] |= (
                            self.amplicon_interval_connections[connection]
                        )
                    del self.amplicon_interval_connections[connection]
                    if connection_map[connection][0] == connection_map[connection][1]:
                        del self.amplicon_interval_connections[connection_map[connection]]
            # Delete intervals
            for ai in range(int_range[0] + 1, int_range[1] + 1)[::-1]:
                ai_unsorted = sorted_ai_indices[ai]
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "Delete amplicon intervals %d - %s." % (ai_unsorted, self.amplicon_intervals[ai_unsorted]),
                )
                del amplicon_intervals_sorted[ai]
                del sorted_ai_indices[ai]

        self.amplicon_intervals = [amplicon_intervals_sorted[ai] for ai in range(len(amplicon_intervals_sorted))]
        ind_map = {sorted_ai_indices[i]: i for i in range(len(sorted_ai_indices))}
        connection_map = {
            connection: (
                min(ind_map[connection[0]], ind_map[connection[1]]),
                max(ind_map[connection[0]], ind_map[connection[1]]),
            )
            for connection in self.amplicon_interval_connections.keys()
        }
        self.amplicon_interval_connections = {
            connection_map[connection]: self.amplicon_interval_connections[connection]
            for connection in self.amplicon_interval_connections.keys()
        }
        # Reset ccids
        ai_explored = np.zeros(len(self.amplicon_intervals))
        for ai in range(len(self.amplicon_intervals)):
            ai_ccid = self.amplicon_intervals[ai][3]
            if ai_explored[ai] == 0:
                L = [ai]  # BFS queue
                while len(L) > 0:
                    ai_ = L.pop(0)
                    ai_explored[ai_] = 1
                    if self.amplicon_intervals[ai_][3] != ai_ccid:
                        self.amplicon_intervals[ai_][3] = ai_ccid
                    for ai1, ai2 in self.amplicon_interval_connections.keys():
                        if ai1 == ai_ and ai_explored[ai2] == 0:
                            L.append(ai2)
                        elif ai2 == ai_ and ai_explored[ai1] == 0:
                            L.append(ai1)

        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "There are %d amplicon intervals after merging." % len(self.amplicon_intervals),
        )
        for ai in range(len(self.amplicon_intervals)):
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "\tAmplicon interval %s after merging." % self.amplicon_intervals[ai],
            )

    def addbp(self, bp_, bpr_, bp_stats_, ccid):
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
                # if bp[0] == bp[3] and bp[2] == '-' and bp[5] == '+' and abs(bp[1] - bp[4]) < self.small_del_cutoff:
                # self.small_del_indices.append(bpi)
                return bpi
        bpi = len(self.new_bp_list)
        self.new_bp_list.append(bp_ + [bpr_])
        self.new_bp_ccids.append(ccid)
        self.new_bp_stats.append(bp_stats_)
        # if bp_[0] == bp_[3] and bp_[2] == '-' and bp_[5] == '+' and abs(bp_[1] - bp_[4]) < self.small_del_cutoff:
        # self.small_del_indices.append(bpi)
        return bpi

    def find_interval_i(self, ai, ccid):
        """Given an amplification interval I indexed by ai, search for amplification intervals connected with I iteratively (with BFS)
                by a breakpoint edge
        Assign I a connected component id ccid if not already assigned
        """
        logging.debug(
            "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tStart BFS on amplicon interval %d." % ai,
        )
        L = [ai]  # BFS queue
        while len(L) > 0:
            logging.debug(
                "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\t\tBFS queue: %s" % L,
            )
            ai_ = L.pop(0)
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "\t\tNext amplicon interval %d: %s." % (ai_, self.amplicon_intervals[ai_]),
            )
            chr = self.amplicon_intervals[ai_][0]
            s = self.amplicon_intervals[ai_][1]
            e = self.amplicon_intervals[ai_][2]
            if self.amplicon_intervals[ai_][3] == -1:
                self.amplicon_intervals[ai_][3] = ccid
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "\t\tReset connected component ID to %d" % ccid,
            )

            # Identify all amplification intervals connected to interval indexed by ai_ with a breakpoint edge
            try:
                si = list(self.pos2cni(chr, s))[0].data
                ei = list(self.pos2cni(chr, e))[0].data
            except:
                continue

            d1_segs = dict()  # All CN segments which share a chimeric alignment to the given interval
            for i in range(si, ei + 1):
                if i in self.chimeric_alignments_seg[chr]:
                    for r in self.chimeric_alignments_seg[chr][i]:
                        rint = self.chimeric_alignments[r][1]
                        for int_ in rint:
                            for i_ in int_[-1]:
                                if (int_[0] != chr) or (i_ <= si or i_ >= ei):
                                    try:
                                        if i_ != -1:
                                            d1_segs[int_[0]][i_].add(r)
                                    except:
                                        if int_[0] not in d1_segs:
                                            d1_segs[int_[0]] = dict()
                                        if i_ != -1:
                                            d1_segs[int_[0]][i_] = set([r])
            # Initial filtering of potential breakpoints of insufficient support
            chr_del_list = []
            for chr_ in d1_segs:
                segs = d1_segs[chr_].keys()
                del_list = []
                for segi in segs:
                    if len(d1_segs[chr_][segi]) < self.min_cluster_cutoff:
                        del_list.append(segi)
                for segi in del_list:
                    del d1_segs[chr_][segi]
                if len(d1_segs[chr_]) == 0:
                    chr_del_list.append(chr_)
            for chr_ in chr_del_list:
                del d1_segs[chr_]

            new_intervals_refined = []
            new_intervals_connections = []
            for chr_ in d1_segs:
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "\t\tFound new intervals on chr %s" % chr_,
                )
                new_intervals = []  # Initial list of new amplicon intervals
                sorted_bins_chr_ = sorted(d1_segs[chr_].keys())
                nir = set([])
                lasti = 0
                for i in range(len(sorted_bins_chr_) - 1):
                    nil = self.cns_intervals_by_chr[chr_][sorted_bins_chr_[i + 1]][1]
                    lir = self.cns_intervals_by_chr[chr_][sorted_bins_chr_[i]][2]
                    if sorted_bins_chr_[i + 1] - sorted_bins_chr_[i] > 2 or nil - lir > self.max_seq_len:
                        nir |= d1_segs[chr_][sorted_bins_chr_[i]]
                        new_intervals.append(
                            [chr_, sorted_bins_chr_[lasti], sorted_bins_chr_[i], nir],
                        )
                        lasti = i + 1
                        nir = set([])
                    else:
                        nir |= d1_segs[chr_][sorted_bins_chr_[i]]
                nir |= d1_segs[chr_][sorted_bins_chr_[-1]]
                new_intervals.append([chr_, sorted_bins_chr_[lasti], sorted_bins_chr_[-1], nir])

                # Refine initial intervals
                for nint_ in new_intervals:
                    ns = self.cns_intervals_by_chr[nint_[0]][nint_[1]][1]
                    ne = self.cns_intervals_by_chr[nint_[0]][nint_[2]][2]
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\t\tRefining new interval %s." % [chr_, ns, ne],
                    )
                    new_bp_list = []
                    if self.nm_filter:
                        for r in nint_[-1]:
                            new_bp_list += alignment2bp_nm(
                                r,
                                self.chimeric_alignments[r],
                                self.min_bp_match_cutoff_,
                                20,
                                self.nm_stats[0] + 3 * self.nm_stats[1],
                                [nint_[0], ns, ne],
                                self.amplicon_intervals[ai_],
                            )
                    else:
                        for r in nint_[-1]:
                            new_bp_list += alignment2bp(
                                r,
                                self.chimeric_alignments[r],
                                self.min_bp_match_cutoff_,
                                20,
                                [nint_[0], ns, ne],
                                self.amplicon_intervals[ai_],
                            )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\t\tFound %d reads connecting the two intervals." % len(new_bp_list),
                    )
                    new_bp_clusters = cluster_bp_list(
                        new_bp_list,
                        self.min_cluster_cutoff,
                        self.max_breakpoint_distance_cutoff,
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\t\tThese reads formed %d clusters." % (len(new_bp_clusters)),
                    )
                    new_bp_refined = []
                    for c in new_bp_clusters:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\t\t\tNew cluster of size %d." % (len(c)),
                        )
                        if len(c) >= self.min_cluster_cutoff:
                            num_subcluster = 0
                            bp_cluster_r = c
                            while len(bp_cluster_r) >= self.min_cluster_cutoff:
                                bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(
                                    bp_cluster_r,
                                    self.min_bp_match_cutoff_,
                                )
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\tSubcluster %d" % (num_subcluster),
                                )
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\t\t\tbp = %s" % (bp),
                                )
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\t\t\tNum long read support = %d" % (len(set(bpr))),
                                )
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\t\t\tbp_stats = %s" % (bp_stats_),
                                )
                                if (num_subcluster == 0 and len(set(bpr)) >= self.min_cluster_cutoff) or (
                                    len(set(bpr)) >= max(self.normal_cov * self.min_bp_cov_factor, 3.0)
                                ):
                                    bpi = self.addbp(bp, set(bpr), bp_stats_, ccid)
                                    if bpi not in new_bp_refined:
                                        new_bp_refined.append(bpi)
                                    logging.debug(
                                        "#TIME "
                                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                                        + "\t\t\tKeeped the cluster.",
                                    )
                                else:
                                    logging.debug(
                                        "#TIME "
                                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                                        + "\t\t\tDiscarded the cluster.",
                                    )
                        else:
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\t\t\tDiscarded the cluster.",
                            )

                    nint_segs = []
                    nint_segs_ = []
                    if len(new_bp_refined) > 0:
                        for bpi in new_bp_refined:
                            bp = self.new_bp_list[bpi][:6]
                            try:
                                if interval_overlap(
                                    [bp[0], bp[1], bp[1]],
                                    self.amplicon_intervals[ai_],
                                ) and interval_overlap([bp[3], bp[4], bp[4]], [nint_[0], ns, ne]):
                                    nint_segs.append(
                                        [list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi],
                                    )
                                elif interval_overlap(
                                    [bp[3], bp[4], bp[4]],
                                    self.amplicon_intervals[ai_],
                                ) and interval_overlap([bp[0], bp[1], bp[1]], [nint_[0], ns, ne]):
                                    nint_segs.append(
                                        [list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi],
                                    )
                                else:
                                    logging.warning(
                                        "#TIME "
                                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                                        + "\t\tExact breakpoint outside amplicon interval.",
                                    )
                                    logging.warning(
                                        "#TIME "
                                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                                        + "\t\tBreakpoint %s." % bp,
                                    )
                                    logging.warning(
                                        "#TIME "
                                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                                        + "\t\tCurrent interval %s." % self.amplicon_intervals[ai_],
                                    )
                                    logging.warning(
                                        "#TIME "
                                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                                        + "\t\tNew interval %s." % [nint_[0], ns, ne],
                                    )
                                    o1 = interval_overlap([bp[0], bp[1], bp[1]], [nint_[0], ns, ne])
                                    o2 = interval_overlap([bp[3], bp[4], bp[4]], [nint_[0], ns, ne])
                                    if o1 and o2:
                                        nint_segs.append(
                                            [list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi],
                                        )
                                        nint_segs.append(
                                            [list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi],
                                        )
                                    elif o1:
                                        nint_segs.append(
                                            [list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi],
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
                                            [list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi],
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
                        nint_segs = sorted(nint_segs, key=lambda item: (item[0], item[1]))
                        nint_segs_ = sorted(
                            nint_segs_,
                            key=lambda item: (CHR_TAG_TO_IDX[item[0]], item[1], item[2]),
                        )
                        lasti = 0
                        for i in range(len(nint_segs) - 1):
                            nil = self.cns_intervals_by_chr[chr_][nint_segs[i + 1][0]][1]
                            ncn = self.cns_intervals_by_chr[chr_][nint_segs[i + 1][0]][3]
                            lir = self.cns_intervals_by_chr[chr_][nint_segs[i][0]][2]
                            lcn = self.cns_intervals_by_chr[chr_][nint_segs[i][0]][3]
                            amp_flag = (ncn >= self.cn_gain) or (lcn >= self.cn_gain)
                            if (
                                (nint_segs[i + 1][0] - nint_segs[i][0] > 2)
                                or (nil - lir > self.max_seq_len / 2)
                                or (nint_segs[i + 1][1] - nint_segs[i][1] > self.max_seq_len)
                                or (not amp_flag and nil - lir > 2 * self.interval_delta)
                                or (not amp_flag and nint_segs[i + 1][1] - nint_segs[i][1] > 3 * self.interval_delta)
                            ):
                                amp_flag_l = self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][3] >= self.cn_gain
                                amp_flag_r = self.cns_intervals_by_chr[chr_][nint_segs[i][0]][3] >= self.cn_gain
                                if not amp_flag_l:
                                    l = max(
                                        nint_segs[lasti][1] - self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][0][1],
                                    )
                                else:
                                    l = max(
                                        self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1] - self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][0][1],
                                    )
                                if not amp_flag_r:
                                    r = min(
                                        nint_segs[i][1] + self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][-1][2],
                                    )
                                else:
                                    r = min(
                                        lir + self.interval_delta,
                                        self.cns_intervals_by_chr[chr_][-1][2],
                                    )
                                if (
                                    self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][3]
                                    and nint_segs[lasti][1] - int(self.max_seq_len / 2) > l
                                ):
                                    l = nint_segs[lasti][1] - int(self.max_seq_len / 2)
                                r = min(nint_segs[i][1] + int(self.max_seq_len / 2), r)
                                if len(list(self.pos2cni(chr_, l))) == 0:
                                    l = self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1]
                                if len(list(self.pos2cni(chr_, r))) == 0:
                                    r = lir
                                new_intervals_refined.append([chr_, l, r, -1])
                                new_intervals_connections.append([])
                                for i_ in range(lasti, i + 1):
                                    new_intervals_connections[-1].append(nint_segs[i_][2])
                                lasti = i + 1
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\tFixed new interval: %s." % new_intervals_refined[-1],
                                )
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\tList of breakpoints connected to the new interval:",
                                )
                                for bpi in new_intervals_connections[-1]:
                                    logging.debug(
                                        "#TIME "
                                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                                        + "\t\t\t%s" % self.new_bp_list[bpi][:6],
                                    )
                        if len(nint_segs) > 0:
                            amp_flag_l = self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][3] >= self.cn_gain
                            amp_flag_r = self.cns_intervals_by_chr[chr_][nint_segs[-1][0]][3] >= self.cn_gain
                            if not amp_flag_l:
                                l = max(
                                    nint_segs[lasti][1] - self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][0][1],
                                )
                            else:
                                l = max(
                                    self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1] - self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][0][1],
                                )
                            if not amp_flag_r:
                                r = min(
                                    nint_segs[-1][1] + self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][-1][2],
                                )
                            else:
                                r = min(
                                    self.cns_intervals_by_chr[chr_][nint_segs[-1][0]][2] + self.interval_delta,
                                    self.cns_intervals_by_chr[chr_][-1][2],
                                )
                            if nint_segs[lasti][1] - int(self.max_seq_len / 2) > l:
                                l = nint_segs[lasti][1] - int(self.max_seq_len / 2) > l
                            r = min(nint_segs[-1][1] + int(self.max_seq_len / 2), r)
                            if len(list(self.pos2cni(chr_, l))) == 0:
                                l = self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1]
                            if len(list(self.pos2cni(chr_, r))) == 0:
                                r = self.cns_intervals_by_chr[chr_][nint_segs[-1][0]][2]
                            new_intervals_refined.append([chr_, l, r, -1])
                            new_intervals_connections.append([])
                            for i_ in range(lasti, len(nint_segs)):
                                new_intervals_connections[-1].append(nint_segs[i_][2])
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\t\tFixed new interval: %s." % new_intervals_refined[-1],
                            )
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\t\tList of breakpoints connected to the new interval:",
                            )
                            for bpi in new_intervals_connections[-1]:
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\t\t%s" % self.new_bp_list[bpi][:6],
                                )
                        lasti = 0
                        for i in range(len(nint_segs_) - 1):
                            # two intervals in nint_segs_ might be on different chrs
                            nil = self.cns_intervals_by_chr[nint_segs_[i + 1][0]][nint_segs_[i + 1][1]][1]
                            ncn = self.cns_intervals_by_chr[nint_segs_[i + 1][0]][nint_segs_[i + 1][1]][3]
                            lir = self.cns_intervals_by_chr[nint_segs_[i][0]][nint_segs_[i][1]][2]
                            lcn = self.cns_intervals_by_chr[nint_segs_[i][0]][nint_segs_[i][1]][3]
                            amp_flag = (ncn >= self.cn_gain) or (lcn >= self.cn_gain)
                            if (
                                (nint_segs_[i + 1][0] != nint_segs_[i][0])
                                or (nint_segs_[i + 1][1] - nint_segs_[i][1] > 2)
                                or (nil - lir > self.max_seq_len / 2)
                                or (nint_segs_[i + 1][2] - nint_segs_[i][2] > self.max_seq_len)
                                or (not amp_flag and nil - lir > 2 * self.interval_delta)
                                or (not amp_flag and nint_segs_[i + 1][2] - nint_segs_[i][2] > 3 * self.interval_delta)
                            ):
                                amp_flag_l = (
                                    self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][3]
                                    >= self.cn_gain
                                )
                                amp_flag_r = (
                                    self.cns_intervals_by_chr[nint_segs_[i][0]][nint_segs_[i][1]][3] >= self.cn_gain
                                )
                                if not amp_flag_l:
                                    l = max(
                                        nint_segs_[lasti][2] - self.interval_delta,
                                        self.cns_intervals_by_chr[nint_segs_[lasti][0]][0][1],
                                    )
                                else:
                                    l = max(
                                        self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1]
                                        - self.interval_delta,
                                        self.cns_intervals_by_chr[nint_segs_[lasti][0]][0][1],
                                    )
                                if not amp_flag_r:
                                    r = min(
                                        nint_segs_[i][2] + self.interval_delta,
                                        self.cns_intervals_by_chr[nint_segs_[i][0]][-1][2],
                                    )
                                else:
                                    r = min(
                                        lir + self.interval_delta,
                                        self.cns_intervals_by_chr[nint_segs_[i][0]][-1][2],
                                    )
                                l = max(nint_segs_[lasti][2] - int(self.max_seq_len / 2), l)
                                r = min(nint_segs_[i][2] + int(self.max_seq_len / 2), r)
                                if len(list(self.pos2cni(nint_segs_[lasti][0], l))) == 0:
                                    l = self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1]
                                if len(list(self.pos2cni(nint_segs_[i][0], r))) == 0:
                                    r = lir
                                new_intervals_refined.append([nint_segs_[lasti][0], l, r, -1])
                                new_intervals_connections.append([])
                                lasti = i + 1
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\tFixed new interval: %s." % new_intervals_refined[-1],
                                )
                                logging.debug(
                                    "#TIME "
                                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                                    + "\t\tSkip breakpoints connected to the new interval.",
                                )
                        if len(nint_segs_) > 0:
                            amp_flag_l = (
                                self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][3] >= self.cn_gain
                            )
                            amp_flag_r = (
                                self.cns_intervals_by_chr[nint_segs_[-1][0]][nint_segs_[-1][1]][3] >= self.cn_gain
                            )
                            if not amp_flag_l:
                                l = max(
                                    nint_segs_[lasti][2] - self.interval_delta,
                                    self.cns_intervals_by_chr[nint_segs_[lasti][0]][0][1],
                                )
                            else:
                                l = max(
                                    self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1]
                                    - self.interval_delta,
                                    self.cns_intervals_by_chr[nint_segs_[lasti][0]][0][1],
                                )
                            if not amp_flag_r:
                                r = min(
                                    nint_segs_[-1][2] + self.interval_delta,
                                    self.cns_intervals_by_chr[nint_segs_[-1][0]][-1][2],
                                )
                            else:
                                r = min(
                                    self.cns_intervals_by_chr[nint_segs_[-1][0]][nint_segs_[-1][1]][2]
                                    + self.interval_delta,
                                    self.cns_intervals_by_chr[nint_segs_[-1][0]][-1][2],
                                )
                            l = max(nint_segs_[lasti][2] - int(self.max_seq_len / 2), l)
                            r = min(nint_segs_[-1][2] + int(self.max_seq_len / 2), r)
                            if len(list(self.pos2cni(nint_segs_[lasti][0], l))) == 0:
                                l = self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1]
                            if len(list(self.pos2cni(nint_segs_[lasti][0], r))) == 0:
                                r = self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[-1][1]][2]
                            new_intervals_refined.append([nint_segs_[lasti][0], l, r, -1])
                            new_intervals_connections.append([])
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\t\tFixed new interval: %s." % new_intervals_refined[-1],
                            )
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\t\tSkip breakpoints connected to the new interval:",
                            )

            # BFS
            logging.debug(
                "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\t\tProcessing new intervals.",
            )
            for ni in range(len(new_intervals_refined)):
                ei, intl = interval_exclusive(new_intervals_refined[ni], self.amplicon_intervals)
                if len(intl) == 0:
                    ei_str = ""
                    for ei_ in ei:
                        ei_str += "%s " % self.amplicon_intervals[ei_]
                    ei_str = ei_str.rstrip()
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\t\tNew interval %s overlaps with existing interval %s."
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
                                    self.amplicon_interval_connections[connection].add(bpi)
                                except:
                                    self.amplicon_interval_connections[connection] = set([bpi])
                    for ei_ in ei:
                        if ei_ != ai_ and self.amplicon_intervals[ei_][3] < 0:
                            L.append(ei_)
                else:
                    for int_ in intl:
                        nai = len(self.amplicon_intervals)
                        self.amplicon_intervals.append(int_)
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\t\tAdded new interval %s to the amplicon interval list." % int_,
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\t\tNew interval index: %d." % nai,
                        )
                        self.amplicon_interval_connections[(ai_, nai)] = set([])
                        if len(ei) == 0:
                            for bpi in new_intervals_connections[ni]:
                                self.amplicon_interval_connections[(ai_, nai)].add(bpi)
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
                                            self.amplicon_interval_connections[connection].add(bpi)
                                        except:
                                            self.amplicon_interval_connections[connection] = set(
                                                [bpi],
                                            )
                                    else:
                                        self.amplicon_interval_connections[(ai_, nai)].add(bpi)
                        L.append(nai)

    def find_breakpoints(self):
        """Search for breakpoints from chimeric alignments within identified amplified intervals
        Then cluster the breakpoints from chimeric alignments
        """
        new_bp_list_ = []
        if self.nm_filter:
            for r in self.chimeric_alignments.keys():
                new_bp_list_ += alignment2bp_nm_l(
                    r,
                    self.chimeric_alignments[r],
                    self.min_bp_match_cutoff_,
                    20,
                    self.nm_stats[0] + 3 * self.nm_stats[1],
                    100,
                    self.amplicon_intervals,
                )
        else:
            for r in self.chimeric_alignments.keys():
                new_bp_list_ += alignment2bp_l(
                    r,
                    self.chimeric_alignments[r],
                    self.min_bp_match_cutoff_,
                    20,
                    100,
                    self.amplicon_intervals,
                )
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Found %d reads with new breakpoints." % (len(new_bp_list_)),
        )

        new_bp_clusters = cluster_bp_list(
            new_bp_list_,
            self.min_cluster_cutoff,
            self.max_breakpoint_distance_cutoff,
        )
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "These reads formed %d clusters." % (len(new_bp_clusters)),
        )
        for c in new_bp_clusters:
            logging.debug(
                "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "New cluster of size %d." % (len(c)),
            )
            if len(c) >= self.min_cluster_cutoff:
                num_subcluster = 0
                bp_cluster_r = c
                while len(bp_cluster_r) >= self.min_cluster_cutoff:
                    bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(
                        bp_cluster_r,
                        self.min_bp_match_cutoff_,
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tSubcluster %d" % (num_subcluster),
                    )
                    logging.debug(
                        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tbp = %s" % (bp),
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tNum long read support = %d" % (len(set(bpr))),
                    )
                    logging.debug(
                        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tbp_stats = %s" % (bp_stats_),
                    )
                    if (num_subcluster == 0 and len(set(bpr)) >= self.min_cluster_cutoff) or (
                        len(set(bpr)) >= max(self.normal_cov * self.min_bp_cov_factor, 3.0)
                    ):
                        io1 = interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals)
                        io2 = interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals)
                        if io1 >= 0 and io2 >= 0:
                            assert self.amplicon_intervals[io1][3] == self.amplicon_intervals[io2][3]
                            bpi = self.addbp(
                                bp,
                                set(bpr),
                                bp_stats_,
                                self.amplicon_intervals[io1][3],
                            )
                            try:
                                self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))].add(bpi)
                            except:
                                self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))] = set([bpi])
                    else:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tDiscarded the subcluster %d." % (num_subcluster),
                        )
                    num_subcluster += 1
            else:
                logging.debug(
                    "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tDiscarded the cluster.",
                )

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
                        if abs(blocks[bi + 1][0] - blocks[bi][1]) > self.min_del_len:
                            agg_del += abs(blocks[bi + 1][0] - blocks[bi][1])
                    rnm = (rnm - agg_del) / read.query_length
                    if rnm < self.nm_stats[0] + 3 * self.nm_stats[1]:
                        for bi in range(len(blocks) - 1):
                            if abs(blocks[bi + 1][0] - blocks[bi][1]) > self.min_del_len:
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
                        if abs(blocks[bi + 1][0] - blocks[bi][1]) > self.min_del_len:
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
        logging.info(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Fetched %d reads with large indels in CIGAR." % (len(self.large_indel_alignments)),
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
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Found %d reads with new small del breakpoints." % (len(new_bp_list_)),
        )

        new_bp_clusters = cluster_bp_list(
            new_bp_list_,
            self.min_cluster_cutoff,
            self.max_breakpoint_distance_cutoff,
        )
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "These reads formed %d clusters." % (len(new_bp_clusters)),
        )
        for c in new_bp_clusters:
            logging.debug(
                "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "New cluster of size %d." % (len(c)),
            )
            if len(c) >= self.min_cluster_cutoff:
                num_subcluster = 0
                bp_cluster_r = c
                while len(bp_cluster_r) >= self.min_cluster_cutoff:
                    bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(
                        bp_cluster_r,
                        self.min_bp_match_cutoff_,
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tSubcluster %d" % (num_subcluster),
                    )
                    logging.debug(
                        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tbp = %s" % (bp),
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tNum long read support = %d" % (len(set(bpr))),
                    )
                    logging.debug(
                        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tbp_stats = %s" % (bp_stats_),
                    )
                    if (num_subcluster == 0 and len(set(bpr)) >= self.min_cluster_cutoff) or (
                        len(set(bpr)) >= max(self.normal_cov * self.min_bp_cov_factor, 3.0)
                    ):
                        io1 = interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals)
                        io2 = interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals)
                        if io1 >= 0 and io2 >= 0:
                            assert self.amplicon_intervals[io1][3] == self.amplicon_intervals[io2][3]
                            bpi = self.addbp(
                                bp,
                                set(bpr),
                                bp_stats_,
                                self.amplicon_intervals[io1][3],
                            )
                            try:
                                self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))].add(bpi)
                            except:
                                self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))] = set([bpi])
                    else:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tDiscarded the subcluster %d." % (num_subcluster),
                        )
                    num_subcluster += 1
            else:
                logging.debug(
                    "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tDiscarded the cluster.",
                )

    def find_cn_breakpoints(self, b=300, n=50):
        """Search for breakpoints corresponding to CN boundaries"""
        cns_boundaries = []
        for ai in range(len(self.amplicon_intervals)):
            chr = self.amplicon_intervals[ai][0]
            s = self.amplicon_intervals[ai][1]
            e = self.amplicon_intervals[ai][2]
            si = list(self.pos2cni(chr, s))[0].data
            ei = list(self.pos2cni(chr, e))[0].data
            for i in range(si, ei):
                cns_boundaries.append(
                    (
                        ai,
                        chr,
                        self.cns_intervals_by_chr[chr][i][1],
                        self.cns_intervals_by_chr[chr][i][2],
                        self.cns_intervals_by_chr[chr][i + 1][2],
                    ),
                )
        cb_del_list = []
        for cbi in range(len(cns_boundaries)):
            cb = cns_boundaries[cbi]
            for bpi in range(len(self.new_bp_list)):
                bp = self.new_bp_list[bpi]
                if (bp[0] == cb[1] and bp[1] - 6001 < cb[3] < bp[1] + 6000) or (
                    bp[3] == cb[1] and bp[4] - 6001 < cb[3] < bp[4] + 6000
                ):
                    cb_del_list.append(cbi)
                    break
        for cbi in range(len(cns_boundaries)):
            if cbi not in cb_del_list:
                cb = cns_boundaries[cbi]
                nl, nr = n, n
                nl = min(nl, (cb[3] - cb[2] + 1) // b)
                nr = min(nr, (cb[4] - cb[3]) // b)
                cov_profile = [
                    sum(
                        [
                            sum(nc)
                            for nc in self.lr_bamfh.count_coverage(
                                cb[1],
                                cb[3] - (i + 1) * b + 1,
                                cb[3] - i * b + 1,
                            )
                        ],
                    )
                    * 1.0
                    / b
                    for i in range(nl)
                ][::-1] + [
                    sum(
                        [
                            sum(nc)
                            for nc in self.lr_bamfh.count_coverage(
                                cb[1],
                                cb[3] + i * b + 1,
                                cb[3] + (i + 1) * b + 1,
                            )
                        ],
                    )
                    * 1.0
                    / b
                    for i in range(nr)
                ]
                cb_refined = [-1, 0.0]
                for i in range(max(1, nl - 6000 // b), nl + min(nr - 1, 6000 // b)):
                    dmu = np.mean(cov_profile[:i]) - np.mean(cov_profile[i:])
                    if abs(dmu) > abs(cb_refined[1]):
                        cb_refined[0] = i
                        cb_refined[1] = dmu
                pval = 1.0
                cov_profile_split = [cov_profile[: cb_refined[0]], cov_profile[cb_refined[0] :]]
                if len(cov_profile_split[0]) > 1 and len(cov_profile_split[1]) > 1:
                    pval = stats.ttest_ind(
                        cov_profile_split[0],
                        cov_profile_split[1],
                        equal_var=False,
                    )[1]
                elif len(cov_profile_split[0]) == 1:
                    zscore = abs(cov_profile_split[0][0] - np.mean(cov_profile)) / np.std(
                        cov_profile,
                    )
                    pval = stats.norm.sf(zscore)
                elif len(cov_profile_split[1]) == 1:
                    zscore = abs(cov_profile_split[1][0] - np.mean(cov_profile)) / np.std(
                        cov_profile,
                    )
                    pval = stats.norm.sf(zscore)
                # only keep obvious cn boundaries
                if pval <= 0.01 and abs(cb_refined[1]) >= 3 * self.normal_cov:
                    if cb_refined[0] < nl:
                        self.source_edges.append(
                            [
                                "source",
                                -1,
                                "-",
                                cb[1],
                                cb[3] - (nl - cb_refined[0]) * b,
                                "+",
                                abs(cb_refined[1]),
                            ],
                        )
                    else:
                        self.source_edges.append(
                            [
                                "source",
                                -1,
                                "-",
                                cb[1],
                                cb[3] + (cb_refined[0] - nl) * b,
                                "+",
                                abs(cb_refined[1]),
                            ],
                        )
                    if cb_refined[1] < 0:
                        self.source_edges[-1][4] += 1
                        self.source_edges[-1][5] = "-"
                    self.source_edge_ccids.append(self.amplicon_intervals[cb[0]][3])

    def build_graph(self):
        """Organize the identified discordant edges into a list of breakpoint graphs, stored in lr_graph
        Each graph represent a connected component of amplified intervals, i.e., amplicon
        """
        # Split amplified intervals according to discordant edges
        split_int = dict()
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            for ai in range(len(self.amplicon_intervals)):
                seg = self.amplicon_intervals[ai]
                if bp[0] == seg[0] and seg[1] < bp[1] < seg[2]:
                    if bp[2] == "+":
                        try:
                            split_int[ai].append((bp[1], bp[1] + 1, bpi, 1, "+"))
                        except:
                            split_int[ai] = [(bp[1], bp[1] + 1, bpi, 1, "+")]
                    if bp[2] == "-":
                        try:
                            split_int[ai].append((bp[1] - 1, bp[1], bpi, 1, "-"))
                        except:
                            split_int[ai] = [(bp[1] - 1, bp[1], bpi, 1, "-")]
                if bp[3] == seg[0] and seg[1] < bp[4] < seg[2]:
                    if bp[5] == "+":
                        try:
                            split_int[ai].append((bp[4], bp[4] + 1, bpi, 4, "+"))
                        except:
                            split_int[ai] = [(bp[4], bp[4] + 1, bpi, 4, "+")]
                    if bp[5] == "-":
                        try:
                            split_int[ai].append((bp[4] - 1, bp[4], bpi, 4, "-"))
                        except:
                            split_int[ai] = [(bp[4] - 1, bp[4], bpi, 4, "-")]

        # Split amplified intervals according to source edges
        for srci in range(len(self.source_edges)):
            srce = self.source_edges[srci]
            for ai in range(len(self.amplicon_intervals)):
                seg = self.amplicon_intervals[ai]
                if srce[3] == seg[0] and seg[1] < srce[4] < seg[2]:
                    if srce[5] == "+":
                        try:
                            split_int[ai].append(
                                (srce[4], srce[4] + 1, len(self.new_bp_list) + srci, 4, "+"),
                            )
                        except:
                            split_int[ai] = [
                                (srce[4], srce[4] + 1, len(self.new_bp_list) + srci, 4, "+"),
                            ]
                    if srce[5] == "-":
                        try:
                            split_int[ai].append(
                                (srce[4] - 1, srce[4], len(self.new_bp_list) + srci, 4, "-"),
                            )
                        except:
                            split_int[ai] = [
                                (srce[4] - 1, srce[4], len(self.new_bp_list) + srci, 4, "-"),
                            ]

        # Construct graphs with sequence and concordant edges
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Will split the following %d amplicon intervals into sequence edges and build breakpoint graphs."
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
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "Will split the amplicon interval at index %d." % (ai),
            )
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "\tSplit interval at %s." % (split_int[ai]),
            )
        for ai in split_int:
            split_int[ai].sort(key=lambda item: item[0])
            # Trim amplified intervals if possible
            """
			if self.amplicon_intervals[ai][2] - self.amplicon_intervals[ai][1] > self.max_seq_len:
				if split_int[ai][0][0] - self.amplicon_intervals[ai][1] > self.interval_delta:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - state_provider.TSTART) + \
							"Modified the start coordinate of interval at index %d." %(ai))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - state_provider.TSTART) + "\t%s:%d->%d." %(self.amplicon_intervals[ai][0],
							self.amplicon_intervals[ai][1], split_int[ai][0][0] - self.interval_delta))
					self.amplicon_intervals[ai][1] = split_int[ai][0][0] - self.interval_delta
				if self.amplicon_intervals[ai][2] - split_int[ai][-1][1] > self.interval_delta:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - state_provider.TSTART) + \
							"Modified the end coordinate of interval at index %d." %(ai))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - state_provider.TSTART) + "\t%s:%d->%d." %(self.amplicon_intervals[ai][0],
							self.amplicon_intervals[ai][2], split_int[ai][-1][1] + self.interval_delta))
					self.amplicon_intervals[ai][2] = split_int[ai][-1][1] + self.interval_delta
			"""
            sseg = self.amplicon_intervals[ai]
            amplicon_idx = self.ccid2id[sseg[3]] - 1
            for ssi in range(len(split_int[ai])):
                if ssi == 0:
                    self.lr_graph[amplicon_idx].add_node((sseg[0], sseg[1], "-"))
                    self.lr_graph[amplicon_idx].add_node((sseg[0], split_int[ai][ssi][0], "+"))
                    self.lr_graph[amplicon_idx].add_node((sseg[0], split_int[ai][ssi][1], "-"))
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
                    self.lr_graph[amplicon_idx].add_node((sseg[0], split_int[ai][ssi - 1][1], "-"))
                    self.lr_graph[amplicon_idx].add_node((sseg[0], split_int[ai][ssi][0], "+"))
                    self.lr_graph[amplicon_idx].add_node((sseg[0], split_int[ai][ssi][1], "-"))
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
            self.lr_graph[amplicon_idx].add_node((sseg[0], split_int[ai][-1][1], "-"))
            self.lr_graph[amplicon_idx].add_node((sseg[0], sseg[2], "+"))
            self.lr_graph[amplicon_idx].add_sequence_edge(sseg[0], split_int[ai][-1][1], sseg[2])
        for ai in range(len(self.amplicon_intervals)):
            if ai not in split_int:
                sseg = self.amplicon_intervals[ai]
                amplicon_idx = self.ccid2id[sseg[3]] - 1
                self.lr_graph[amplicon_idx].add_node((sseg[0], sseg[1], "-"))
                self.lr_graph[amplicon_idx].add_node((sseg[0], sseg[2], "+"))
                self.lr_graph[amplicon_idx].add_sequence_edge(sseg[0], sseg[1], sseg[2])
        for amplicon_idx in range(len(self.lr_graph)):
            self.lr_graph[amplicon_idx].sort_edges()

        # Add nodes corresponding to interval ends
        for ai in range(len(self.amplicon_intervals)):
            sseg = self.amplicon_intervals[ai]
            amplicon_idx = self.ccid2id[sseg[3]] - 1
            self.lr_graph[amplicon_idx].amplicon_intervals.append([sseg[0], sseg[1], sseg[2]])
            self.lr_graph[amplicon_idx].add_endnode((sseg[0], sseg[1], "-"))
            self.lr_graph[amplicon_idx].add_endnode((sseg[0], sseg[2], "+"))
        logging.debug(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "The following nodes correspond to interval ends.",
        )
        for amplicon_idx in range(len(self.lr_graph)):
            for node in self.lr_graph[amplicon_idx].endnodes.keys():
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "\tAmplicon %d, node %s." % (amplicon_idx + 1, str(node)),
                )

        # Construct graphs with discordant and source edges
        for bpi in range(len(self.new_bp_list)):
            bp = self.new_bp_list[bpi]
            bp_ccid = self.new_bp_ccids[bpi]
            io1 = interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals)
            io2 = interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals)
            assert self.amplicon_intervals[io1][3] == self.amplicon_intervals[io2][3]
            amplicon_idx = self.ccid2id[self.amplicon_intervals[io1][3]] - 1
            if self.amplicon_intervals[io1][3] != bp_ccid:
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "Reset the ccid for breakpoint %s at index %d from %d to %d."
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
            self.lr_graph[amplicon_idx].add_source_edge(srce[3], srce[4], srce[5])

        # Print summary statistics for each amplicon
        for amplicon_idx in range(len(self.lr_graph)):
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "Num sequence edges in amplicon %d = %d."
                % (amplicon_idx + 1, len(self.lr_graph[amplicon_idx].sequence_edges)),
            )
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "Num concordant edges in amplicon %d = %d."
                % (amplicon_idx + 1, len(self.lr_graph[amplicon_idx].concordant_edges)),
            )
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "Num discordant edges in amplicon %d = %d."
                % (amplicon_idx + 1, len(self.lr_graph[amplicon_idx].discordant_edges)),
            )
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "Num source edges in amplicon %d = %d."
                % (amplicon_idx + 1, len(self.lr_graph[amplicon_idx].source_edges)),
            )

    def assign_cov(self):
        """Extract the long read coverage from bam file, if missing, for each sequence edge"""
        for amplicon_idx in range(len(self.lr_graph)):
            for seqi in range(len(self.lr_graph[amplicon_idx].sequence_edges)):
                seg = self.lr_graph[amplicon_idx].sequence_edges[seqi]
                if seg[5] == -1:
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Finding LR cov for sequence edge %s." % (seg[:3]),
                    )
                    """
					For long read, use the total number of nucleotides
					"""
                    rl_list = [
                        read for read in self.lr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if read.infer_read_length()
                    ]
                    self.lr_graph[amplicon_idx].sequence_edges[seqi][5] = len(rl_list)
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
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "LR cov assigned for sequence edge %s." % (seg[:3]),
                    )
        """
		Extract the long read coverage from bam file, if missing, for each concordant edge 
		"""
        for amplicon_idx in range(len(self.lr_graph)):
            for eci in range(len(self.lr_graph[amplicon_idx].concordant_edges)):
                ec = self.lr_graph[amplicon_idx].concordant_edges[eci]
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "Finding cov for concordant edge %s." % (ec[:6]),
                )
                rls = set(
                    [read.query_name for read in self.lr_bamfh.fetch(contig=ec[0], start=ec[1], stop=ec[1] + 1)],
                )
                rrs = set(
                    [read.query_name for read in self.lr_bamfh.fetch(contig=ec[3], start=ec[4], stop=ec[4] + 1)],
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
                for bpi in self.lr_graph[amplicon_idx].nodes[(ec[0], ec[1], ec[2])][2]:
                    for r in self.lr_graph[amplicon_idx].discordant_edges[bpi][10]:
                        rbps.add(r[0])
                for bpi in self.lr_graph[amplicon_idx].nodes[(ec[3], ec[4], ec[5])][2]:
                    for r in self.lr_graph[amplicon_idx].discordant_edges[bpi][10]:
                        rbps.add(r[0])
                self.lr_graph[amplicon_idx].concordant_edges[eci][9] = rls | rrs
                self.lr_graph[amplicon_idx].concordant_edges[eci][8] = len(
                    (rls & rrs & rls1 & rrs1) - rbps,
                )
                logging.debug(
                    "#TIME "
                    + "%.4f\t" % (time.time() - state_provider.TSTART)
                    + "LR cov assigned for concordant edge %s." % (ec[:6]),
                )

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
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "There are %d reads in total covering at least one breakpoint in amplicon %d."
                % (len(bp_reads), amplicon_idx + 1),
            )

            for rn in bp_reads:
                bp_reads_rn = bp_reads[rn][0]
                bp_reads_rn_sdel = bp_reads[rn][1]
                paths = []
                if len(bp_reads_rn) == 1 and len(bp_reads_rn_sdel) == 0:
                    rints = [aint[:4] for aint in self.chimeric_alignments[rn][1]]
                    ai1 = bp_reads_rn[0][0]
                    ai2 = bp_reads_rn[0][1]
                    bpi = bp_reads_rn[0][2]
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Read %s covers a single breakpoint." % rn,
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tAlignment intervals on reference = %s; mapq = %s; bp = (%d, %d, %d)"
                        % (rints, self.chimeric_alignments[rn][2], ai1, ai2, bpi),
                    )
                    path = chimeric_alignment_to_path_i(
                        self.lr_graph[amplicon_idx],
                        rints,
                        ai1,
                        ai2,
                        bpi,
                    )
                    paths.append(path)
                    logging.debug(
                        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tResulting subpath = %s" % path,
                    )
                elif len(bp_reads_rn) > 1 and len(bp_reads_rn_sdel) == 0:
                    bp_reads_rn = sorted(bp_reads_rn, key=lambda item: min(item[0], item[1]))
                    bp_reads_rn_split = [[0]]
                    last_ai = max(bp_reads_rn[0][0], bp_reads_rn[0][1])
                    for i in range(1, len(bp_reads_rn)):
                        if min(bp_reads_rn[i][0], bp_reads_rn[i][1]) == last_ai:
                            bp_reads_rn_split[-1].append(i)
                        else:
                            bp_reads_rn_split.append([i])
                        last_ai = max(bp_reads_rn[i][0], bp_reads_rn[i][1])
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Read %s covers multiple breakpoints." % rn,
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tBlocks of local alignments: %s" % bp_reads_rn_split,
                    )
                    qints = self.chimeric_alignments[rn][0]
                    skip = 0
                    for qi in range(len(qints) - 1):
                        if qints[qi + 1][0] - qints[qi][1] < -self.min_bp_match_cutoff_:
                            skip = 1
                            break
                    if skip == 1:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tDiscarded the read due to overlapping local alignments.",
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on reference = %s; mapq = %s."
                            % (self.chimeric_alignments[rn][1], self.chimeric_alignments[rn][2]),
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on the read = %s." % qints,
                        )
                        continue
                    for ai_block in bp_reads_rn_split:
                        rints = [aint[:4] for aint in self.chimeric_alignments[rn][1]]
                        ai_list = [bp_reads_rn[bi][:2] for bi in ai_block]
                        bp_list = [bp_reads_rn[bi][2] for bi in ai_block]
                        if len(set(bp_list)) < len(bp_list):
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\tDiscarded the block due to repeated breakpoints.",
                            )
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\tBlocks of local alignments: %s" % ai_block,
                            )
                            continue
                        path = chimeric_alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            rints,
                            ai_list,
                            bp_list,
                        )
                        paths.append(path)
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on reference = %s; mapq = %s; bps = %s"
                            % (rints, self.chimeric_alignments[rn][2], bp_reads_rn),
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tResulting subpath = %s" % path,
                        )
                elif len(bp_reads_rn) == 0 and len(bp_reads_rn_sdel) == 1:
                    rints = self.large_indel_alignments[rn][0]
                    rq = rints[-1]
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Read %s covers a single small del breakpoint." % rn,
                    )
                    if rints[3] < rints[4]:
                        if rints[2] < rints[1]:
                            rints = [
                                [rints[0], rints[3], rints[2], "+"],
                                [rints[0], rints[1], rints[4], "+"],
                            ]
                        else:
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\tDiscarded the read due to inconsistent alignment information.",
                            )
                            continue
                    elif rints[2] > rints[1]:
                        rints = [
                            [rints[0], rints[3], rints[2], "-"],
                            [rints[0], rints[1], rints[4], "-"],
                        ]
                    else:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tDiscarded the read due to inconsistent alignment information.",
                        )
                        continue
                    bpi = bp_reads_rn_sdel[0][2]
                    path = []
                    if rints[0][3] == "+":
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on reference = %s; mapq = %s; bp = (1, 0, %d)" % (rints, rq, bpi),
                        )
                        path = chimeric_alignment_to_path_i(
                            self.lr_graph[amplicon_idx],
                            rints,
                            1,
                            0,
                            bpi,
                        )
                        paths.append(path)
                    else:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on reference = %s; mapq = %s; bp = (0, 1, %d)" % (rints, rq, bpi),
                        )
                        path = chimeric_alignment_to_path_i(
                            self.lr_graph[amplicon_idx],
                            rints,
                            0,
                            1,
                            bpi,
                        )
                        paths.append(path)
                    logging.debug(
                        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "\tResulting subpath = %s" % path,
                    )
                elif len(bp_reads_rn) == 0 and len(bp_reads_rn_sdel) > 1:
                    rints = self.large_indel_alignments[rn]
                    rq = rints[0][-1]
                    rints_ = set(
                        [(rint[0], min(rint[3], rint[4]), max(rint[3], rint[4])) for rint in rints],
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Read %s covers multiple small del breakpoints." % rn,
                    )
                    if len(rints_) > 1 or len(rints) <= 1:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tDiscarded the read due to inconsistent alignment information.",
                        )
                        continue
                    rints_ = [[rint[0], min(rint[3], rint[4]), max(rint[3], rint[4]), "+"] for rint in rints]
                    rints = sorted(rints, key=lambda item: min(item[1], item[2]))
                    for ri in range(len(rints)):
                        rint = rints[ri]
                        rints_.append([rint[0], min(rint[3], rint[4]), max(rint[3], rint[4]), "+"])
                        rints_[ri][2] = min(rint[1], rint[2])
                        rints_[ri + 1][1] = max(rint[1], rint[2])
                    bp_reads_rn_sdel_split = [[]]
                    bp_reads_rn_sdel = sorted(bp_reads_rn_sdel, key=lambda item: item[0])
                    last_ai = 0
                    for i in range(len(bp_reads_rn_sdel)):
                        if i == 0 or bp_reads_rn_sdel[i][0] == last_ai + 1:
                            bp_reads_rn_sdel_split[-1].append(i)
                        else:
                            bp_reads_rn_sdel_split.append([i])
                        last_ai = bp_reads_rn_sdel[i][0]
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tBlocks of local alignments: %s" % bp_reads_rn_sdel_split,
                    )
                    for ai_block in bp_reads_rn_sdel_split:
                        ai_list = [[bp_reads_rn_sdel[bi][0], bp_reads_rn_sdel[bi][0] + 1] for bi in ai_block]
                        bp_list = [bp_reads_rn_sdel[bi][2] for bi in ai_block]
                        if len(set(bp_list)) < len(bp_list):
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\tDiscarded the block due to repeated breakpoints.",
                            )
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\tBlocks of local alignments: %s" % ai_block,
                            )
                            continue
                        path = chimeric_alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            rints_,
                            ai_list,
                            bp_list,
                        )
                        paths.append(path)
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on reference = %s; mapq = %s; bps = %s"
                            % (rints_, rq, bp_reads_rn_sdel),
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tResulting subpath = %s" % path,
                        )
                else:
                    rints = [aint[:4] for aint in self.chimeric_alignments[rn][1]]
                    rints_ = self.large_indel_alignments[rn]
                    rint_split = []
                    skip = 0
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "Read %s covers breakpoints and small del breakpoints." % rn,
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tAlignment intervals on reference = %s; mapq = %s"
                        % (rints, self.chimeric_alignments[rn][2]),
                    )
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tSmall del alignment intervals on reference = %s" % rints_,
                    )
                    for rint_ in rints_:
                        fount_split_rint = 0
                        for ri in range(len(rints)):
                            rint = rints[ri]
                            if (
                                rint_[0] == rint[0]
                                and min(rint_[1], rint_[2]) > min(rint[1], rint[2])
                                and max(rint_[1], rint_[2]) < max(rint[1], rint[2])
                            ):
                                fount_split_rint = 1
                                rint_split.append(ri)
                                break
                        if fount_split_rint == 0:
                            skip = 1
                            break
                    if skip == 1:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tDiscarded the read due to inconsistent alignment information.",
                        )
                        continue
                    for rsi in range(len(rint_split)):
                        ri = rint_split[rsi]
                        rints.insert(ri, rints[ri][:])
                        if rints[ri][3] == "+":
                            rints[ri][2] = min(rints_[rsi][1], rints_[rsi][2])
                            rints[ri + 1][1] = max(rints_[rsi][1], rints_[rsi][2])
                        else:
                            rints[ri][2] = max(rints_[rsi][1], rints_[rsi][2])
                            rints[ri + 1][1] = min(rints_[rsi][1], rints_[rsi][2])
                        for i in range(len(bp_reads_rn)):
                            if bp_reads_rn[i][0] >= ri and bp_reads_rn[i][1] >= ri:
                                bp_reads_rn[i][0] += 1
                                bp_reads_rn[i][1] += 1
                        for i in range(len(bp_reads_rn_sdel)):
                            if bp_reads_rn_sdel[i][0] == rsi:
                                if rints[ri][3] == "+":
                                    bp_reads_rn.append([ri + 1, ri, bp_reads_rn_sdel[i][2]])
                                else:
                                    bp_reads_rn.append([ri, ri + 1, bp_reads_rn_sdel[i][2]])
                    bp_reads_rn = sorted(bp_reads_rn, key=lambda item: min(item[0], item[1]))
                    bp_reads_rn_split = [[0]]
                    last_ai = max(bp_reads_rn[0][0], bp_reads_rn[0][1])
                    for i in range(1, len(bp_reads_rn)):
                        if min(bp_reads_rn[i][0], bp_reads_rn[i][1]) == last_ai:
                            bp_reads_rn_split[-1].append(i)
                        else:
                            bp_reads_rn_split.append([i])
                        last_ai = max(bp_reads_rn[i][0], bp_reads_rn[i][1])
                    logging.debug(
                        "#TIME "
                        + "%.4f\t" % (time.time() - state_provider.TSTART)
                        + "\tBlocks of local alignments: %s" % bp_reads_rn_split,
                    )
                    qints = self.chimeric_alignments[rn][0]
                    skip = 0
                    for qi in range(len(qints) - 1):
                        if qints[qi + 1][0] - qints[qi][1] < -self.min_bp_match_cutoff_:
                            skip = 1
                            break
                    if skip == 1:
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tDiscarded the read due to overlapping local alignments.",
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on reference = %s; mapq = %s."
                            % (self.chimeric_alignments[rn][1], self.chimeric_alignments[rn][2]),
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on the read = %s." % qints,
                        )
                        continue
                    for ai_block in bp_reads_rn_split:
                        ai_list = [bp_reads_rn[bi][:2] for bi in ai_block]
                        bp_list = [bp_reads_rn[bi][2] for bi in ai_block]
                        if len(set(bp_list)) < len(bp_list):
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\tDiscarded the block due to repeated breakpoints.",
                            )
                            logging.debug(
                                "#TIME "
                                + "%.4f\t" % (time.time() - state_provider.TSTART)
                                + "\tBlocks of local alignments: %s" % ai_block,
                            )
                            continue
                        path = chimeric_alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            rints,
                            ai_list,
                            bp_list,
                        )
                        paths.append(path)
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tAlignment intervals on reference = %s; mapq (unsplit) = %s; bps = %s"
                            % (rints, self.chimeric_alignments[rn][2], bp_reads_rn),
                        )
                        logging.debug(
                            "#TIME "
                            + "%.4f\t" % (time.time() - state_provider.TSTART)
                            + "\tResulting subpath = %s" % path,
                        )
                for pi in range(len(paths)):
                    path = paths[pi]
                    if len(path) > 5 and valid_path(self.lr_graph[amplicon_idx], path):
                        if path in self.path_constraints[amplicon_idx][0]:
                            pci = self.path_constraints[amplicon_idx][0].index(path)
                            self.path_constraints[amplicon_idx][1][pci] += 1
                        elif path[::-1] in self.path_constraints[amplicon_idx][0]:
                            pci = self.path_constraints[amplicon_idx][0].index(path[::-1])
                            self.path_constraints[amplicon_idx][1][pci] += 1
                        else:
                            self.path_constraints[amplicon_idx][0].append(path)
                            self.path_constraints[amplicon_idx][1].append(1)
                            self.path_constraints[amplicon_idx][2].append(amplicon_idx)
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "There are %d distinct subpaths due to reads involving breakpoints in amplicon %d."
                % (len(self.path_constraints[amplicon_idx][0]), amplicon_idx + 1),
            )

            # Extract reads in concordant_edges_reads
            lc = len(self.lr_graph[amplicon_idx].concordant_edges)
            for ci in range(lc):
                for rn in self.lr_graph[amplicon_idx].concordant_edges[ci][9]:
                    if rn not in self.large_indel_alignments and rn not in self.chimeric_alignments:
                        concordant_reads[rn] = amplicon_idx
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "There are %d concordant reads within amplicon intervals in amplicon %d."
                % (len(concordant_reads), amplicon_idx + 1),
            )
            for aint in self.amplicon_intervals:
                if amplicon_idx != self.ccid2id[aint[3]] - 1:
                    continue
                for read in self.lr_bamfh.fetch(aint[0], aint[1], aint[2] + 1):
                    rn = read.query_name
                    q = read.mapq
                    if q >= 20 and rn in concordant_reads:
                        path = alignment_to_path(
                            self.lr_graph[amplicon_idx],
                            [read.reference_name, read.reference_start, read.reference_end],
                        )
                        if len(path) > 5 and valid_path(self.lr_graph[amplicon_idx], path):
                            if path in self.path_constraints[amplicon_idx][0]:
                                pci = self.path_constraints[amplicon_idx][0].index(path)
                                self.path_constraints[amplicon_idx][1][pci] += 1
                            elif path[::-1] in self.path_constraints[amplicon_idx][0]:
                                pci = self.path_constraints[amplicon_idx][0].index(path[::-1])
                                self.path_constraints[amplicon_idx][1][pci] += 1
                            else:
                                self.path_constraints[amplicon_idx][0].append(path)
                                self.path_constraints[amplicon_idx][1].append(1)
                                self.path_constraints[amplicon_idx][2].append(concordant_reads[rn])
            logging.debug(
                "#TIME "
                + "%.4f\t" % (time.time() - state_provider.TSTART)
                + "There are %d distinct subpaths in total in amplicon %d."
                % (len(self.path_constraints[amplicon_idx][0]), amplicon_idx + 1),
            )

    def closebam(self):
        """Close the short read and long read bam file"""
        self.lr_bamfh.close()


def reconstruct_graph(
    lr_bam: str,
    cnv_seed: str,
    output_prefix: str,
    cn_seg: str,
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
    start_time = time.time()
    state_provider.TSTART = start_time
    log_file = output_prefix or "infer_breakpoint_graph.log"
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.DEBUG,
        format="[%(name)s:%(levelname)s]\t%(message)s",
    )
    logging.info("Python version " + sys.version + "\n")
    commandstring = "Commandline: "
    for arg in sys.argv:
        if " " in arg:
            commandstring += f'"{arg}" '
        else:
            commandstring += f"{arg} "
    logging.info("#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + commandstring)

    b2bn = BamToBreakpointNanopore(lr_bam, cnv_seed)
    # if filter_bp_by_edit_distance:
    # b2bn.nm_filter = True
    b2bn.min_bp_cov_factor = min_bp_support
    logging.info("#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Opened LR bam files.")
    b2bn.read_cns(cn_seg)
    logging.info(
        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Completed parsing CN segment files.",
    )
    b2bn.fetch()
    logging.info(
        "#TIME "
        + "%.4f\t" % (time.time() - state_provider.TSTART)
        + "Completed fetching reads containing breakpoints.",
    )
    b2bn.hash_alignment_to_seg()
    logging.info(
        "#TIME "
        + "%.4f\t" % (time.time() - state_provider.TSTART)
        + "Completed hashing chimeric reads to CN segments.",
    )
    b2bn.find_amplicon_intervals()
    logging.info(
        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Completed finding amplicon intervals.",
    )
    b2bn.find_smalldel_breakpoints()
    logging.info(
        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Completed finding small del breakpoints.",
    )
    b2bn.find_breakpoints()
    logging.info(
        "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Completed finding all discordant breakpoints.",
    )
    if output_bp:
        b2bn.build_graph()
        logging.info(
            "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Breakpoint graph built for all amplicons.",
        )
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
            breakpoint_graph.output_breakpoint_info_lr(
                b2bn.lr_graph[gi],
                output_prefix + "_amplicon" + str(gi + 1) + "_breakpoints.txt",
                bp_stats_i,
            )
        logging.info(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Wrote breakpoint information, for all amplicons, to %s."
            % (output_prefix + "_amplicon*_breakpoints.txt"),
        )
    else:
        # b2bn.find_cn_breakpoints()
        # logging.info("#TIME " + '%.4f\t' %(time.time() - state_provider.TSTART) + "Completed finding breakpoints corresponding to CN changes.")
        b2bn.build_graph()
        logging.info(
            "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Breakpoint graph built for all amplicons.",
        )
        b2bn.assign_cov()
        logging.info(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Fetched read coverage for all sequence and concordant edges.",
        )
        for gi in range(len(b2bn.lr_graph)):
            b2bn.lr_graph[gi].compute_cn_lr(b2bn.normal_cov)
        logging.info(
            "#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Computed CN for all edges.",
        )
        for gi in range(len(b2bn.lr_graph)):
            breakpoint_graph.output_breakpoint_graph_lr(
                b2bn.lr_graph[gi],
                output_prefix + "_amplicon" + str(gi + 1) + "_graph.txt",
            )
        logging.info(
            "#TIME "
            + "%.4f\t" % (time.time() - state_provider.TSTART)
            + "Wrote breakpoint graph for all complicons to %s." % (output_prefix + "_amplicon*_graph.txt"),
        )
    return b2bn


def print_complete_message():
    logging.info("#TIME " + "%.4f\t" % (time.time() - state_provider.TSTART) + "Total runtime.")
