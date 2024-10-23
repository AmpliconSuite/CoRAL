from __future__ import annotations

import io
import logging
from dataclasses import dataclass
from typing import Any, Dict, List

import intervaltree
import numpy as np
import pysam

from coral.types import CnsInterval

logger = logging.getLogger(__name__)


@dataclass
class RawCNData:
    tree: intervaltree.IntervalTree  # Interval tree structure for each chromosome
    intervals: List[CnsInterval]  # CN segments
    intervals_by_chr: Dict[str, List[List[int]]]
    log2_cn: List[float]  # log2 CN for each CN segment

    @classmethod
    def from_file(cls, file: io.TextIOWrapper) -> RawCNData:
        tree, log2_cns = {}, []
        intervals = []
        intervals_by_chr: Dict[str, List[List[Any]]] = {}
        is_cns = file.name.endswith(".cns")

        next(file)  # Skip header row
        idx = 0
        for line in file:
            fields = line.strip().split()
            chr_tag, start, end = fields[0], int(fields[1]), int(fields[2])
            intervals.append([chr_tag, start, end - 1])
            if chr_tag not in tree:
                tree[chr_tag] = intervaltree.IntervalTree()
                intervals_by_chr[chr_tag] = []
                idx = 0
            tree[chr_tag][start:end] = idx
            idx += 1

            # Calc log2 CN for .cns
            log2_cn = float(fields[4]) if is_cns else np.log2(float(fields[3]) / 2.0)
            raw_cn = 2 * (2**log2_cn) if is_cns else float(fields[3])

            log2_cns.append(log2_cn)
            intervals_by_chr[chr_tag].append([chr_tag, start, end - 1, raw_cn])

        return cls(tree, intervals, intervals_by_chr, log2_cns)


@dataclass
class CoverageInfo:
    cns: RawCNData

    cn_gain: float = 5.0
    normal_cov: float = 0.0  # Normal long read coverage - used for CN assignment

    @classmethod
    def from_cns(
        cls,
        cns_file: io.TextIOWrapper,
        bam_file: pysam.AlignmentFile,
        min_cluster_cutoff: int = 3,
        min_bp_cov_factor: float = 1.0,
        min_bp_match_cutoff: int = 100,
    ):
        """Read in (cnvkit) *.cns file and estimate the normal long read coverage"""
        cns = RawCNData.from_file(cns_file)

        logger.debug(f"Total num LR copy number segments:{len(cns.log2_cn)}.")
        log2_cn_order = np.argsort(cns.log2_cn)
        cns_intervals_median = []
        log2_cn_median = []
        im = int(len(log2_cn_order) / 2.4)
        ip = im + 1
        total_int_len = 0
        cns_intervals_median.append(cns.intervals[log2_cn_order[ip]])
        cns_intervals_median.append(cns.intervals[log2_cn_order[im]])
        log2_cn_median.append(cns.log2_cn[log2_cn_order[ip]])
        log2_cn_median.append(cns.log2_cn[log2_cn_order[im]])
        total_int_len += cns.intervals[log2_cn_order[ip]][2] - cns.intervals[log2_cn_order[ip]][1] + 1
        total_int_len += cns.intervals[log2_cn_order[im]][2] - cns.intervals[log2_cn_order[im]][1] + 1
        i = 1
        while total_int_len < 10000000:
            cns_intervals_median.append(cns.intervals[log2_cn_order[ip + i]])
            cns_intervals_median.append(cns.intervals[log2_cn_order[im - i]])
            log2_cn_median.append(cns.log2_cn[log2_cn_order[ip]])
            log2_cn_median.append(cns.log2_cn[log2_cn_order[im]])
            total_int_len += cns.intervals[log2_cn_order[ip + i]][2] - cns.intervals[log2_cn_order[ip + i]][1] + 1
            total_int_len += cns.intervals[log2_cn_order[im - i]][2] - cns.intervals[log2_cn_order[im - i]][1] + 1
            i += 1
        logger.debug(f"Use {len(cns_intervals_median)} LR copy number segments.")
        logger.debug(f"Total length of LR copy number segments: {total_int_len}.")
        logger.debug(f"Average LR copy number: {np.average(log2_cn_median)}.")
        nnc = 0
        for i in range(len(cns_intervals_median)):
            nnc += sum(
                [
                    sum(nc)
                    for nc in bam_file.count_coverage(
                        cns_intervals_median[i][0],
                        cns_intervals_median[i][1],
                        cns_intervals_median[i][2] + 1,
                        quality_threshold=0,
                        read_callback="nofilter",
                    )
                ]
            )
        normal_cov = nnc * 1.0 / total_int_len
        logger.info(f"LR normal cov ={normal_cov}, {nnc=}, {total_int_len=}.")
        min_cluster_cutoff = max(min_cluster_cutoff, int(min_bp_cov_factor * normal_cov))
        logger.debug(f"Reset min_cluster_cutoff to {min_cluster_cutoff}.")
