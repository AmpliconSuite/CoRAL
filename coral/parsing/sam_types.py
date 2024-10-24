from __future__ import annotations

import functools
import io
import logging
import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, NamedTuple

import intervaltree
import numpy as np
import pysam

import coral.cigar_parsing
import coral.parsing.sam_utils
from coral.parsing import sam_utils
from coral.parsing.cns_types import CnsInterval

logger = logging.getLogger(__name__)


@dataclass
class EditDistTracker:
    total: float = 0.0
    sum_of_squares: float = 0.0
    count: int = 0

    def observe(self, value: float) -> None:
        self.total += value
        self.sum_of_squares += value**2
        self.count += 1


class EditDistStats(NamedTuple):
    mean: float
    std_dev: float
    length: int

    @staticmethod
    def from_tracker(tracker: EditDistTracker) -> EditDistStats:
        mean = tracker.total / tracker.count
        return EditDistStats(
            mean=tracker.total / tracker.count,
            std_dev=math.sqrt(tracker.sum_of_squares / tracker.count - mean**2),
            length=tracker.count,
        )


@dataclass
class ChimericBamData:
    # Map read name -> chimeric alignment (containing two or more records for one read)
    chimeric_alignments_by_name: dict[str, sam_utils.ChimericAlignment]

    # NM = "Edit distance to the reference, including ambiguous bases but excluding clipping"
    # https://www.samformat.info/sam-format-alignment-tags
    nm_stats: EditDistStats

    @classmethod
    def from_bam(cls, lr_bamfh: pysam.AlignmentFile) -> ChimericBamData:
        read_name_to_length: Dict[str, int] = {}
        raw_chimeric_reads: Dict[str, list[sam_utils.ChimericAlignmentElement]] = defaultdict(list)
        nm_tracker = EditDistTracker()

        for read in lr_bamfh.fetch():
            if not (read_name := read.query_name):
                continue
            if read.flag < 256 and read_name not in read_name_to_length:
                read_name_to_length[read_name] = read.query_length
            if tag_alignments := sam_utils.parse_chimeric_tag(read):
                raw_chimeric_reads[read_name].append(tag_alignments)
            elif read.mapping_quality == 60:  # Corresponds to unique read when using minimap # https://github.com/lh3/minimap2/issues/447
                edit_dist = read.get_cigar_stats()[0][-1] / read.query_length
                nm_tracker.observe(edit_dist)
        nm_stats = EditDistStats.from_tracker(nm_tracker)
        logger.info(f"Fetched {len(raw_chimeric_reads)} chimeric reads.")
        chimeric_alignments_by_name: Dict[str, sam_utils.ChimericAlignment] = {}

        # reads_wo_primary_alignment = []
        for read_name in raw_chimeric_reads:
            if read_name not in read_name_to_length:
                logger.warning(f"Found chimeric read name {read_name} without primary alignment; Read length: N/A.")
                logger.warning(f"All CIGAR strings: {raw_chimeric_reads[read_name]}.")
                # reads_wo_primary_alignment.append(read_name)
                continue
            read_length = read_name_to_length[read_name]
            chimeric_alignments_by_name[read_name] = sam_utils.chimeric_alignment_from_raw_elements(
                raw_chimeric_reads[read_name], read_length
            )
        # for read_name in reads_wo_primary_alignment:
        #     del chimeric_alignments_by_name[read_name]
        logger.info("Computed alignment intervals on all chimeric reads.")

        return cls(chimeric_alignments_by_name=chimeric_alignments_by_name, nm_stats=nm_stats)
