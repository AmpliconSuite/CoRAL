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

from coral import cigar_parsing
from coral.types import CnsInterval

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
class ParsedBam:
    # Map read name -> chimeric alignments (two or more records for one read)
    chimeric_alignments: dict[str, list[pysam.AlignedSegment]]

    # For edit distance filter of breakpoints
    nm_stats: EditDistStats
    nm_filter: bool = False  # TODO: not set anywhere currently. what is NM?

    @classmethod
    def fetch(cls, lr_bamfh: pysam.AlignmentFile) -> ParsedBam:
        read_name_to_length: Dict[str, int] = {}
        chimeric_alignments: Dict[str, List[pysam.AlignedSegment]] = defaultdict(list)
        nm_tracker = EditDistTracker()

        for read in lr_bamfh.fetch():
            if not (read_name := read.query_name):
                continue
            if read.flag < 256 and read_name not in read_name_to_length:
                read_name_to_length[read_name] = read.query_length
            # SA tag corresponds to chimeric alignment, https://www.samformat.info/sam-format-alignment-tags
            chimeric_tag_value = read.get_tag("SA:Z")
            try:
                tag_chimeric_alignments: List[pysam.AlignedSegment] = chimeric_tag_value[:-1].split(";")
                for chimeric_alignment in tag_chimeric_alignments:
                    chimeric_alignments[read_name].append(chimeric_alignment)
            except:
                if read.mapping_quality == 60:
                    edit_dist = read.get_cigar_stats()[0][-1] / read.query_length
                    nm_tracker.observe(edit_dist)
        nm_stats = EditDistStats.from_tracker(nm_tracker)
        logger.info(f"Fetched {len(chimeric_alignments)} chimeric reads.")
        reads_wo_primary_alignment = []
        for read_name in chimeric_alignments:
            if read_name not in read_name_to_length:
                logger.warning(f"Found chimeric read name {read_name} without primary alignment; Read length: N/A.")
                logger.warning(f"All CIGAR strings: {chimeric_alignments[read_name]}.")
                reads_wo_primary_alignment.append(read_name)
                continue
            read_length = read_name_to_length[read_name]
            chimeric_alignments[read_name] = cigar_parsing.alignment_from_satags(self.chimeric_alignments[read_name], read_length)
        for read_name in reads_wo_primary_alignment:
            del chimeric_alignments[read_name]
        logger.info(
            "Computed alignment intervals on all chimeric reads.",
        )
        return cls(chimeric_alignments=chimeric_alignments, nm_stats=nm_stats)
