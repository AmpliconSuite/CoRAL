from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import DefaultDict, Dict, Literal, NamedTuple, cast

import pysam

from coral import types
from coral.cigar_parsing import CIGAR2POS_OPS, logger
from coral.parsing.cns_types import CNSegIntervals


class ChimericAlignmentElement(NamedTuple):
    "Represents a single canonical alignment element within an SA tag for a chimeric alignment."

    rname: str
    pos: int
    strand: types.Strand
    cigar: str
    mapq: int
    nm: float


class Range(NamedTuple):
    start: int
    end: int


@dataclass
class ReadCoordinate:
    name: str
    start: int
    end: int
    strand: types.Strand
    cn_idx_range: Range | None = None


@dataclass
class RawChimericSAMData:
    reads: list[pysam.AlignedSegment]
    read_length: int
    raw_cns: CNSegIntervals


def parse_chimeric_tag(read: pysam.AlignedSegment) -> list[ChimericAlignmentElement] | None:
    # SA tag corresponds to chimeric alignment, https://www.samformat.info/sam-format-alignment-tags
    # Elements in semicolon-delimited list are (rname, pos, strand, CIGAR, mapQ, NM;) that correspond to canonical alignment partitions
    try:
        chimeric_tag_value = read.get_tag("SA:Z")
        canonical_alignment_elements = [record.split(",") for record in chimeric_tag_value.split(";")[:-1]]  # type:ignore[union-attr]
        return [
            ChimericAlignmentElement(elem[0], int(elem[1]), cast(types.Strand, elem[2]), elem[3], int(elem[4]), float(elem[5]))
            for elem in canonical_alignment_elements
        ]
    except KeyError:  # KeyError raised by pysam on missing tag, faster than checking `has_tag` + then accessing
        return None


@dataclass
class ChimericAlignment:
    ranges: list[Range] = field(default_factory=list)
    coords: list[ReadCoordinate] = field(default_factory=list)
    quals: list[int] = field(default_factory=list)
    nms: list[float] = field(default_factory=list)

    def sort(self) -> None:
        sorted_items = sorted(zip(self.ranges, self.coords, self.quals, self.nms), key=lambda x: (x[0].start, x[0].end))
        self.ranges = [item[0] for item in sorted_items]
        self.coords = [item[1] for item in sorted_items]
        self.quals = [item[2] for item in sorted_items]
        self.nms = [item[3] / (item[0].end - item[0].start) for item in sorted_items]

    def __len__(self) -> int:
        return len(self.ranges)


def chimeric_alignment_from_raw_elements(canonical_elements: list[ChimericAlignmentElement], read_length: int) -> ChimericAlignment:
    """
    Convert list of ChimericAlignmentElements (sourced from canonical partitions of "SA:Z" tag) into a new chimeric alignment.
    Require at least one (soft) clip and one match for each canonical alignment record in a chimeric alignment
            If not, trigger a warning message in logger

    Args:
            sa_list: A list of ChimericAlignmentElements generated from a SAM "SA:Z" tag
            read_length: Read length
    Returns:
            Chimeric Alignment in the form of qint (alignment start/end), rint (read coordinates), qual (map quality), and nm (edit dist) lists
            Alignments sorted according to the starting positions on the read on positive strand

    """
    chimeric_alignment = ChimericAlignment()
    for element in canonical_elements:
        if "S" not in element.cigar or "M" not in element.cigar:
            # Require a chimeric alignment record having at least some (soft)clips and matches
            logger.warning("Found chimeric alignment without match or soft clips.")
            logger.warning(f"\tAll CIGAR strings: {canonical_elements}")
            return chimeric_alignment
        op = "".join(c for c in element.cigar if not c.isdigit())
        start, end, alignment_length = CIGAR2POS_OPS[op](element.cigar, element.strand, read_length)
        chimeric_alignment.ranges.append(Range(start, end))

        # Convert to 0 based coordinates
        if element.strand == "+":
            chimeric_alignment.coords.append(ReadCoordinate(element.rname, element.pos - 1, element.pos + alignment_length - 2, "+"))
        else:
            chimeric_alignment.coords.append(ReadCoordinate(element.rname, element.pos + alignment_length - 2, element.pos - 1, "-"))

        chimeric_alignment.quals.append(element.mapq)
        chimeric_alignment.nms.append(float(element.nm))
    chimeric_alignment.sort()
    return chimeric_alignment


def map_alignment_hash_to_cn_seg(
    chimeric_alignments_by_read_name: dict[str, ChimericAlignment], raw_cns: CNSegIntervals
) -> dict[str, dict[int, list[str]]]:
    """Speed up amplified interval search by hashing chimeric alignments from each long read to CN segments"""
    cn_seg_to_chimeric_reads: dict[str, dict[int, list[str]]] = defaultdict(lambda: defaultdict(list))
    for read_name, chimeric_alignment in chimeric_alignments_by_read_name.items():
        # breakpoint()
        for coordinates in chimeric_alignment.coords:
            if (chr_tag := coordinates.name) in raw_cns.tree:
                start_cn_intervals = raw_cns.get_cn_intervals(chr_tag, min(coordinates.start, coordinates.end))
                end_cn_intervals = raw_cns.get_cn_intervals(chr_tag, max(coordinates.start, coordinates.end))
                assert len(start_cn_intervals) <= 1 and len(end_cn_intervals) <= 1
                start_cn_idx, end_cn_idx = -1, -1
                if len(start_cn_intervals) == 1:
                    start_cn_idx = next(iter(start_cn_intervals)).data
                    cn_seg_to_chimeric_reads[chr_tag][start_cn_idx].append(read_name)
                if len(end_cn_intervals) == 1:
                    end_cn_idx = next(iter(start_cn_intervals)).data
                    cn_seg_to_chimeric_reads[chr_tag][end_cn_idx].append(read_name)
                cniset = {start_cn_idx, end_cn_idx}
                if len(cniset) > 1 and -1 in cniset:
                    cniset.remove(-1)
                coordinates.cn_idx_range = Range(start_cn_idx, end_cn_idx)
    logger.info("Completed hashing chimeric reads to CN segments.")
    return cn_seg_to_chimeric_reads
