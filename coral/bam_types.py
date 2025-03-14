from __future__ import annotations

from typing import Generator, Protocol

import pysam

from coral.datatypes import Interval, Node


class BAMReadProtocol(Protocol):
    pass


class BAMRead(BAMReadProtocol, pysam.AlignedRead):  # type: ignore
    pass


class SAMSegmentProtocol(Protocol):
    # Abstract class used to add more informative typing to pysam BAM ops.
    # Also, .pyx docstrings don't appear to work natively with VSCode/Pylance.
    query_name: str  # Force non-optional (unlike pysam.AlignedSegment)

    def get_cigar_stats(self) -> tuple[list[int], list[int]]:
        """Summary of operations in cigar string.

        The output order in the array is "MIDNSHP=X" followed by a
        field for the NM tag. If the NM tag is not present, this
        field will always be 0.

        +-----+--------------+-----+
        |M    |BAM_CMATCH    |0    |
        +-----+--------------+-----+
        |I    |BAM_CINS      |1    |
        +-----+--------------+-----+
        |D    |BAM_CDEL      |2    |
        +-----+--------------+-----+
        |N    |BAM_CREF_SKIP |3    |
        +-----+--------------+-----+
        |S    |BAM_CSOFT_CLIP|4    |
        +-----+--------------+-----+
        |H    |BAM_CHARD_CLIP|5    |
        +-----+--------------+-----+
        |P    |BAM_CPAD      |6    |
        +-----+--------------+-----+
        |=    |BAM_CEQUAL    |7    |
        +-----+--------------+-----+
        |X    |BAM_CDIFF     |8    |
        +-----+--------------+-----+
        |B    |BAM_CBACK     |9    |
        +-----+--------------+-----+
        |NM   |NM tag        |10   |
        +-----+--------------+-----+

        If no cigar string is present, empty arrays will be returned.

        Returns:
            arrays :
                two arrays. The first contains the nucleotide counts within
                each cigar operation, the second contains the number of blocks
                for each cigar operation.

        """
        ...


class SAMSegment(SAMSegmentProtocol, pysam.AlignedSegment):  # type: ignore
    pass


class BAMWrapper(pysam.AlignmentFile):
    def fetch_interval(
        self, interval: Interval
    ) -> Generator[SAMSegment, None, None]:
        return self.fetch(interval.chr, interval.start, interval.end + 1)  # type: ignore

    def fetch_node(self, node: Node) -> Generator[SAMSegment, None, None]:
        return self.fetch(node.chr, node.pos, node.pos + 1)  # type: ignore

    def count_raw_coverage(
        self,
        interval: Interval,
        quality_threshold: int = 0,
        read_callback_type: str = "nofilter",
    ) -> float:
        # Much faster than pysam.count_coverage, ignore per-base coverage.
        if read_callback_type == "all":
            # Filter out unmapped, secondary, supplementary, and PCR/optical
            # duplicates.
            read_cb = lambda read: not (
                read.flag & (0x4 | 0x100 | 0x200 | 0x400)
            )
        else:
            read_cb = lambda _: True
        reads = [
            read
            for read in self.fetch_interval(interval)
            if read.mapping_quality >= quality_threshold and read_cb(read)
        ]
        total_read_coverage = 0
        for read in reads:
            for block_start, block_end in read.get_blocks():
                if block_end < interval.start or block_start > interval.end:
                    continue
                total_read_coverage += min(block_end, interval.end) - max(
                    block_start, interval.start
                )
        return total_read_coverage / (interval.end - interval.start)


# class TypedBAM(pysam.AlignmentFile, BAMProtocol):
#     pass
