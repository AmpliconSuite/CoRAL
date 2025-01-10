from __future__ import annotations

from typing import Generator, Protocol

import pysam

from coral.datatypes import Interval


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


# class TypedBAM(pysam.AlignmentFile, BAMProtocol):
#     pass
