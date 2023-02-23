"""
Container class for sequence edges.
"""
from typing import Dict, Optional, Tuple


class SequenceEdge:
    def __init__(
        self,
        chromosome: str,
        start: int,
        end: int,
        support: Optional[Dict[str, int]] = None,
        copy_number: Optional[Dict[str, int]] = None,
        average_read_length: Optional[int] = None,
    ):
        """Class for sequence edges.

        Sequence edges are contiguous nucleotide sequences spanning different
        breakpoints.

        Args:
            chromosome: Chromosome
            start: Start of sequence edge
            end: End of sequence edge
            support: Read count support dictionary. Can store
                support for arbitrarily many sequencing libraries.
            average_read_length: Average read length. Used for copy-number
                estimation.
            copy_number: Estimated copy number.
        """
        if end < start:
            raise Exception("Sequence end start is after end.")
        if support and any(support[key] < 0 for key in support):
            raise Exception("Support should be positive.")
        if copy_number and any(copy_number[key] < 0 for key in copy_number):
            raise Exception("Copy number should be positive.")

        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.support = support
        self.copy_number = copy_number
        self.average_read_length = average_read_length

    @property
    def length(self) -> int:
        """Return length of sequence.

        Note that both start and end node will be inclusive on a sequence edge.
        """
        return self.end - self.start + 1

    @property
    def left_breakpoint(self) -> Tuple[str, int, str]:
        """Quick accessor of the left breakpoint."""
        return (self.chromosome, self.start, "+")

    @property
    def right_breakpoint(self) -> Tuple[str, int, str]:
        """Quick accessor of the right breakpoint."""
        return (self.chromosome, self.end, "+")

    def __str__(self):
        """Prints SequenceEdge."""
        return f"{self.chromosome}:{self.start}-{self.end}"
