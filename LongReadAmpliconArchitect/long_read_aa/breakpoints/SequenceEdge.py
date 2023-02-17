"""
Container class for sequence edges.
"""
from typing import Optional

class SequenceEdge:
    def __init__(
        self,
        chromosome: str,
        start: int,
        end: int,
        support: Optional[int] = None,
        copy_number: Optional[float] = None,
    ):
        """Class for sequence edges.

        Sequence edges are contiguous nucleotide sequences spanning different
        breakpoints.

        Args:
            chromosome: Chromosome
            start: Start of sequence edge
            end: End of sequence edge
            support: Read count support
            copy_number: Estimated copy number.
        """
        if end < start:
            raise Exception("Sequence end start is after end.")
        if support and support < 0:
            raise Exception("Support should be positive.")
        if copy_number and copy_number < 0:
            raise Exception("Copy number should be positive.")
        
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.support = support
        self.copy_number = copy_number

    @property
    def length(self) -> int:
        """Return length of sequence."""
        return self.end - self.start

    def __str__(self):
        """Prints SequenceEdge."""
        return f"{self.chromosome}:{self.start}-{self.end}"
