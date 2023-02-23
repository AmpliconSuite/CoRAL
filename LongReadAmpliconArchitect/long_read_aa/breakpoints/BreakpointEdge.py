"""
A simple class to encode a breakpoint edge.
"""
from typing import Dict, Literal, Optional


class BreakpointEdge:
    """A container class for breakpoint edges."""

    def __init__(
        self,
        chromosome1: str,
        position1: int,
        orientation1: str,
        chromosome2: str,
        position2: int,
        orientation2: str,
        annotation: Optional[
            Literal["concordant", "discordant", "source"]
        ] = None,
        support: Optional[Dict[str, int]] = None,
        copy_number: Optional[Dict[str, int]] = None,
    ):
        """Constructor of BreakpointEdge.

        Args:
            chromosome1: Chromosome of left breakpoint
            position1: Position of left breakpoint
            orientation1: Orientation of left breakpoint ("+" or "-")
            chromosome2: Chromosome of right breakpoint
            position2: Position of right breakpoint
            orientation2: Orientation of right breakpoint ("+" or "-")
            annotation: String indicating the type of breakpoint
            support: Read number supporting this breakpoint
            copy_number: Estimated copy number of breakpoint
            amplicon_intervals: Amplicon intervals covered by this breakpoint graph.
        """
        if support and support < 0:
            raise Exception("Support should be positive.")
        if copy_number and copy_number < 0:
            raise Exception("Copy number should be positive.")

        self.chromosome1 = chromosome1
        self.position1 = position1
        self.orientation1 = orientation1
        self.chromosome2 = chromosome2
        self.position2 = position2
        self.orientation2 = orientation2
        self.annotation = annotation

        self.support = support
        self.copy_number = copy_number

        # store read names supporting each breakpoint edge. 
        self.read_names = []

    @property
    def left_breakpoint(self):
        """Quick accessor of the left breakpoint."""
        return (self.chromosome1, self.position1, self.orientation1)

    @property
    def right_breakpoint(self):
        """Quick accessor of the right breakpoint."""
        return (self.chromosome2, self.position2, self.orientation2)

    def __str__(self):
        """Print BreakpointEdge."""
        return(f"{self.left_breakpoint} -> {self.right_breakpoint}")
