"""
Implements a class for BreakpointGraph. Will serve as a container for the 
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file. 
"""
from typing import List

from . import BreakpointEdge, SequenceEdge

class BreakpointGraph:
    """A container object for the breakpoint graphs. 
    """
    def __init__(
        self,
        sequence_edges: List[SequenceEdge.SequenceEdge],
        discordant_edges: List[BreakpointEdge.BreakpointEdge],
        concordant_edges: List[BreakpointEdge.BreakpointEdge],
    ):
        
        self.sequence_edges = sequence_edges
        self.discordant_edges = discordant_edges
        self.concordant_edges = concordant_edges

    def infer(self):
        pass
