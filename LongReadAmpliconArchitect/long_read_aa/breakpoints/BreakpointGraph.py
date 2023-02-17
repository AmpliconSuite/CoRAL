"""
Implements a class for BreakpointGraph. Will serve as a container for the 
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file. 
"""
from typing import List, Union

import pandas as pd

from . import BreakpointEdge, SequenceEdge


class BreakpointGraph:
    """A container object for the breakpoint graphs. 
    """
    def __init__(
        self,
        sequence_edges: List[SequenceEdge.SequenceEdge] = [],
        discordant_edges: List[BreakpointEdge.BreakpointEdge] = [],
        concordant_edges: List[BreakpointEdge.BreakpointEdge] = [],
        source_edges: List[BreakpointEdge.BreakpointEdge] = []
    ):
        
        # Keep these private so they're protected
        self.__sequence_edges = sequence_edges
        self.__discordant_edges = discordant_edges
        self.__concordant_edges = concordant_edges
        self.__source_edges = source_edges
    
    @property
    def sequence_edges(self):
        return self.__sequence_edges

    @property
    def discordant_edges(self):
        return self.__discordant_edges

    @property
    def concordant_edges(self):
        return self.__concordant_edges

    @property
    def source_edges(self):
        return self.__source_edges

    def add_edge(self, edge: Union[BreakpointEdge.BreakpointEdge, SequenceEdge.SequenceEdge], annotation: str):
        """Adds an edge to the appropriate list.
        """
        pass

    def get_breakpoint_graph_dataframe(self) -> pd.DataFrame:
        """Obtain a DataFrame rendering of the breakpoint graph.
        """
        pass

