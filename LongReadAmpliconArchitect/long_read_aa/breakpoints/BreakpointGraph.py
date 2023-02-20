"""
Implements a class for BreakpointGraph. Will serve as a container for the 
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file. 
"""
from typing import List, Union

import pandas as pd

from long_read_aa.breakpoints import BreakpointEdge, SequenceEdge


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
        self.__nodes = dict()
    
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

    def __getitem__(self, key: str):
        """Access the adjacency list of the breakpoint graph.

        Args:
            key: A breakpoint node.

        Returns:
            The edges connected to the node.
        """
        if key not in self.__nodes:
            raise Exception("Breakpoint node does not exist.")
        return self.__nodes[key]

    def __setitem__(self, key, value):
        """Add a new node to the breakpoint graph.

        Args:
            key: A breakpoint node.
            value: The edges connected to the breakpoint node.
        """
        if type(value) != list:
            raise Exception("Value to be added to the adjacency list must be a list.")

        if key not in self.__nodes:
            self.__nodes[key] = []
        self.__nodes[key] = value
        
    def add_edge(self, edge: Union[BreakpointEdge.BreakpointEdge, SequenceEdge.SequenceEdge], annotation: str):
        """Adds an edge to the appropriate list.
        """
        pass

    def get_breakpoint_graph_dataframe(self) -> pd.DataFrame:
        """Obtain a DataFrame rendering of the breakpoint graph.
        """
        pass

