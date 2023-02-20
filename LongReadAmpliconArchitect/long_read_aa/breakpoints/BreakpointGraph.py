"""
Implements a class for BreakpointGraph. Will serve as a container for the 
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file. 
"""
from typing import List, Optional, Tuple

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
        source_edges: List[BreakpointEdge.BreakpointEdge] = [],
        amplicon_intervals: Optional[List[str]] = None,

    ):
        
        # Keep these private so they're protected
        self.__sequence_edges = sequence_edges
        self.__discordant_edges = discordant_edges
        self.__concordant_edges = concordant_edges
        self.__source_edges = source_edges
        self.amplicon_intervals = amplicon_intervals

        # backend for adjacency list.
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

    def __getitem__(self, key: Tuple[str, int]) -> list:
        """Access the adjacency list of the breakpoint graph.

        Args:
            key: A breakpoint node.

        Returns:
            The edges connected to the node.
        """
        if type(key) != tuple or len(key) != 2:
            raise Exception("Breakpoint node must be of form (chr, pos)")
        
        if key not in self.__nodes:
            return None

        return self.__nodes[key]

    def __delitem__(self, key: Tuple[str, int]):
        """Delete entry from adjacency list.

        Args:
            key: A breakpoint node.
        """
        if key not in self.__nodes:
            raise Exception("Breakpoint node does not exist.")

        if type(key) != tuple or len(key) != 2:
            raise Exception("Breakpoint node must be of form (chr, pos)")

        del self.__nodes[key]

    def __setitem__(self, key: Tuple[str, int], value: list):
        """Add a new node to the breakpoint graph.

        Args:
            key: A breakpoint node.
            value: The edges connected to the breakpoint node.
        """
        if type(value) != list:
            raise Exception("Value to be added to the adjacency list must be a list.")

        if type(key) != tuple or len(key) != 2:
            raise Exception("Breakpoint node must be of form (chr, pos)")

        if key not in self.__nodes:
            self.__nodes[key] = []
        self.__nodes[key] = value

    def __contains__(self, key: Tuple[str, int]) -> bool:
        """Returns whether or not the breakpoint node is in the graph.
        """
        if type(key) != tuple or len(key) != 2:
            raise Exception("Breakpoint node must be of form (chr, pos)")
        
        return key in self.__nodes.keys()
        
    def add_sequence_edge(self, edge:  SequenceEdge.SequenceEdge):
        """Adds a SequenceEdge to the graph.
        """
        pass

    def add_concordant_edge(self, edge: BreakpointEdge.BreakpointEdge):
        """Adds a concordant edge to the graph.
        """
        pass

    def add_discordant_edge(self, edge: BreakpointEdge.BreakpointEdge):
        """Adds a discordant edge to the graph.
        """
        pass

    def add_source_edge(self, edge: BreakpointEdge.BreakpointEdge):
        """Adds a source edge to the graph.
        """
        pass

    def get_breakpoint_graph_dataframe(self) -> pd.DataFrame:
        """Obtain a DataFrame rendering of the breakpoint graph.
        """
        pass

