"""
Implements a class for BreakpointGraph. Will serve as a container for the 
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file. 
"""
from collections import defaultdict
from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from long_read_aa.breakpoints import BreakpointEdge, SequenceEdge


class BreakpointGraph:
    """A container object for the breakpoint graphs."""

    def __init__(
        self,
        sequence_edges: List[SequenceEdge.SequenceEdge] = [],
        discordant_edges: List[BreakpointEdge.BreakpointEdge] = [],
        concordant_edges: List[BreakpointEdge.BreakpointEdge] = [],
        source_edges: List[BreakpointEdge.BreakpointEdge] = [],
        amplicon_intervals: Optional[List[Tuple[str, int, int]]] = None,
    ):
        # Keep these private so they're protected
        self.__sequence_edges = set(sequence_edges)
        self.__discordant_edges = set(discordant_edges)
        self.__concordant_edges = set(concordant_edges)
        self.__source_edges = set(source_edges)
        self.__nodes = set([])
        self.amplicon_intervals = amplicon_intervals

        # backend for adjacency list.
        self.__node_to_edge = defaultdict(list)

        # mapping of edges to nodes
        self.__edge_to_node = defaultdict(list)

        # update graph
        for edge in sequence_edges:
            self.add_sequence_edge(edge)
        for edge in discordant_edges:
            self.add_discordant_edge(edge)
        for edge in concordant_edges:
            self.add_concordant_edge(edge)
        for edge in source_edges:
            self.add_source_edge(edge)

    @property
    def sequence_edges(self):
        return list(self.__sequence_edges)

    @property
    def discordant_edges(self):
        return list(self.__discordant_edges)

    @property
    def concordant_edges(self):
        return list(self.__concordant_edges)

    @property
    def source_edges(self):
        return list(self.__source_edges)

    @property
    def nodes(self):
        return list(self.__nodes)

    def get_edges(self, key: Tuple[str, int, str]) -> List[Any]:
        """Access the edges incident on a breakpoint node.

        Args:
            key: A breakpoint node.

        Returns:
            The edges connected to the node.
        """
        if type(key) != tuple or len(key) != 3:
            raise Exception(
                "Breakpoint node must be of form (chr, pos, orientation)"
            )

        if key not in self.__nodes:
            raise Exception("Breakpoint node not found in graph.")

        return self.__node_to_edge[key]

    def remove_edge(self, breakpoint_edge: BreakpointEdge.BreakpointEdge):
        """Delete edge from graph.

        Deletes an edge from the graph. Specifically, this function will remove
        a concordant, discordant, or source edge from the appropriate list and
        as needed from the adjacency list. The function will also merge together
        sequence edges should a breakpoint node no longer exist after the edge
        removal.

        Args:
            breakpoint_edge: A breakpoint edge.
        """
        if type(breakpoint_edge) != BreakpointEdge.BreakpointEdge:
            raise Exception("Breakpoint edge must be of type BreakpointEdge.")

        if breakpoint_edge in self.__discordant_edges:
            self.__discordant_edges.remove(breakpoint_edge)
        elif breakpoint_edge in self.__concordant_edges:
            self.__concordant_edges.remove(breakpoint_edge)
        elif breakpoint_edge in self.__source_edges:
            self.__source_edges.remove(breakpoint_edge)
        else:
            raise Exception("Breakpoint edge does not exist.")

        incident_nodes = self.__edge_to_node[breakpoint_edge]
        del self.__edge_to_node[breakpoint_edge]
        
        # update node mappings and merge sequence edges as needed
        for node in incident_nodes:
            self.__node_to_edge[node].remove(breakpoint_edge)

            # if node now has no discordant edge, merge sequence edges as needed
            node_edges = self.__node_to_edge[node]
            if len(node_edges) == 2 and all(
                [
                    edge not in self.__discordant_edges
                    and edge not in self.__source_edges for edge in node_edges
                ]
            ):
                if node_edges[0] in self.__concordant_edges:
                    concordant_edge, sequence_edge1 = (
                        node_edges[0],
                        node_edges[1],
                    )
                else:
                    concordant_edge, sequence_edge1 = (
                        node_edges[1],
                        node_edges[0],
                    )

                # get the node connecting to the concordant edge
                pair_node = (
                    concordant_edge.right_breakpoint
                    if node == concordant_edge.left_breakpoint
                    else concordant_edge.left_breakpoint
                )
                print(node, pair_node)

                # don't do anything if the pair node connected to the concordant
                # edge still has a discordant edge.
                if any(
                    [
                        edge in self.__discordant_edges
                        or edge in self.__source_edges
                        for edge in self.__node_to_edge[pair_node]
                    ]
                ):
                    continue

                # remove concordant edge and update mappings
                self.__concordant_edges.remove(concordant_edge)
                for _node in self.__edge_to_node[concordant_edge]:
                    self.__node_to_edge[_node].remove(concordant_edge)
                del self.__edge_to_node[concordant_edge]

                sequence_edge2 = self.__node_to_edge[pair_node][0]

                # remove old sequence edges and update mappings
                self.__sequence_edges.remove(sequence_edge1)
                self.__sequence_edges.remove(sequence_edge2)
                for _node in self.__edge_to_node[sequence_edge1]:
                    self.__node_to_edge[_node].remove(sequence_edge1)
                for _node in self.__edge_to_node[sequence_edge2]:
                    self.__node_to_edge[_node].remove(sequence_edge2)
                del self.__edge_to_node[sequence_edge1]
                del self.__edge_to_node[sequence_edge2]

                # remove node
                nodes_to_remove = [
                    sequence_edge1.right_breakpoint,
                    sequence_edge2.left_breakpoint,
                ]
                for _node in nodes_to_remove:
                    self.__nodes.remove(_node)
                    for edge in self.__node_to_edge[_node]:
                        self.__edge_to_node[edge].remove(_node)
                    del self.__node_to_edge[_node]

                # create new edge
                chromosome = sequence_edge1.chromosome
                start = sequence_edge1.start
                end = sequence_edge2.end
                
                support = {
                    k: sequence_edge1.support[k] + sequence_edge2.support[k]
                    for k in sequence_edge1.support
                }
                copy_number = {
                    k: np.mean(
                        [
                            sequence_edge1.copy_number[k],
                            sequence_edge2.copy_number[k],
                        ]
                    )
                    for k in sequence_edge1.copy_number
                }
                average_read_length = np.mean(
                    [
                        sequence_edge1.average_read_length,
                        sequence_edge2.average_read_length,
                    ]
                )

                new_sequence_edge = SequenceEdge.SequenceEdge(
                    chromosome,
                    start,
                    end,
                    support,
                    copy_number,
                    average_read_length,
                )
                self.add_sequence_edge(new_sequence_edge)

    def add_node(self, key: Tuple[str, int, str]):
        """Add a new node to the breakpoint graph.

        Args:
            key: A breakpoint node.
        """
        if type(key) != tuple or len(key) != 3:
            raise Exception(
                "Breakpoint node must be of form (chr, pos, orientation)"
            )

        if key in self.__nodes:
            pass

        self.__nodes.add(node)
        self.__node_to_edge[node] = []

    def __contains__(self, key: Tuple[str, int, str]) -> bool:
        """Returns whether or not the breakpoint node is in the graph."""
        if type(key) != tuple or len(key) != 3:
            raise Exception(
                "Breakpoint node must be of form (chr, pos, orientation)"
            )

        return key in self.__nodes

    def add_sequence_edge(self, edge: SequenceEdge.SequenceEdge):
        """Adds a SequenceEdge to the graph."""
        self.__sequence_edges.add(edge)

        self.__update_graph(edge)

    def add_concordant_edge(self, edge: BreakpointEdge.BreakpointEdge):
        """Adds a concordant edge to the graph."""

        self.__concordant_edges.add(edge)

        self.__update_graph(edge)

    def add_discordant_edge(self, edge: BreakpointEdge.BreakpointEdge):
        """Adds a discordant edge to the graph."""
        self.__discordant_edges.add(edge)

        self.__update_graph(edge)

    def add_source_edge(self, edge: BreakpointEdge.BreakpointEdge):
        """Adds a source edge to the graph."""
        self.__source_edges.add(edge)

        self.__update_graph(edge)

    def get_breakpoint_dataframe(self) -> pd.DataFrame:
        """Obtain a DataFrame rendering of the breakpoints."""
        pass

    def __update_graph(
        self,
        edge: Union[SequenceEdge.SequenceEdge, BreakpointEdge.BreakpointEdge],
    ):
        """A helper function to update the graph.

        Synchronizes the node-to-edge and edge-to-node mappings after adding
        an edge.

        Args:
            edge: An edge.
        """
        # add breakpoint termini to node list.
        self.__nodes.add(edge.left_breakpoint)
        self.__nodes.add(edge.right_breakpoint)

        self.__node_to_edge[edge.left_breakpoint].append(edge)
        self.__node_to_edge[edge.right_breakpoint].append(edge)

        self.__edge_to_node[edge].append(edge.left_breakpoint)
        self.__edge_to_node[edge].append(edge.right_breakpoint)
