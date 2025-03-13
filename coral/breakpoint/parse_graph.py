from __future__ import annotations

import io
import logging
import pathlib
import re
from typing import NamedTuple

import typer

from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    AmpliconInterval,
    Breakpoint,
    ConcordantEdge,
    DiscordantEdge,
    EdgeId,
    EdgeType,
    FinalizedPathConstraint,
    Node,
    PathConstraint,
    SequenceEdge,
    Strand,
    Walk,
)

logger = logging.getLogger(__name__)


def parse_sequence_edge(line: str) -> SequenceEdge:
    s = line.strip().split("\t")
    return SequenceEdge(
        chr=s[1].split(":")[0],
        start=int(s[1].split(":")[1][:-1]),
        end=int(s[2].split(":")[1][:-1]),
        cn=float(s[3]),
        lr_nc=float(s[4]),
        lr_count=int(float(s[6])),
    )


def parse_concordant_edge(line: str) -> ConcordantEdge:
    s = line.strip().split("\t")
    node1_str = s[1].split("->")[0]
    node2_str = s[1].split("->")[1]
    if s[2] == "":
        s.pop(2)  # Some older simulations had a superfluous \t character

    return ConcordantEdge(
        node1=Node(
            chr=node1_str.split(":")[0],
            pos=int(node1_str.split(":")[1][:-1]),
            strand=Strand(node1_str.split(":")[1][-1]),
        ),
        node2=Node(
            chr=node2_str.split(":")[0],
            pos=int(node2_str.split(":")[1][:-1]),
            strand=Strand(node2_str.split(":")[1][-1]),
        ),
        lr_count=int(float(s[3])),  # Some older graphs outputted float
        cn=float(s[2]),
    )


def parse_discordant_edge(line: str) -> DiscordantEdge:
    s = line.strip().split("\t")
    node1_str = s[1].split("->")[0]
    node2_str = s[1].split("->")[1]
    if s[2] == "":
        s.pop(2)  # Some older simulations had a superfluous \t character

    return DiscordantEdge(
        node1=Node(
            chr=node1_str.split(":")[0],
            pos=int(node1_str.split(":")[1][:-1]),
            strand=Strand(node1_str.split(":")[1][-1]),
        ),
        node2=Node(
            chr=node2_str.split(":")[0],
            pos=int(node2_str.split(":")[1][:-1]),
            strand=Strand(node2_str.split(":")[1][-1]),
        ),
        lr_count=int(float(s[3])),  # Some older graphs outputted float
        cn=float(s[2]),
    )


def parse_path(
    bp_graph: BreakpointGraph, s: str
) -> tuple[Walk, dict[EdgeId, int]]:
    path_pieces = s.split(",")
    path: Walk = []
    edge_counts: dict[EdgeId, int] = {}

    # Each path component (or edge) is of form
    # <edge_type><seq_edge_idx><strand>:<edge_count>

    id_str, count_str = path_pieces[0].split(":")
    edge_id = EdgeId(EdgeType(id_str[0]), int(id_str[1:-1]) - 1)
    edge_counts[edge_id] = int(count_str)
    strand = Strand(id_str[-1])

    seq_edge = bp_graph.sequence_edges[edge_id.idx]
    first_node = (
        seq_edge.start_node if strand == Strand.REVERSE else seq_edge.end_node
    )
    path.extend([edge_id, first_node])

    # Iterate through path_pieces in edge pairs, ignore final edge
    for edge_str in path_pieces[1:-1]:
        id_str, count_str = edge_str.split(":")
        edge_id = EdgeId(EdgeType(id_str[0]), int(id_str[1:-1]) - 1)
        edge_counts[edge_id] = int(count_str)

        if edge_id.type == EdgeType.SEQUENCE:
            seq_edge = bp_graph.sequence_edges[edge_id.idx]
            prev_node = (
                seq_edge.start_node
                if strand == Strand.REVERSE  # Use strand from previous BP edge
                else seq_edge.end_node
            )
            path.append(prev_node)
            path.append(edge_id)
            strand = Strand(id_str[-1])

            next_node = (
                seq_edge.start_node
                if strand == Strand.REVERSE  # Use strand from previous BP edge
                else seq_edge.end_node
            )
            path.append(next_node)
        else:
            path.append(edge_id)

    id_str, count_str = path_pieces[-1].split(":")
    edge_id = EdgeId(EdgeType(id_str[0]), int(id_str[1:-1]) - 1)
    edge_counts[edge_id] = int(count_str)
    return path, edge_counts


def parse_breakpoint_graph(graph_file: io.TextIOWrapper) -> BreakpointGraph:
    bp_graph = BreakpointGraph()
    for line in graph_file:
        s = line.strip().split("\t")
        if s[0] == "sequence":
            seq_edge = parse_sequence_edge(line)
            bp_graph.add_sequence_edge(
                seq_edge.chr,
                seq_edge.start,
                seq_edge.end,
                seq_edge.lr_count,
                seq_edge.lr_nc,
                seq_edge.cn,
            )
            bp_graph.max_cn = max(seq_edge.cn, bp_graph.max_cn)
        elif s[0] == "concordant":
            concordant_edge = parse_concordant_edge(line)
            bp_graph.add_concordant_edge(
                concordant_edge.node1,
                concordant_edge.node2,
                lr_count=concordant_edge.lr_count,
                cn=concordant_edge.cn,
            )
        elif s[0] == "discordant":
            discordant_edge = parse_discordant_edge(line)
            bp_graph.add_discordant_edge(
                discordant_edge.node1,
                discordant_edge.node2,
                read_support=discordant_edge.lr_count,
                cn=discordant_edge.cn,
            )
        elif s[0] == "path_constraint":
            # This requires the bp_graph to be built first,
            # as path_constraints are represented using edge indices
            path, edge_counts = parse_path(bp_graph, s[1])
            bp_graph.path_constraints.append(
                PathConstraint(path=path, support=int(s[2]))
            )
            bp_graph.longest_path_constraints.append(
                FinalizedPathConstraint(
                    edge_counts=edge_counts,
                    pc_idx=len(bp_graph.path_constraints) - 1,
                    support=int(s[2]),
                )
            )
        elif s[0] == "interval":
            intv = AmpliconInterval(s[1], int(s[2]), int(s[3]))
            bp_graph.amplicon_intervals.append(intv)
            bp_graph.add_endnode(Node(intv.chr, intv.start, Strand.REVERSE))
            bp_graph.add_endnode(Node(intv.chr, intv.end, Strand.FORWARD))

    # Match line 397 in breakpoint_graph.py
    # TODO: Understand why this is needed
    bp_graph.max_cn -= 1.0

    if not bp_graph.amplicon_intervals:
        logger.warning(
            "No amplicon intervals found in breakpoint graph. Are you "
            "using a *_graph.txt file generated with CoRAL v2.1.0+?"
        )

    return bp_graph


def get_all_graphs_from_dir(bp_dir: pathlib.Path) -> list[BreakpointGraph]:
    bp_graphs = []
    if not (graph_paths := list(bp_dir.glob("*_graph.txt"))):
        raise ValueError(
            "No valid breakpoint graph files found in provided directory."
        )
    for bp_filepath in graph_paths:
        with bp_filepath.open("r") as f:
            parsed_bp_graph = parse_breakpoint_graph(f)
            if not parsed_bp_graph.amplicon_intervals:
                raise ValueError(
                    "Breakpoint graph file does not contain any amplicon "
                    "intervals, are you sure it was generated by "
                    "CoRAL version 2.1.0 or later?"
                )

            # We 1-index on outputting graph files
            amplicon_idx = (
                int(bp_filepath.name.split("_")[-2].split("amplicon")[1]) - 1
            )
            parsed_bp_graph.amplicon_idx = amplicon_idx
            for interval in parsed_bp_graph.amplicon_intervals:
                interval.amplicon_id = amplicon_idx
            bp_graphs.append(parsed_bp_graph)
    bp_graphs.sort(key=lambda x: x.amplicon_idx)
    return bp_graphs


class ReconstructionPaths(NamedTuple):
    graph_path: pathlib.Path
    cycle_path: pathlib.Path | None = None


def get_all_reconstruction_paths_from_dir(
    reconstruction_dir: pathlib.Path,
) -> list[ReconstructionPaths]:
    reconstruction_paths = []
    for bp_filepath in reconstruction_dir.glob("*_graph.txt"):
        cycle_filename = re.sub(
            r"amplicon(\d+)_graph.txt",
            r"amplicon\1_cycles.txt",
            bp_filepath.name,
        )
        cycle_filepath = reconstruction_dir / cycle_filename
        reconstruction_paths.append(
            ReconstructionPaths(
                graph_path=bp_filepath,
                cycle_path=cycle_filepath if cycle_filepath.exists() else None,
            )
        )
    return reconstruction_paths


def get_all_cycle_paths_from_dir(
    reconstruction_dir: pathlib.Path,
) -> dict[int, pathlib.Path]:
    if not (cycle_paths := list(reconstruction_dir.glob("*_cycle.txt"))):
        return {}
    amplicon_idx_to_cycle_path: dict[int, pathlib.Path] = {}
    for cycle_path in cycle_paths:
        amplicon_idx = (
            int(cycle_path.name.split("_")[-2].split("amplicon")[1]) - 1
        )
        amplicon_idx_to_cycle_path[amplicon_idx] = cycle_path
    return amplicon_idx_to_cycle_path
