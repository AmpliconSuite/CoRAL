import typer

from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    Breakpoint,
    DiscordantEdge,
    EdgeId,
    EdgeType,
    Node,
    PathConstraint,
    SequenceEdge,
    Strand,
    Walk,
)


def parse_sequence_edge(line: str) -> SequenceEdge:
    s = line.strip().split("\t")
    return SequenceEdge(
        chr=s[1].split(":")[0],
        start=int(s[1].split(":")[1][:-1]),
        end=int(s[2].split(":")[1][:-1]),
        lr_nc=float(s[3]),
        lr_count=int(s[6]),
        cn=float(s[5]),
    )


def parse_discordant_edge(line: str) -> DiscordantEdge:
    s = line.strip().split("\t")
    node1_str = s[1].split("->")[0]
    node2_str = s[1].split("->")[1]
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
        lr_count=int(s[3]),
        cn=float(s[2]),
    )


def parse_path(bp_graph: BreakpointGraph, s: str) -> Walk:
    path_pieces = s.split(",")
    path: Walk = []
    # Of form <seq_edge_idx><strand>
    for i, directed_edge in enumerate(path_pieces):
        seq_edge_idx = int(directed_edge[:-1])
        strand = Strand(directed_edge[-1])
        seq_edge = bp_graph.sequence_edges[seq_edge_idx]
        path.append(EdgeId(EdgeType.SEQUENCE, seq_edge_idx))

        if i < 
        node = bp_graph.get_node(seq_edge.end, strand)

        path.append(
            Node(
                chr=s[i + 1].split(":")[0],
                pos=int(s[i + 1].split(":")[1][:-1]),
                strand=Strand(s[i + 1].split(":")[1][-1]),
            )
        )
    return Walk(path)


def parse_path_constraint(
    bp_graph: BreakpointGraph, line: str
) -> PathConstraint:
    s = line.strip().split("\t")
    return PathConstraint(
        path=parse_path(bp_graph, s[1]),
        support=int(s[2]),
        amplicon_id=int(s[3]),
    )


def parse_breakpoint_graph(graph_file: typer.FileText) -> BreakpointGraph:
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
            # as path_constraints are represented using sequence edge indices
            path_constraint = parse_path_constraint(bp_graph, line)
            bp_graph.path_constraints.append(path_constraint)
    return bp_graph
