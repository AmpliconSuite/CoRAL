"""
Tests that d<N> tokens in written path constraints refer to the N-th
discordant line of the *_graph.txt file.
"""

from __future__ import annotations

import io
import re
from typing import TYPE_CHECKING

import pytest

from coral import core_utils
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.breakpoint_utilities import output_breakpoint_graph_lr
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.datatypes import (
    EdgeId,
    EdgeType,
    FinalizedPathConstraint,
    Node,
    OutputPCOptions,
    PathConstraint,
    Strand,
)
from coral.output.cycle_output import (
    output_amplicon_walks,
    write_path_constraint_to_file,
)

if TYPE_CHECKING:
    import pathlib

TOKEN_RE = re.compile(r"([ecd])(\d+)([+-]):\d+")


def build_graph() -> tuple[BreakpointGraph, PathConstraint]:
    """Four segments on chr2 and two SVs joining different segment pairs, whose
    in-memory order (by the stored higher-coordinate node1, as CoRAL builds
    them) is the reverse of their output order (by the lower-coordinate
    endpoint)."""
    g = BreakpointGraph()
    for start in (100, 200, 300, 400):
        g.add_sequence_edge("chr2", start, start + 99, lr_count=10, cn=5.0)
    # Added first: s0 end -> s2 start (lower endpoint chr2:199)
    g.add_discordant_edge(
        Node("chr2", 300, Strand.REVERSE),
        Node("chr2", 199, Strand.FORWARD),
        read_support=5,
    )
    # Added second: s3 end -> s0 start (lower endpoint chr2:100)
    g.add_discordant_edge(
        Node("chr2", 499, Strand.FORWARD),
        Node("chr2", 100, Strand.REVERSE),
        read_support=5,
    )
    # s3+ -> s0+ -> s2+
    d_s3_s0 = g.node_adjacencies[Node("chr2", 499, Strand.FORWARD)].discordant[0]
    d_s0_s2 = g.node_adjacencies[Node("chr2", 199, Strand.FORWARD)].discordant[0]
    path = [
        EdgeId(EdgeType.SEQUENCE, 3),
        Node("chr2", 499, Strand.FORWARD),
        EdgeId(EdgeType.DISCORDANT, d_s3_s0),
        Node("chr2", 100, Strand.REVERSE),
        EdgeId(EdgeType.SEQUENCE, 0),
        Node("chr2", 199, Strand.FORWARD),
        EdgeId(EdgeType.DISCORDANT, d_s0_s2),
        Node("chr2", 300, Strand.REVERSE),
        EdgeId(EdgeType.SEQUENCE, 2),
    ]
    pc = PathConstraint(path=path, support=3)
    g.path_constraints = [pc]
    g.longest_path_constraints = [
        FinalizedPathConstraint(
            edge_counts=core_utils.path_to_edge_count(path),
            pc_idx=0,
            support=3,
        )
    ]
    return g, pc


def parse_graph_file(text: str) -> tuple[list, list, list[str]]:
    seq, disc, pcs = [], [], []
    for line in text.splitlines():
        fields = line.split("\t")
        if fields[0] == "sequence":
            seq.append({fields[1][:-1], fields[2][:-1]})
        elif fields[0] == "discordant":
            a, b = fields[1].split("->")
            disc.append({a[:-1], b[:-1]})
        elif fields[0] == "path_constraint":
            pcs.append(fields[1])
    return seq, disc, pcs


def bad_discordant_tokens(pc_str: str, seq: list, disc: list) -> list[str]:
    """d<N> tokens whose N-th discordant line does not join the sequence
    edges on either side of it."""
    toks = TOKEN_RE.findall(pc_str)
    bad = []
    for i in range(1, len(toks) - 1):
        edge_type, num, _ = toks[i]
        if edge_type != "d":
            continue
        ends = seq[int(toks[i - 1][1]) - 1] | seq[int(toks[i + 1][1]) - 1]
        if not disc[int(num) - 1] <= ends:
            bad.append(f"d{num}")
    return bad


@pytest.mark.parametrize(
    "pc_option", [OutputPCOptions.LONGEST, OutputPCOptions.ALL]
)
def test_graph_file_pc_indices_match_discordant_lines(
    tmp_path: pathlib.Path, pc_option: OutputPCOptions
) -> None:
    g, _ = build_graph()
    out = tmp_path / "test_amplicon1_graph.txt"
    output_breakpoint_graph_lr(g, str(out), pc_option)
    seq, disc, pcs = parse_graph_file(out.read_text())
    assert len(disc) == 2
    assert len(pcs) == 1
    assert bad_discordant_tokens(pcs[0], seq, disc) == []


def test_cycles_file_pc_matches_graph_file(tmp_path: pathlib.Path) -> None:
    """The cycles file must use the same d<N> numbering as the graph file."""
    g, pc = build_graph()
    out = tmp_path / "test_amplicon1_graph.txt"
    output_breakpoint_graph_lr(g, str(out), OutputPCOptions.LONGEST)
    _, _, graph_pcs = parse_graph_file(out.read_text())

    fp = io.StringIO()
    write_path_constraint_to_file(
        0,
        g.longest_path_constraints[0],
        pc.path,
        fp,
    )
    cycles_pc = fp.getvalue().split("\t")[1]
    assert cycles_pc == graph_pcs[0]


OLD_GRAPH = """\
SequenceEdge: StartPosition, EndPosition, PredictedCN, NumberOfLongReads, Size
sequence\tchr1:1-\tchr1:100+\t2.0\t10.0\t100\t10
sequence\tchr1:101-\tchr1:300+\t4.0\t10.0\t200\t10
sequence\tchr1:301-\tchr1:500+\t6.0\t10.0\t200\t10
sequence\tchr1:501-\tchr1:900+\t4.0\t10.0\t400\t10
sequence\tchr1:901-\tchr1:1000+\t2.0\t10.0\t100\t10
BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads
concordant\tchr1:100+->chr1:101-\t2.0\t5
concordant\tchr1:300+->chr1:301-\t4.0\t5
concordant\tchr1:500+->chr1:501-\t4.0\t5
concordant\tchr1:900+->chr1:901-\t2.0\t5
discordant\tchr1:500+->chr1:301-\t2.0\t5
discordant\tchr1:900+->chr1:101-\t2.0\t5
PathConstraint: Path, Support
path_constraint\te4+:1,d2-:1,e2+:1\t3
"""


def test_cycle_mode_keeps_input_numbering(tmp_path: pathlib.Path) -> None:
    """`coral cycle` reads the graph from file and must write the cycles-file
    path constraints with the input's d<N> numbering, even when the input's
    discordant lines are not in lower-end order."""
    g = parse_breakpoint_graph(io.StringIO(OLD_GRAPH))
    g.sort_edges()  # as infer_amplicon_intervals() does in coral cycle
    g.amplicon_idx = 0
    output_amplicon_walks(
        g, str(tmp_path / "test"), pc_output_option=OutputPCOptions.ALL
    )
    cycles_text = (tmp_path / "test_amplicon1_cycles.txt").read_text()
    cycles_pcs = [
        line.split("\t")[1]
        for line in cycles_text.splitlines()
        if line.startswith("Path constraint")
    ]
    seq, disc, _ = parse_graph_file(OLD_GRAPH)
    assert len(cycles_pcs) == 1
    assert "d2" in cycles_pcs[0]
    assert bad_discordant_tokens(cycles_pcs[0], seq, disc) == []
