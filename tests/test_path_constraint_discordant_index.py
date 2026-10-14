"""Tests that d<N> tokens in written path constraints refer to the N-th
discordant line of the *_graph.txt file.

The graph writer sorts discordant edges by their lower-coordinate endpoint
(see ``BreakpointEdge.ordered_nodes``), which can differ from the in-memory
order that path constraints index into. Before the fix, path constraints
were written with in-memory indices, so they named the wrong SV whenever the
two orders differed (e.g. CHP212 / COLO320DM ONT reconstructions).
"""

from __future__ import annotations

import io
import re
from typing import TYPE_CHECKING

import pytest

from coral import core_utils
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.breakpoint_utilities import output_breakpoint_graph_lr
from coral.datatypes import (
    EdgeId,
    EdgeType,
    FinalizedPathConstraint,
    Node,
    OutputPCOptions,
    PathConstraint,
    Strand,
)
from coral.output.cycle_output import write_path_constraint_to_file

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
    # In-memory index 0: s0 end -> s2 start (lower endpoint chr2:199)
    g.add_discordant_edge(
        Node("chr2", 300, Strand.REVERSE),
        Node("chr2", 199, Strand.FORWARD),
        read_support=5,
    )
    # In-memory index 1: s3 end -> s0 start (lower endpoint chr2:100)
    g.add_discordant_edge(
        Node("chr2", 499, Strand.FORWARD),
        Node("chr2", 100, Strand.REVERSE),
        read_support=5,
    )
    # s3+ -> s0+ -> s2+
    path = [
        EdgeId(EdgeType.SEQUENCE, 3),
        Node("chr2", 499, Strand.FORWARD),
        EdgeId(EdgeType.DISCORDANT, 1),
        Node("chr2", 100, Strand.REVERSE),
        EdgeId(EdgeType.SEQUENCE, 0),
        Node("chr2", 199, Strand.FORWARD),
        EdgeId(EdgeType.DISCORDANT, 0),
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


def test_output_order_differs_from_memory_order() -> None:
    """Sanity check that the fixture exercises the reordering."""
    g, _ = build_graph()
    assert g.discordant_output_order() == [1, 0]
    assert g.discordant_output_index() == {1: 0, 0: 1}


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
        g.discordant_output_index(),
    )
    cycles_pc = fp.getvalue().split("\t")[1]
    assert cycles_pc == graph_pcs[0]
