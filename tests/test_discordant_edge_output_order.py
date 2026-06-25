"""Tests for discordant-edge directionality in the graph output.

A discordant edge (SV) must be written with its lower-coordinate end first
(earlier chromosome first when inter-chromosomal), matching the
AmpliconArchitect ``_graph.txt`` convention. Internally CoRAL stores
``node1`` as the higher-coordinate end (see ``interval2bp`` / the ``bp_match``
invariant), so the writer relies on ``DiscordantEdge.ordered_nodes`` to flip the
endpoints when needed.

The cases below mirror the three discordant edges from sample GBM39, whose
directionality was being emitted reversed relative to AA before the fix.
"""

from coral.datatypes import DiscordantEdge, Node, Strand


def fmt(node: Node) -> str:
    return f"{node.chr}:{node.pos}{node.strand}"


def edge_str(de: DiscordantEdge) -> str:
    first, second = de.ordered_nodes
    return f"{fmt(first)}->{fmt(second)}"


# ---------------------------------------------------------------------------
# GBM39 regression: stored node1 is the higher coordinate; output must lead
# with the lower coordinate and keep each strand attached to its position.
# ---------------------------------------------------------------------------

def test_intrachromosomal_edges_written_lower_first() -> None:
    cases = [
        # (stored node1 [higher], stored node2 [lower], expected output)
        (
            Node("chr7", 55155021, Strand.REVERSE),
            Node("chr7", 55127266, Strand.FORWARD),
            "chr7:55127266+->chr7:55155021-",
        ),
        (
            Node("chr7", 55610095, Strand.REVERSE),
            Node("chr7", 55609190, Strand.FORWARD),
            "chr7:55609190+->chr7:55610095-",
        ),
        (
            Node("chr7", 56049369, Strand.FORWARD),
            Node("chr7", 54763282, Strand.REVERSE),
            "chr7:54763282-->chr7:56049369+",
        ),
    ]
    for node1, node2, expected in cases:
        de = DiscordantEdge(node1, node2, lr_count=1)
        assert edge_str(de) == expected


def test_already_lower_first_is_unchanged() -> None:
    """An edge already stored lower-first must be emitted verbatim."""
    de = DiscordantEdge(
        Node("chr7", 54763282, Strand.REVERSE),
        Node("chr7", 56049369, Strand.FORWARD),
        lr_count=1,
    )
    assert edge_str(de) == "chr7:54763282-->chr7:56049369+"


def test_interchromosomal_orders_by_chromosome() -> None:
    """Earlier chromosome leads, regardless of position within each chrom."""
    de = DiscordantEdge(
        Node("chr5", 100, Strand.FORWARD),
        Node("chr2", 900, Strand.REVERSE),
        lr_count=1,
    )
    assert edge_str(de) == "chr2:900-->chr5:100+"


def test_foldback_same_position_reverse_first() -> None:
    """At an equal position the reverse-strand end is written first, matching
    AA's ``absPos + 0.4 * strand`` vertex ordering."""
    de = DiscordantEdge(
        Node("chr7", 500, Strand.FORWARD),
        Node("chr7", 500, Strand.REVERSE),
        lr_count=1,
    )
    assert edge_str(de) == "chr7:500-->chr7:500+"
