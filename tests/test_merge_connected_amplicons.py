"""Tests for LongReadBamToBreakpointMetadata.merge_connected_amplicons.

All tests operate purely on in-memory data structures; no BAM file is needed.
"""

from unittest.mock import MagicMock

import pytest

from coral.breakpoint.infer_breakpoint_graph import LongReadBamToBreakpointMetadata
from coral.core_types import Strand
from coral.datatypes import AmpliconInterval, BPAlignments, Breakpoint, Node


def make_metadata(
    amplicon_intervals: list[AmpliconInterval],
    new_bp_list: list[Breakpoint],
) -> LongReadBamToBreakpointMetadata:
    """Construct a minimal metadata object with the given intervals and BPs."""
    obj = LongReadBamToBreakpointMetadata(
        lr_bamfh=MagicMock(),
        bam=MagicMock(),
    )
    obj.amplicon_intervals = amplicon_intervals
    obj.new_bp_list = new_bp_list
    return obj


def make_bp(chr1: str, pos1: int, chr2: str, pos2: int) -> Breakpoint:
    dummy_aln = BPAlignments(name="read1", alignment1=0, alignment2=1)
    return Breakpoint(
        node1=Node(chr=chr1, pos=pos1, strand=Strand.FORWARD),
        node2=Node(chr=chr2, pos=pos2, strand=Strand.REVERSE),
        alignment_info=dummy_aln,
        gap=0,
        was_reversed=False,
        mapq1=60,
        mapq2=60,
    )


# ---------------------------------------------------------------------------
# Basic merge: reproduces the crash case from SOG134_Pre_TM
# ---------------------------------------------------------------------------

def test_cross_amplicon_bp_merges_amplicons() -> None:
    """Two intervals on chr2 with different amplicon_ids, connected by a BP,
    should be unified into a single amplicon."""
    intervals = [
        AmpliconInterval(chr="chr2", start=185_370_001, end=186_410_682, amplicon_id=0),
        AmpliconInterval(chr="chr2", start=188_815_619, end=191_281_309, amplicon_id=1),
    ]
    # BP mirrors the actual crash: chr2:188.9M → chr2:186.3M
    bp = make_bp("chr2", 188_915_618, "chr2", 186_316_524)
    meta = make_metadata(intervals, [bp])

    meta.merge_connected_amplicons()

    ids = {intv.amplicon_id for intv in meta.amplicon_intervals}
    assert len(ids) == 1, "Both intervals should share one amplicon_id after merge"


# ---------------------------------------------------------------------------
# No merge when no connecting BP exists
# ---------------------------------------------------------------------------

def test_no_bp_keeps_amplicons_separate() -> None:
    """Intervals with no connecting BP should retain distinct amplicon_ids."""
    intervals = [
        AmpliconInterval(chr="chr2", start=100_000, end=200_000, amplicon_id=0),
        AmpliconInterval(chr="chr5", start=300_000, end=400_000, amplicon_id=1),
    ]
    meta = make_metadata(intervals, [])

    meta.merge_connected_amplicons()

    assert meta.amplicon_intervals[0].amplicon_id == 0
    assert meta.amplicon_intervals[1].amplicon_id == 1


# ---------------------------------------------------------------------------
# BP within the same amplicon — no change
# ---------------------------------------------------------------------------

def test_intra_amplicon_bp_unchanged() -> None:
    """A BP whose both endpoints fall in the same amplicon should not
    change any amplicon_id."""
    intervals = [
        AmpliconInterval(chr="chr19", start=29_600_001, end=29_908_701, amplicon_id=2),
        AmpliconInterval(chr="chr5",  start=500_000,    end=600_000,    amplicon_id=3),
    ]
    # Both ends of the BP land in the chr19 interval (amplicon_id=2)
    bp = make_bp("chr19", 29_724_664, "chr19", 29_704_462)
    meta = make_metadata(intervals, [bp])

    meta.merge_connected_amplicons()

    assert meta.amplicon_intervals[0].amplicon_id == 2
    assert meta.amplicon_intervals[1].amplicon_id == 3


# ---------------------------------------------------------------------------
# Chain merge: A→B and B→C should all collapse into one amplicon
# ---------------------------------------------------------------------------

def test_chain_merge() -> None:
    """Three amplicons connected A→B and B→C should all be merged."""
    intervals = [
        AmpliconInterval(chr="chr1", start=100_000, end=200_000, amplicon_id=0),
        AmpliconInterval(chr="chr1", start=300_000, end=400_000, amplicon_id=1),
        AmpliconInterval(chr="chr1", start=500_000, end=600_000, amplicon_id=2),
    ]
    bp_ab = make_bp("chr1", 150_000, "chr1", 350_000)
    bp_bc = make_bp("chr1", 350_001, "chr1", 550_000)
    meta = make_metadata(intervals, [bp_ab, bp_bc])

    meta.merge_connected_amplicons()

    ids = {intv.amplicon_id for intv in meta.amplicon_intervals}
    assert len(ids) == 1, "All three intervals should collapse to one amplicon_id"


# ---------------------------------------------------------------------------
# BP endpoint outside all intervals — should skip gracefully, no crash
# ---------------------------------------------------------------------------

def test_bp_outside_intervals_is_skipped() -> None:
    """A BP with one endpoint outside all known intervals should be skipped
    without raising an exception."""
    intervals = [
        AmpliconInterval(chr="chr3", start=100_000, end=200_000, amplicon_id=0),
    ]
    # pos2 is far outside any interval
    bp = make_bp("chr3", 150_000, "chr3", 999_999_999)
    meta = make_metadata(intervals, [bp])

    meta.merge_connected_amplicons()  # must not raise

    assert meta.amplicon_intervals[0].amplicon_id == 0
