from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import (
    TYPE_CHECKING,
    Dict,
    Generic,
    NamedTuple,
    Set,
    TypeVar,
)

import intervaltree
import numpy as np
import pyomo.environ as pyo

from coral import types

if TYPE_CHECKING:
    from coral.breakpoint.breakpoint_graph import BreakpointGraph

T = TypeVar("T")


class Strand(str, enum.Enum):
    FORWARD = "+"
    REVERSE = "-"

    @property
    def inverse(self) -> Strand:
        return Strand.FORWARD if self == Strand.REVERSE else Strand.REVERSE


class CNSInterval(NamedTuple):
    chr_tag: str
    start: int
    end: int
    raw_cn: float


@dataclass
class EditDistanceStats:
    """Stores edit distance from reference genome (NM tag) statistics."""

    count: int = 0
    total: int = 0
    sum_of_squares: int = 0

    def observe(self, value: int) -> None:
        self.count += 1
        self.total += value
        self.sum_of_squares += np.pow(value, 2)

    @property
    def mean(self) -> float:
        return self.total / self.count

    @property
    def std_dev(self) -> float:
        return np.sqrt(self.sum_of_squares / self.count - np.pow(self.mean, 2))


class CigarEnds(NamedTuple):
    start: int
    end: int


class CigarAlignment(NamedTuple):
    """
    The start and end position on the read on positive strand, and the alignment
    length on the reference genome.
    """

    start: int
    end: int
    ref_length: int  # Length on reference genome (!= read query length above)


@dataclass
class Interval:
    chr_tag: str
    start: int
    end: int

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __lt__(self, other: Interval) -> bool:
        return (self.chr_tag, self.start, self.end) < (
            other.chr_tag,
            other.start,
            other.end,
        )

    @property
    def left(self) -> int:
        return min(self.start, self.end)

    @property
    def right(self) -> int:
        return max(self.start, self.end)

    def does_overlap(self, other: Interval) -> bool:
        """Check if two chromosome intervals overlap (share a subsequence).

        Intervals are given in the form [chr, start, end], where:
            chr: chromosome number
            s: start position/index
            e: end position/index
        """
        return (
            self.chr_tag == other.chr_tag
            and self.start <= other.end
            and self.end >= other.start
        )

    def is_adjacent(self, other: Interval) -> bool:
        """Check if two intervals are adjacent to each other."""
        if self.chr_tag != other.chr_tag:
            return False
        if self.start <= other.start:
            return other.start == self.end + 1
        return self.start == other.end + 1

    def intersects(self, y: Interval, extend=0, margin=0.0) -> bool:
        if margin > 0.0:
            margin_offset = (1 - margin) * (y.end - y.start)
            margin_interval = (
                (y.start + margin_offset, y.end - margin_offset),
                (self.start + margin_offset, self.end - margin_offset),
            )

            for start, end in margin_interval:
                if (
                    self.chr_tag == y.chr_tag
                    and self.start <= end
                    and self.end >= start
                ):
                    return True
            return False

        # Adjust the intervals with the extension
        self_start, self_end = max(0, self.start - extend), self.end + extend
        n_start, n_end = y.start, y.end

        # Check for chromosome match and interval overlap
        if self.chr_tag != y.chr_tag:
            return False

        return self_start <= n_end and self_end >= n_start

    # Sourced from AA
    def intersection(self, y: Interval) -> Interval | None:
        if not self.intersects(y):
            return None
        return Interval(
            self.chr_tag, max(self.start, y.start), min(self.end, y.end)
        )

    # Sourced from AA
    def merge(self, y: Interval, extend=0) -> Interval | None:
        if not self.intersects(y, extend):
            return None
        return Interval(
            self.chr_tag, min(self.start, y.start), max(self.end, y.end)
        )


@dataclass
class ReadInterval(Interval):
    """Reference genome interval that a read has been matched to."""

    strand: Strand
    name: str


@dataclass
class ChimericAlignment:
    query_ends: CigarEnds
    ref_interval: ReadInterval
    mapq: int
    edit_dist: float

    # Matching CN segment indices, if found
    cns: set[int] = field(default_factory=set)

    def __lt__(self, other: ChimericAlignment) -> bool:
        return self.query_ends < other.query_ends


@dataclass
class AmpliconInterval(Interval):
    amplicon_id: int = -1

    def __str__(self) -> str:
        return f"{self.amplicon_id} {self.chr_tag}:{self.start}-{self.end}"


@dataclass
class IntervalWithReads(Interval):
    reads: set[str] = field(default_factory=set)


class BPReads(NamedTuple):
    """Container for storing reads that support a breakpoint."""

    name: str
    read1: int
    read2: int


class ChrPairOrientation(NamedTuple):
    chr1: str
    chr2: str
    strand1: Strand
    strand2: Strand


@dataclass
class Breakpoint:
    chr1: str
    start: int  # Start position of the breakpoint
    strand1: Strand
    chr2: str
    end: int  # End position of the breakpoint
    strand2: Strand
    read_info: BPReads
    gap: int  # Gap between interval endpoints
    was_reversed: bool
    mapq1: int
    mapq2: int


@dataclass
class EdgeToCN:
    """Container class for mapping edge indices (by type) to CN (copy number) values."""

    sequence: Dict[int, float] = field(default_factory=dict)
    concordant: Dict[int, float] = field(default_factory=dict)
    discordant: Dict[int, float] = field(default_factory=dict)
    source: Dict[int, float] = field(default_factory=dict)

    @staticmethod
    def from_graph(bp_graph: BreakpointGraph):
        return EdgeToCN(
            sequence={
                i: edge[-1] for i, edge in enumerate(bp_graph.sequence_edges)
            },
            concordant={
                i: edge[-1] for i, edge in enumerate(bp_graph.concordant_edges)
            },
            discordant={
                i: edge[-1] for i, edge in enumerate(bp_graph.discordant_edges)
            },
            source={
                i: edge[-1] for i, edge in enumerate(bp_graph.source_edges)
            },
        )


class WalkData(NamedTuple, Generic[T]):
    """Container for storing graph walk data, separated into cycles and paths."""

    cycles: list[T]
    paths: list[T]


@dataclass
class CycleSolution:
    """Container for storing MIQCP Pyomo solution state.

    This consists of parsing solved cycles, their weights, and corresponding
    satisfied path constraints.

    """

    solver_status: pyo.SolverStatus
    termination_condition: pyo.TerminationCondition
    total_weights_included: float = 0.0
    walks: WalkData[types.AmpliconWalk] = field(
        default_factory=lambda: WalkData([], [])
    )
    walk_weights: WalkData[float] = field(
        default_factory=lambda: WalkData([], [])
    )
    satisfied_pc: WalkData[list[int]] = field(
        default_factory=lambda: WalkData([], [])
    )  # Each list contains indices of satisfied path constraints for relevant walk.
    satisfied_pc_set: Set[int] = field(default_factory=set)

    @property
    def num_pc_satisfied(self) -> int:
        return len(self.satisfied_pc_set)

    @property
    def num_cycles(self) -> int:
        return len(self.walks.cycles)

    @property
    def num_paths(self) -> int:
        return len(self.walks.paths)


class InitialSolution(NamedTuple):
    """Container for storing initial solution state (generated by a heuristic method or incomplete solver run), to be passed to another model."""

    walks: WalkData[types.AmpliconWalk]
    walk_weights: WalkData[float]
    satisfied_pc: WalkData[list[int]]


class Solver(str, enum.Enum):
    """Enum for specifying the solver to pass to the MIQCP model used for cycle identification."""

    GUROBI = "gurobi"
    SCIP = "scip"
    # Coin-OR solvers
    BONMIN = "bonmin"
    COUENNE = "couenne"


@dataclass
class SolverOptions:
    num_threads: int = -1
    time_limit_s: int = 7200
    output_dir: str = "models"
    model_prefix: str = "pyomo"
    solver: Solver = Solver.GUROBI


class CNSIntervalTree(intervaltree.IntervalTree):
    def pos2cni(self, chr_tag: str, pos: int) -> set[intervaltree.Interval]:
        return self[chr_tag][pos]

    def get_single_cns_idx(self, chr_tag: str, pos: int) -> int | None:
        """Get the IntervalTree index (CN segment) for a given position."""
        intervals = self.pos2cni(chr_tag, pos)
        if len(intervals) > 1:
            raise ValueError(f"Expected 1 interval, not {intervals=}")
        if not intervals:
            return None  # No matching segment found
        return intervals.pop().data

    def get_cns_ends(self, query_interval: Interval) -> tuple[int, int]:
        """Assumes query interval fits within CN segments."""
        left_seg_idx = self.get_single_cns_idx(
            query_interval.chr_tag, query_interval.left
        )
        right_seg_idx = self.get_single_cns_idx(
            query_interval.chr_tag, query_interval.right
        )
        if not left_seg_idx or not right_seg_idx:
            raise KeyError(
                f"Unable to match CNS ends for {query_interval=}, \
                {left_seg_idx=}, {right_seg_idx=}"
            )
        return (left_seg_idx, right_seg_idx)

    def get_cn_segment_indices(self, read_interval: ReadInterval) -> set[int]:
        seg_idxs = set()
        left_seg_idx = self.get_single_cns_idx(
            read_interval.chr_tag, read_interval.left
        )
        right_seg_idx = self.get_single_cns_idx(
            read_interval.chr_tag, read_interval.right
        )
        if left_seg_idx:
            seg_idxs.add(left_seg_idx)
        if right_seg_idx:
            seg_idxs.add(right_seg_idx)
        return seg_idxs
