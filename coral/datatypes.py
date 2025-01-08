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
    def is_forward(self) -> bool:
        return self == Strand.FORWARD

    @property
    def is_reverse(self) -> bool:
        return self == Strand.REVERSE

    @property
    def inverse(self) -> Strand:
        return Strand.FORWARD if self == Strand.REVERSE else Strand.REVERSE


@dataclass
class BasicStatTracker:
    """Calculates running statistics for a stream of numerical observations."""

    # TODO: fix algo to avoid catastrophic cancellation

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
        if not self.count:
            return 0.0
        return np.sqrt(self.sum_of_squares / self.count - np.pow(self.mean, 2))


class CigarBounds(NamedTuple):
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

    def __str__(self) -> str:
        return f"{self.chr_tag}:{self.start}-{self.end}"

    @property
    def left(self) -> int:
        return min(self.start, self.end)

    @property
    def right(self) -> int:
        return max(self.start, self.end)

    def contains(self, chr_tag: str, pos: int) -> bool:
        if self.chr_tag != chr_tag:
            return False
        return self.start <= pos <= self.end

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
class CNSInterval(Interval):
    cn: float


@dataclass
class ReadInterval(Interval):
    """Reference genome interval that a read has been matched to."""

    strand: Strand
    name: str


@dataclass
class ChimericAlignment:
    query_bounds: CigarBounds
    ref_interval: ReadInterval
    mapq: int
    edit_dist: float

    # Matching CN segment indices, if found
    cns: set[int] = field(default_factory=set)

    def __lt__(self, other: ChimericAlignment) -> bool:
        return self.query_bounds < other.query_bounds


@dataclass
class AmpliconInterval(Interval):
    amplicon_id: int = -1  # TODO: update to None for more obvious behavior

    def __str__(self) -> str:
        return f"{self.amplicon_id} {self.chr_tag}:{self.start}-{self.end}"


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
    all_reads: set[BPReads] = field(default_factory=set)

    def is_close(self, other: Breakpoint) -> bool:
        if self.chr1 != other.chr1 or self.chr2 != other.chr2:
            return False
        if self.strand1 != other.strand1 or self.strand2 != other.strand2:
            return False
        if abs(self.start - other.start) >= 200:
            return False
        if abs(self.end - other.end) >= 200:
            return False
        return True

    def __str__(self) -> str:
        return (
            f"{self.chr1}-{self.start}({self.strand1})___"
            f"{self.chr2}-{self.end}({self.strand2})"
        )


@dataclass
class BreakpointStats:
    bp_distance_cutoff: float
    start: BasicStatTracker = field(default_factory=BasicStatTracker)
    end: BasicStatTracker = field(default_factory=BasicStatTracker)
    mapq1: float = 0.0
    mapq2: float = 0.0

    def observe(self, bp: Breakpoint) -> None:
        self.start.observe(bp.start)
        self.end.observe(bp.end)
        if not bp.was_reversed:
            self.mapq1 += bp.mapq1
            self.mapq2 += bp.mapq2
        else:
            self.mapq1 += bp.mapq2
            self.mapq2 += bp.mapq1

    @property
    def start_window(self) -> float:
        default_window = self.bp_distance_cutoff / 2.99
        if self.start.count == 0:
            return default_window
        return max(self.start.std_dev, default_window)

    @property
    def end_window(self) -> float:
        default_window = self.bp_distance_cutoff / 2.99
        if self.end.count == 0:
            return default_window
        return max(self.end.std_dev, default_window)


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
                i: edge.cn for i, edge in enumerate(bp_graph.discordant_edges)
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
    def get_single_cns_idx(self, pos: int) -> int | None:
        """Get the IntervalTree index (CN segment) for a given position."""
        intervals = self[pos]
        if len(intervals) > 1:
            raise ValueError(f"Expected 1 interval, not {intervals=}")
        if not intervals:
            return None  # No matching segment found
        return intervals.pop().data

    def get_cns_ends(self, query_interval: Interval) -> tuple[int, int]:
        """Assumes query interval fits within CN segments."""
        left_seg_idx = self.get_single_cns_idx(query_interval.left)
        right_seg_idx = self.get_single_cns_idx(query_interval.right)
        if not left_seg_idx or not right_seg_idx:
            raise KeyError(
                f"Unable to match CNS ends for {query_interval=}, \
                {left_seg_idx=}, {right_seg_idx=}"
            )
        return (left_seg_idx, right_seg_idx)

    def get_cn_segment_indices(self, read_interval: ReadInterval) -> set[int]:
        seg_idxs = set()
        left_seg_idx = self.get_single_cns_idx(read_interval.left)
        right_seg_idx = self.get_single_cns_idx(read_interval.right)
        if left_seg_idx:
            seg_idxs.add(left_seg_idx)
        if right_seg_idx:
            seg_idxs.add(right_seg_idx)
        return seg_idxs


@dataclass
class SingleArmInfo:
    """Container for storing chromosome arm-specific CN segment information."""

    interval: Interval
    size: int
    segs: list[CNSInterval] = field(default_factory=list)
    ccn: float = 2.0  # TODO: verify meaning of CCN here + in `aggregate_arm_cn`

    @property
    def total_length(self) -> float:
        return sum(len(seg) for seg in self.segs)


@dataclass
class ChrArmInfo:
    interval: Interval  # Combined interval spanning p + q arms
    p_arm: SingleArmInfo
    q_arm: SingleArmInfo


class BPToCNI(NamedTuple):
    """Container for storing breakpoint-to-CN segment index mappings."""

    cni: types.CNSIdx  # CN segment index within CNS Interval Tree
    pos: int  # Breakpoint position (start/end) that falls within the CN segment
    bp_idx: types.BPIdx  # Breakpoint index


class BPToChrCNI(NamedTuple):
    """Container for storing breakpoint-to-CN segment index mappings, along with
    chromosome specification. Necessary for breakpoints that span multiple
    chromosomes."""

    chr: types.ChrTag  # Chromosome tag
    cni: types.CNSIdx  # CN segment index within CNS Interval Tree
    pos: int  # Breakpoint position (start/end) that falls within the CN segment
    bp_idx: types.BPIdx  # Breakpoint index


class Node(NamedTuple):
    """Container for storing info about a specific genomic point."""

    chr: types.ChrTag
    pos: int
    strand: Strand

    def __str__(self) -> str:
        return f"{self.chr}-{self.pos}({self.strand})"


@dataclass
class DiscordantEdge:
    """Container for storing discordant edge information."""

    chr1: str
    pos1: int
    strand1: Strand
    chr2: str
    pos2: int
    strand2: Strand

    sr_count: int  # Short read count
    sr_flag: str  # Short read flag
    sr_cn: float  # Short read Copy Number

    lr_count: int  # Long read count
    reads: set[BPReads] = field(default_factory=set)
    cn: float = 0.0  # Edge Copy Number

    @property
    def is_self_loop(self) -> bool:
        return (
            self.chr1 == self.chr2
            and self.strand1 == self.strand2
            and self.pos1 == self.pos2
        )

    def matches_bp(self, bp: Breakpoint) -> bool:
        return (
            self.chr1 == bp.chr1
            and self.chr2 == bp.chr2
            and self.strand1 == bp.strand1
            and self.strand2 == bp.strand2
            and self.pos1 == bp.start
            and self.pos2 == bp.end
        )
