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
    Union,
)

import intervaltree
import numpy as np
import pyomo.environ as pyo

from coral import types
from coral.constants import CHR_TAG_TO_IDX

if TYPE_CHECKING:
    from coral.breakpoint.breakpoint_graph import BreakpointGraph

T = TypeVar("T")


EdgeIdx = int
EdgeCount = int


# EdgeType = Literal["e", "c", "d", "s", "t", "ns", "nt", "$"]
class EdgeType(enum.StrEnum):
    SEQUENCE = "e"
    CONCORDANT = "c"
    DISCORDANT = "d"
    SOURCE = "s"
    SINK = "t"
    SYNTHETIC_SOURCE = "ns"
    SYNTHETIC_SINK = "nt"
    TERMINAL = "$"

    @property
    def is_synthetic(self) -> bool:
        return self in (EdgeType.SYNTHETIC_SOURCE, EdgeType.SYNTHETIC_SINK)


class EdgeId(NamedTuple):
    type: EdgeType
    idx: EdgeIdx


class DirectedEdge(NamedTuple):
    idx: EdgeIdx
    strand: Strand


OptimizationWalk = dict[EdgeId, EdgeCount]


class Strand(enum.StrEnum):
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
    total: float = 0
    sum_of_squares: float = 0

    def observe(self, value: float) -> None:
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


class Node(NamedTuple):
    """Container for storing info about a specific genomic point."""

    chr: types.ChrTag
    pos: int
    strand: Strand

    def __str__(self) -> str:
        return f"{self.chr}-{self.pos:,d}[{self.strand.value}]"

    def __repr__(self) -> str:
        return self.__str__()


@dataclass
class SplitInterval:
    start: int
    end: int
    strand: Strand


@dataclass
class Interval:
    chr: str
    start: int
    end: int

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __lt__(self, other: Interval) -> bool:
        return (CHR_TAG_TO_IDX[self.chr], self.start, self.end) < (
            CHR_TAG_TO_IDX[other.chr],
            other.start,
            other.end,
        )

    def __str__(self) -> str:
        return f"{self.chr}:{self.start:,d}-{self.end:,d}"

    @property
    def left(self) -> int:
        return min(self.start, self.end)

    @property
    def right(self) -> int:
        return max(self.start, self.end)

    @property
    def reverse(self) -> Interval:
        return Interval(self.chr, self.end, self.start)

    def contains(self, chr_tag: str, pos: int) -> bool:
        if self.chr != chr_tag:
            return False
        return self.start <= pos <= self.end

    def contains_node(self, node: Node) -> bool:
        if self.chr != node.chr:
            return False
        return self.start < node.pos < self.end  # Non-inclusive

    def does_overlap(self, other: Interval) -> bool:
        """Check if two chromosome intervals overlap (share a subsequence).

        Intervals are given in the form [chr, start, end], where:
            chr: chromosome number
            s: start position/index
            e: end position/index
        """
        return (
            self.chr == other.chr
            and self.start <= other.end
            and self.end >= other.start
        )

    def is_adjacent(self, other: Interval) -> bool:
        """Check if two intervals are adjacent to each other."""
        if self.chr != other.chr:
            return False
        if self.start <= other.start:
            return other.start == self.end + 1
        return self.start == other.end + 1

    def intersects(
        self, y: Interval, extend: int = 0, margin: float = 0.0
    ) -> bool:
        if margin > 0.0:
            margin_offset = (1 - margin) * (y.end - y.start)
            margin_interval = (
                (y.start + margin_offset, y.end - margin_offset),
                (self.start + margin_offset, self.end - margin_offset),
            )

            for start, end in margin_interval:
                if (
                    self.chr == y.chr
                    and self.start <= end
                    and self.end >= start
                ):
                    return True
            return False

        # Adjust the intervals with the extension
        self_start, self_end = max(0, self.start - extend), self.end + extend
        n_start, n_end = y.start, y.end

        # Check for chromosome match and interval overlap
        if self.chr != y.chr:
            return False

        return self_start <= n_end and self_end >= n_start

    # Sourced from AA
    def intersection(self, y: Interval) -> Interval | None:
        if not self.intersects(y):
            return None
        return Interval(
            self.chr, max(self.start, y.start), min(self.end, y.end)
        )

    # Sourced from AA
    def merge(self, y: Interval, extend: int = 0) -> Interval | None:
        if not self.intersects(y, extend):
            return None
        return Interval(
            self.chr, min(self.start, y.start), max(self.end, y.end)
        )


@dataclass
class DirectedInterval(Interval):
    strand: Strand


@dataclass
class CNSInterval(Interval):
    cn: float


@dataclass
class ReferenceInterval(Interval):
    """Reference genome interval that a read has been matched to."""

    strand: Strand
    name: str  # Read name


@dataclass
class AmpliconInterval(Interval):
    amplicon_id: int = -1  # TODO: update to None for more obvious behavior

    def __str__(self) -> str:
        return f"Amplicon{self.amplicon_id}>{super().__str__()}"


@dataclass
class ChimericAlignment:
    query_bounds: CigarBounds
    ref_interval: ReferenceInterval
    mapq: int
    edit_dist: float

    # Matching CN segment indices, if found
    cns: set[int] = field(default_factory=set)

    def __lt__(self, other: ChimericAlignment) -> bool:
        return self.query_bounds < other.query_bounds


@dataclass
class LargeIndelAlignment:
    # TODO: add better docstring

    chr_tag: types.ChrTag
    next_start: int
    curr_end: int
    read_start: int
    read_end: int
    mapq: int  # Read mapping quality


class BPAlignments(NamedTuple):
    """Container for storing local alignments from the same read,
    that support a breakpoint."""

    name: str
    alignment1: int
    alignment2: int


class ChrPairOrientation(NamedTuple):
    chr1: str
    chr2: str
    strand1: Strand
    strand2: Strand


@dataclass
class Breakpoint:
    # TODO: refactor as node1, node2 for easier integration with edge types
    node1: Node
    node2: Node

    alignment_info: BPAlignments
    gap: int  # Gap between interval endpoints
    was_reversed: bool
    mapq1: int
    mapq2: int
    all_alignments: set[BPAlignments] = field(default_factory=set)
    amplicon_id: int | None = None
    stats: BreakpointStats | None = None

    @property
    def chr1(self) -> types.ChrTag:
        return self.node1.chr

    @property
    def chr2(self) -> types.ChrTag:
        return self.node2.chr

    @property
    def pos1(self) -> int:
        return self.node1.pos

    @property
    def pos2(self) -> int:
        return self.node2.pos

    @property
    def strand1(self) -> Strand:
        return self.node1.strand

    @property
    def strand2(self) -> Strand:
        return self.node2.strand

    @property
    def read_support(self) -> int:
        return len(self.all_alignments)

    def is_close(self, other: Breakpoint) -> bool:
        if self.chr1 != other.chr1 or self.chr2 != other.chr2:
            return False
        if self.strand1 != other.strand1 or self.strand2 != other.strand2:
            return False
        if abs(self.pos1 - other.pos1) >= 200:
            return False
        return abs(self.pos2 - other.pos2) < 200

    def __str__(self) -> str:
        return f"BP({self.node1}->{self.node2}, {self.alignment_info})"

    def __repr__(self) -> str:
        return (
            f"BP({self.node1}->{self.node2}, {self.alignment_info},"
            f"gap={self.gap}, was_reversed={self.was_reversed},"
            f"mapq1={self.mapq1}, mapq2={self.mapq2},"
            f"read_support={self.read_support})"
        )


@dataclass
class BreakpointStats:
    bp_distance_cutoff: float
    start: BasicStatTracker = field(default_factory=BasicStatTracker)
    end: BasicStatTracker = field(default_factory=BasicStatTracker)
    mapq1: float = 0.0
    mapq2: float = 0.0

    def observe(self, bp: Breakpoint) -> None:
        self.start.observe(bp.pos1)
        self.end.observe(bp.pos2)
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
    def from_graph(bp_graph: BreakpointGraph) -> EdgeToCN:
        return EdgeToCN(
            sequence={
                i: edge.cn for i, edge in enumerate(bp_graph.sequence_edges)
            },
            concordant={
                i: edge.cn for i, edge in enumerate(bp_graph.concordant_edges)
            },
            discordant={
                i: edge.cn for i, edge in enumerate(bp_graph.discordant_edges)
            },
            source={i: edge.cn for i, edge in enumerate(bp_graph.source_edges)},
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
    walks: WalkData[OptimizationWalk] = field(
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

    walks: WalkData[OptimizationWalk]
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
        if left_seg_idx is None or right_seg_idx is None:
            raise KeyError(
                f"Unable to match CNS ends for {query_interval=}, \
                {left_seg_idx=}, {right_seg_idx=}"
            )
        return (left_seg_idx, right_seg_idx)

    def get_cn_segment_indices(
        self, read_interval: ReferenceInterval
    ) -> set[int]:
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


@dataclass
class SourceEdge:
    """Container for storing source edge information."""

    node: Node

    # TODO: previously available, but never used anywhere?
    # sr_count: int  # Short read count
    # sr_flag: str  # Short read flag
    # sr_cn: float  # Short read Copy Number
    # lr_cn: float  # Long read Copy Number
    cn: float = 0.0  # Edge Copy Number


@dataclass
class DiscordantEdge:
    """Container for storing discordant edge information."""

    node1: Node
    node2: Node

    lr_count: int  # Long read count

    alignments: set[BPAlignments] = field(default_factory=set)
    cn: float = 0.0  # Edge Copy Number

    @property
    def is_self_loop(self) -> bool:
        return self.node1 == self.node2

    def matches_bp(self, bp: Breakpoint) -> bool:
        return self.node1 == bp.node1 and self.node2 == bp.node2


@dataclass
class SequenceEdge:
    chr: types.ChrTag
    start: int
    end: int

    lr_nc: float  # Long read Normal Coverage
    lr_count: int  # Long read count

    cn: float = 0.0  # Edge Copy Number

    def __str__(self) -> str:
        return f"SeqEdge({self.chr}-{self.start:,d}-{self.end:,d})"

    def __lt__(self, other: SequenceEdge) -> bool:
        return (CHR_TAG_TO_IDX[self.chr], self.start, self.end) < (
            CHR_TAG_TO_IDX[self.chr],
            other.start,
            other.end,
        )

    def __len__(self) -> int:
        return self.end - self.start + 1

    @property
    def gap(self) -> int:
        return self.end - self.start + 1

    @property
    def start_node(self) -> Node:
        return Node(self.chr, self.start, Strand.REVERSE)

    @property
    def end_node(self) -> Node:
        return Node(self.chr, self.end, Strand.FORWARD)

    @property
    def interval(self) -> Interval:
        return Interval(self.chr, self.start, self.end)


@dataclass
class ConcordantEdge:
    node1: Node
    node2: Node

    lr_count: int  # Long read count
    read_names: set[str] = field(default_factory=set)
    cn: float = 0.0  # Edge Copy Number

    def __lt__(self, other: ConcordantEdge) -> bool:
        return (CHR_TAG_TO_IDX[self.node1.chr], self.node1.pos) < (
            CHR_TAG_TO_IDX[other.node1.chr],
            other.node1.pos,
        )


Edge = Union[SourceEdge, DiscordantEdge, SequenceEdge, ConcordantEdge]


@dataclass
class AdjacencyMatrix:
    """Container for storing adjacency matrix information."""

    sequence: list[int] = field(default_factory=list)
    concordant: list[int] = field(default_factory=list)
    discordant: list[int] = field(default_factory=list)
    source: list[int] = field(default_factory=list)

    @property
    def num_nonsrc_edges(self) -> int:
        return len(self.sequence) + len(self.concordant) + len(self.discordant)

    def get_edges_by_type(self, edge_type: EdgeType) -> list[int]:
        if edge_type == EdgeType.SEQUENCE:
            return self.sequence
        if edge_type == EdgeType.CONCORDANT:
            return self.concordant
        if edge_type == EdgeType.DISCORDANT:
            return self.discordant
        if edge_type == EdgeType.SOURCE:
            return self.source
        if edge_type == EdgeType.SINK:
            return self.source
        raise ValueError(f"Invalid edge type: {edge_type}")


Walk = list[Node | EdgeId]
DirectedWalk = list[DirectedEdge]


@dataclass
class BPIndexedAlignments:
    alignment1: int  # Index of first alignment
    alignment2: int  # Index of second alignment
    discordant_idx: int  # Index of discordant edge


@dataclass
class BPIndexedAlignmentContainer:
    # Alignments whose indices do not match (breakpoints)
    unequal: list[BPIndexedAlignments] = field(default_factory=list)
    # Alignments whose indices match (small del breakpoints)
    equal: list[BPIndexedAlignments] = field(default_factory=list)


@dataclass
class PathConstraint:
    path: Walk  # List of nodes and edges
    support: int = 0  # Number of supporting reads
    amplicon_id: int = -1  # Amplicon index

    def to_file_str(self) -> str:
        """Convert path constraint to string for outputting to file (e.g.,
        *_graph.txt, *_cycle.txt).

            ex: "12+,16+,19+     Support=3"
        """
        path = self.path
        # Reverse path if start node's pos is > than end node's
        if path[0][1] > path[-1][1]:
            path = path[::-1]

        path_str = ""
        for i in range(len(path)):
            if i % 4 == 0:
                edge_id: EdgeId = path[i]  # type: ignore[assignment]
                path_str = f"{edge_id.idx+1}"
                if i < len(path) - 1:
                    next_node: Node = path[i + 1]  # type: ignore[assignment]
                    path_str += f"{next_node.strand.value},"
                else:
                    prev_node: Node = path[i - 1]  # type: ignore[assignment]
                    path_str += f"{prev_node.strand.inverse.value}\t"
        breakpoint()
        path_str += f"Support={self.support}\t"
        return path_str

    def from_file_str(self, line: str) -> PathConstraint:
        """Convert string from file (e.g., *_graph.txt, *_cycle.txt) to
        `PathConstraint` type."""
        pass


@dataclass
class FinalizedPathConstraint:
    edge_counts: dict[EdgeId, int]
    pc_idx: int  # Original PC index (of above `PathConstraint` type)
    support: int


@dataclass
class PathMetric:
    path_idxs: list[int]  # Associated path indices
    path_length: int  # Cumulative path length
    path_support: int  # Cumulative path support
