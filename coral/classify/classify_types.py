from __future__ import annotations

import enum
from dataclasses import dataclass

import intervaltree

from coral import cycle_decomposition
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.classify import classify_utils


class AmpliconType(enum.StrEnum):
    UNKNOWN = "unknown"
    NONE = "none"
    LINEAR = "linear"
    COMPLEX_NON_CYCLIC = "complex-non-cyclic"
    VIRAL = "viral"
    CYCLIC = "cyclic"  # ecDNA
    BFB = "breakage-fusion-bridge"


@dataclass
class AmpliconProperties:
    rearranged_edges: int
    fb_edges: int
    fb_prop: float
    max_cn: int
    tot_amplified_seq_len: int
    bp_read_supports: list[int]
    bp_cns: list[float]
    seq_cns: list[float]
    seq_coverages: list[float]

    @classmethod
    def from_breakpoint_graph(
        cls,
        bp_graph: BreakpointGraph,
        low_complexity_region_map: dict[str, intervaltree.IntervalTree]
        | None = None,
        min_cn_gain: float = 4.5,
        foldback_dist_threshold: int = 25_000,
    ) -> AmpliconProperties:
        if low_complexity_region_map is None:
            low_complexity_region_map = (
                classify_utils.build_low_complexity_region_map()
            )
        max_cn, tot_amplified_seq_len = (
            classify_utils.get_sequence_edge_properties(
                bp_graph, low_complexity_region_map, min_cn_gain
            )
        )
        rearranged_edges, fb_edges, fb_prop = (
            classify_utils.get_discordant_edge_properties(
                bp_graph,
                low_complexity_region_map,
                foldback_dist_threshold,
            )
        )
        return cls(
            rearranged_edges=rearranged_edges,
            fb_edges=fb_edges,
            fb_prop=fb_prop,
            max_cn=max_cn,
            tot_amplified_seq_len=tot_amplified_seq_len,
            bp_read_supports=[
                edge.lr_count for edge in bp_graph.discordant_edges
            ],
            bp_cns=[edge.cn for edge in bp_graph.discordant_edges],
            seq_cns=[edge.cn for edge in bp_graph.sequence_edges],
            seq_coverages=[edge.lr_nc for edge in bp_graph.sequence_edges],
        )
