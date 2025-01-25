from __future__ import annotations

import logging

from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    BPIndexedAlignmentContainer,
    BPIndexedAlignments,
    ChimericAlignment,
    Walk,
)
from coral.models import path_constraints

logger = logging.getLogger(__name__)


def get_single_bp_path(
    bp_graph: BreakpointGraph,
    rn: str,
    chimeras: list[ChimericAlignment],
    bp_alignment: BPIndexedAlignments,
) -> Walk:
    logger.debug(f"Read {rn} covers a single breakpoint.")
    logger.debug(
        f"Alignment intervals on reference = {chimeras}; "
        f"bp = {bp_alignment}"
    )
    path = path_constraints.chimeric_alignment_to_path_i(
        bp_graph,
        [chimera.ref_interval for chimera in chimeras],
        bp_alignment,
    )
    logger.debug(f"Resulting subpath = {path}")
    return path


def get_multiple_bp_paths(
    bp_graph: BreakpointGraph,
    rn: str,
    chimeras: list[ChimericAlignment],
    bp_alignments: list[BPIndexedAlignments],
    min_bp_match_cutoff: int,
) -> list[Walk]:
    if any(
        chimeras[i + 1].query_bounds.start - chimeras[i].query_bounds.end
        < -min_bp_match_cutoff
        for i in range(len(chimeras) - 1)
    ):
        logger.debug(
            f"Discarded the read due to overlapping local alignments: {chimeras}."
        )
        return []

    bp_alignments_split = [[0]]
    last_ai = max(bp_alignments[0].alignment1, bp_alignments[0].alignment2)
    for i in range(1, len(bp_alignments)):
        if (
            min(bp_alignments[i].alignment1, bp_alignments[i].alignment2)
            == last_ai
        ):
            bp_alignments_split[-1].append(i)
        else:
            bp_alignments_split.append([i])
        last_ai = max(bp_alignments[i].alignment1, bp_alignments[i].alignment2)
    logger.debug(f"Read {rn} covers multiple breakpoints.")
    logger.debug(f"Blocks of local alignments: {bp_alignments_split}")

    paths = []
    for ai_block in bp_alignments_split:
        chimeric_intvs = [chimeras[i].ref_interval for i in ai_block]
        ai_list = [
            (bp_alignments[i].alignment1, bp_alignments[i].alignment2)
            for i in ai_block
        ]
        bp_list = [bp_alignments[i].discordant_idx for i in ai_block]
        if len(set(bp_list)) < len(bp_list):
            logger.debug("\tDiscarded the block due to repeated breakpoints.")
            logger.debug(f"\tBlocks of local alignments: {ai_block}.")
            continue

        path = path_constraints.chimeric_alignment_to_path(
            bp_graph,
            chimeric_intvs,
            ai_list,
            bp_list,
        )

        logger.debug(
            f"Alignment intervals on reference = {chimeras}; "
            f"bp = {bp_alignments[i]}"
        )
        logger.debug(f"Resulting subpath = {path}")
        paths.append(path)
    return paths


def get_single_indel_bp_path(
    bp_graph: BreakpointGraph,
    rn: str,
    chimeras: list[ChimericAlignment],
    bp_alignment: BPIndexedAlignments,
) -> Walk:
    return []


def get_bp_graph_paths(
    bp_graph: BreakpointGraph,
    rn: str,
    all_bp_alignments: BPIndexedAlignmentContainer,
    chimeras: list[ChimericAlignment],
    min_bp_match_cutoff: int,
) -> list[Walk]:
    bp_alignments = sorted(
        all_bp_alignments.equal,
        key=lambda bp_alignment: min(
            bp_alignment.alignment1, bp_alignment.alignment2
        ),
    )
    indel_alignments = sorted(
        all_bp_alignments.unequal,
        key=lambda bp_alignment: min(
            bp_alignment.alignment1, bp_alignment.alignment2
        ),
    )

    paths = []

    if len(bp_alignments) == 1 and len(indel_alignments) == 0:
        paths = [get_single_bp_path(bp_graph, rn, chimeras, bp_alignments[0])]
    elif len(bp_alignments) > 1 and len(indel_alignments) == 0:
        paths = get_multiple_bp_paths(
            bp_graph, rn, chimeras, bp_alignments, min_bp_match_cutoff
        )
    return paths
