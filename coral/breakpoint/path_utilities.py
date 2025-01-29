from __future__ import annotations

import copy
import logging

from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.datatypes import (
    BPIndexedAlignmentContainer,
    BPIndexedAlignments,
    ChimericAlignment,
    LargeIndelAlignment,
    ReferenceInterval,
    Strand,
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
    logger.debug(f"Read {rn} covers multiple breakpoints {bp_alignments}.")
    logger.debug(f"Blocks of local alignments: {bp_alignments_split}")

    paths = []
    for ai_block in bp_alignments_split:
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
            [chimera.ref_interval for chimera in chimeras],
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
    indel_alignment: LargeIndelAlignment,
    bp_idx: int,
) -> Walk:
    logger.debug(f"Read {rn} covers a single small del breakpoint.")
    if (
        indel_alignment.read_start < indel_alignment.read_end
        and indel_alignment.next_start >= indel_alignment.curr_end
    ) or (
        indel_alignment.read_start >= indel_alignment.read_end
        and indel_alignment.next_start <= indel_alignment.curr_end
    ):
        logger.debug(
            f"Discarded the read due to inconsistent alignment information: "
            f"{indel_alignment}."
        )
        return []

    strand = (
        Strand.FORWARD
        if indel_alignment.next_start < indel_alignment.curr_end
        else Strand.REVERSE
    )
    ref_intvs = [
        ReferenceInterval(
            indel_alignment.chr_tag,
            indel_alignment.read_start,
            indel_alignment.curr_end,
            strand,
            rn,
        ),
        ReferenceInterval(
            indel_alignment.chr_tag,
            indel_alignment.next_start,
            indel_alignment.read_end,
            strand,
            rn,
        ),
    ]
    ai1, ai2 = (1, 0) if strand == Strand.FORWARD else (0, 1)
    logger.debug(
        f"Alignment intervals on reference = {ref_intvs}; "
        f"mapq = {indel_alignment.mapq}; "
        f"bp = {(ai1, ai2, bp_idx)}"
    )
    path = path_constraints.chimeric_alignment_to_path_i(
        bp_graph,
        ref_intvs,
        BPIndexedAlignments(ai1, ai2, bp_idx),
    )
    logger.debug(f"Resulting subpath = {path}")
    return path


def get_multiple_indel_bp_paths(
    bp_graph: BreakpointGraph,
    rn: str,
    indels: list[LargeIndelAlignment],
    indel_alignments: list[BPIndexedAlignments],
    min_bp_match_cutoff: int,
) -> list[Walk]:
    # TODO: verify logic is preserved
    logger.debug(f"Read {rn} covers multiple small del breakpoints.")
    indel_ref_intvs = sorted(
        [
            ReferenceInterval(
                indel.chr_tag,
                min(indel.read_start, indel.read_end),
                max(indel.read_start, indel.read_end),
                Strand.FORWARD,
                rn,
            )
            for indel in indels
        ]
    )
    if len(indels) <= 1 or len(set(indel_ref_intvs)) > 1:
        logger.debug(
            "\tDiscarded the read due to inconsistent alignment information.",
        )
        return []
    bp_alignments_split: list[list[int]] = [[]]
    last_ai = 0
    for i in range(len(indel_alignments)):
        if i == 0 or indel_alignments[i].alignment1 == last_ai:
            bp_alignments_split[-1].append(i)
        else:
            bp_alignments_split.append([i])
        last_ai = indel_alignments[i].alignment1
    logger.debug(f"Blocks of local alignments: {bp_alignments_split}")

    paths = []
    for bp_block in bp_alignments_split:
        ai_list = [
            (indel_alignments[i].alignment1, indel_alignments[i].alignment2)
            for i in bp_block
        ]
        bp_list = [indel_alignments[i].discordant_idx for i in bp_block]
        if len(set(bp_list)) < len(bp_list):
            logger.debug("\tDiscarded the block due to repeated breakpoints.")
            logger.debug(f"\tBlocks of local alignments: {bp_block}.")
            continue

        path = path_constraints.chimeric_alignment_to_path(
            bp_graph,
            indel_ref_intvs,
            ai_list,
            bp_list,
        )
        logger.debug(
            f"Alignment intervals on reference = {indel_ref_intvs}; "
            f"bp = {indel_alignments[i]}"
        )
        logger.debug(f"Resulting subpath = {path}")
        paths.append(path)

    return paths


def get_multiple_bp_and_indel_bp_paths(
    bp_graph: BreakpointGraph,
    rn: str,
    all_bp_alignments: BPIndexedAlignmentContainer,
    chimeras: list[ChimericAlignment],
    indels: list[LargeIndelAlignment],
    min_bp_match_cutoff: int,
) -> list[Walk]:
    chimeric_ref_intvs = [chimera.ref_interval for chimera in chimeras]
    indel_ref_intvs = [
        ReferenceInterval(
            indel.chr_tag,
            min(indel.read_start, indel.read_end),
            max(indel.read_start, indel.read_end),
            Strand.FORWARD,
            rn,
        )
        for indel in indels
    ]

    logger.debug(
        f"Read {rn} covers multiple breakpoints and small del breakpoints."
    )
    logger.debug(f"Alignment intervals on reference = {chimeric_ref_intvs}")
    logger.debug(
        f"Small del alignment intervals on reference = {indel_ref_intvs}"
    )

    if any(
        chimeras[i + 1].query_bounds.start - chimeras[i].query_bounds.end
        < -min_bp_match_cutoff
        for i in range(len(chimeras) - 1)
    ):
        logger.debug(
            f"Discarded the read due to overlapping local alignments: {chimeras}."
        )
        return []

    should_skip = False
    split_chimeric_ref_idxs: list[int] = []

    for indel_ref_intv in indel_ref_intvs:
        found_split_rint = False
        for ri in range(len(chimeric_ref_intvs)):
            rint = chimeric_ref_intvs[ri]
            if (
                indel_ref_intv.chr == rint.chr
                and min(indel_ref_intv.start, indel_ref_intv.end)
                > min(rint.start, rint.end)
                and max(indel_ref_intv.start, indel_ref_intv.end)
                < max(rint.start, rint.end)
            ):
                found_split_rint = True
                split_chimeric_ref_idxs.append(ri)
                break
        if not found_split_rint:
            should_skip = True
            break
    if should_skip:
        logger.debug(
            "\tDiscarded the read due to inconsistent alignment information.",
        )
        return []

    paths = []
    for i, split_chimeric_ref_idx in enumerate(split_chimeric_ref_idxs):
        chimeric_ref_intv = copy.deepcopy(
            chimeric_ref_intvs[split_chimeric_ref_idx]
        )
        indel_ref_intv = indel_ref_intvs[i]
        chimeric_ref_intvs.insert(
            split_chimeric_ref_idx,
            chimeric_ref_intv,
        )
        next_chimeric_ref_intv = chimeric_ref_intvs[split_chimeric_ref_idx + 1]
        if chimeric_ref_intv.strand == Strand.FORWARD:
            chimeric_ref_intv.end = min(
                indel_ref_intv.start, indel_ref_intv.end
            )
            next_chimeric_ref_intv.start = max(
                indel_ref_intv.start, indel_ref_intv.end
            )
        else:
            chimeric_ref_intv.start = max(
                indel_ref_intv.start, indel_ref_intv.end
            )
            next_chimeric_ref_intv.end = min(
                indel_ref_intv.start, indel_ref_intv.end
            )
        for chimeric_bp_alignment in all_bp_alignments.equal:
            if (
                chimeric_bp_alignment.alignment1 >= split_chimeric_ref_idx
                and chimeric_bp_alignment.alignment2 >= split_chimeric_ref_idx
            ):
                chimeric_bp_alignment.alignment1 += 1
                chimeric_bp_alignment.alignment2 += 1
        for indel_bp_alignment in all_bp_alignments.unequal:
            if indel_bp_alignment.alignment1 == i:
                if chimeric_ref_intv.strand == Strand.FORWARD:
                    all_bp_alignments.equal.append(
                        BPIndexedAlignments(
                            split_chimeric_ref_idx + 1,
                            split_chimeric_ref_idx,
                            indel_bp_alignment.discordant_idx,
                        )
                    )
                else:
                    all_bp_alignments.equal.append(
                        BPIndexedAlignments(
                            split_chimeric_ref_idx,
                            split_chimeric_ref_idx + 1,
                            indel_bp_alignment.discordant_idx,
                        )
                    )
        split_bp_alignment_idxs = [[0]]
        last_ai = max(
            all_bp_alignments.equal[0].alignment1,
            all_bp_alignments.equal[0].alignment2,
        )
        for i in range(1, len(all_bp_alignments.equal)):
            if (
                min(
                    all_bp_alignments.equal[i].alignment1,
                    all_bp_alignments.equal[i].alignment2,
                )
                == last_ai
            ):
                split_bp_alignment_idxs[-1].append(i)
            else:
                split_bp_alignment_idxs.append([i])
            last_ai = max(
                all_bp_alignments.equal[i].alignment1,
                all_bp_alignments.equal[i].alignment2,
            )
        logger.debug(f"Blocks of local alignments: {split_bp_alignment_idxs}")

        for ai_block in split_bp_alignment_idxs:
            ai_list = [
                (
                    all_bp_alignments.equal[i].alignment1,
                    all_bp_alignments.equal[i].alignment2,
                )
                for i in ai_block
            ]
            bp_list = [
                all_bp_alignments.equal[i].discordant_idx for i in ai_block
            ]
            if len(set(bp_list)) < len(bp_list):
                logger.debug(
                    "\tDiscarded the block due to repeated breakpoints.",
                )
                logger.debug(f"\tBlocks of local alignments: {ai_block}")
                continue
            path = path_constraints.chimeric_alignment_to_path(
                bp_graph,
                chimeric_ref_intvs,
                ai_list,
                bp_list,
            )
            logger.debug(f"Resulting subpath = {path}")
            paths.append(path)

    return paths


def get_bp_graph_paths(
    bp_graph: BreakpointGraph,
    rn: str,
    all_bp_alignments: BPIndexedAlignmentContainer,
    chimeras: list[ChimericAlignment],
    indels: list[LargeIndelAlignment],
    min_bp_match_cutoff: int,
) -> list[Walk]:
    bp_alignments = sorted(
        all_bp_alignments.unequal,
        key=lambda bp_alignment: min(
            bp_alignment.alignment1, bp_alignment.alignment2
        ),
    )
    indel_alignments = sorted(
        all_bp_alignments.equal,
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
    elif len(bp_alignments) == 0 and len(indel_alignments) == 1:
        paths = [
            get_single_indel_bp_path(
                bp_graph, rn, indels[0], indel_alignments[0].discordant_idx
            )
        ]
    elif len(bp_alignments) == 0 and len(indel_alignments) > 1:
        paths = get_multiple_indel_bp_paths(
            bp_graph,
            rn,
            indels,
            indel_alignments,
            min_bp_match_cutoff,
        )
    else:
        paths = get_multiple_bp_and_indel_bp_paths(
            bp_graph,
            rn,
            all_bp_alignments,
            chimeras,
            indels,
            min_bp_match_cutoff,
        )

    return paths
