from __future__ import annotations

import copy
import logging
import os
import pathlib

import typer

from coral.constants import CHR_SIZES, CNSIZE_MAX
from coral.datatypes import ChrArmInfo, CNSInterval, Interval, SingleArmInfo

logger = logging.getLogger(__name__)


def parse_centromere_arms() -> dict[str, ChrArmInfo]:
    chr_arms: dict[str, ChrArmInfo] = {}

    # TODO: clean up structure of reference files
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__))
    )
    with open(
        os.path.join(__location__, "annotations", "GRCh38_centromere.bed")
    ) as fp:
        p_line = fp.readline()
        # Iterate through centromere file in pairs of lines to match p + q arms
        while p_line:
            q_line = fp.readline()
            p_pieces, q_pieces = p_line.strip().split(), q_line.strip().split()
            p_intv = Interval(p_pieces[0], int(p_pieces[1]), int(p_pieces[2]))
            q_intv = Interval(q_pieces[0], int(q_pieces[1]), int(q_pieces[2]))
            if p_intv.chr != q_intv.chr:
                raise ValueError("Centromere file is not sorted by chromosome.")

            full_intv = Interval(p_intv.chr, p_intv.start, q_intv.end)
            chr_arms[p_intv.chr] = ChrArmInfo(
                interval=full_intv,
                p_arm=SingleArmInfo(p_intv, size=p_intv.end),
                q_arm=SingleArmInfo(
                    q_intv, size=CHR_SIZES[p_intv.chr] - q_intv.end
                ),
            )
            p_line = fp.readline()

    return chr_arms


def aggregate_arm_cn(arm: SingleArmInfo) -> float:
    ccn = 2.0
    if arm.total_length < 0.5 * arm.size:
        return ccn

    sum_cns_len = 0
    for cns in sorted(arm.segs, key=lambda cns: cns.cn):
        ccn = cns.cn
        sum_cns_len += len(cns)
        # TODO: shouldn't we use 0.5 for both thresholds?
        if sum_cns_len >= 0.49 * arm.total_length:
            break
    return ccn


def run_seeding(
    cn_seg_file: typer.FileText,
    output_prefix: str,
    gain: float,
    min_seed_size: float,
    max_seg_gap: float,
) -> None:
    """Generate seed intervals from file containing WGS CN calls.

    Breakpoint graph reconstruction initially requires a set of focally
    amplified seed intervals, from which breakpoint edges are explored. This
    method produces these intervals using the given parameters as heuristic
    cutoffs.

    Args:
        cn_seg_file: File containing long-read segmented whole genome CN calls.
        output_prefix: Prefix for output file.
        gain: Minimum CN threshold for an interval to be considered as a seed.
        min_seed_size: Minimum size (in base pairs) of a seed interval.
        max_seg_gap: Maximum gap size (in base pairs) between two adjacent
            potential seed intervals for them to be merged into a single seed.
            If merged, min_seed_size is enforced on the combined interval.

    """

    chr_arms = parse_centromere_arms()
    cnv_seeds: list[list[CNSInterval]] = []
    cur_seed: list[CNSInterval] = []
    for line in cn_seg_file:
        s = line.strip().split()
        if s[0] != "chromosome":
            chr_tag, start, end = s[0], int(s[1]), int(s[2])
            arm_info = chr_arms[chr_tag]
            arm_intv = arm_info.interval
            if cn_seg_file.name.endswith(".cns"):
                cn = 2 * (2 ** float(s[4]))
            elif cn_seg_file.name.endswith(".bed"):
                cn = float(s[3])
            else:
                logger.error(cn_seg_file.name + "\n")
                raise SystemExit("Invalid cn_seg file format!\n")
            # Require absolute CN >= max(gain, cn_cutoff_chrarm)
            cn_intv = CNSInterval(chr_tag, start, end, cn)
            if cn >= gain and (end <= arm_intv.start or start >= arm_intv.end):  # type: ignore[possibly-undefined]
                if (
                    len(cur_seed) > 0
                    and chr_tag == cur_seed[-1].chr
                    and start - cur_seed[-1].end <= max_seg_gap
                ):
                    cur_seed.append(cn_intv)
                elif len(cur_seed) == 0:
                    cur_seed = [cn_intv]
                else:
                    cnv_seeds.append(cur_seed)
                    cur_seed = [cn_intv]
            if end <= arm_intv.start:
                arm_info.p_arm.segs.append(cn_intv)
            if start >= arm_intv.end:
                arm_info.q_arm.segs.append(cn_intv)

    # Add final seed if non-empty
    if cur_seed:
        cnv_seeds.append(cur_seed)

    for arm_info in chr_arms.values():
        arm_info.p_arm.ccn = aggregate_arm_cn(arm_info.p_arm)
        arm_info.q_arm.ccn = aggregate_arm_cn(arm_info.q_arm)

    if output_prefix:
        output_filename = f"{output_prefix}_CNV_SEEDS.bed"
    else:
        output_filename = cn_seg_file.name.replace(".cns", "CNV_SEEDS.bed")

    with open(output_filename, "w") as fp:
        for seed_intvs in cnv_seeds:
            sum_seed_len = sum([len(cns) for cns in seed_intvs])
            cn_cutoff_chrarm = gain
            seed_chr_tag = seed_intvs[-1].chr
            arm_info = chr_arms[seed_chr_tag]
            if sum_seed_len > CNSIZE_MAX:
                cn_cutoff_chrarm = 1.2 * gain
            if seed_intvs[-1].end <= arm_info.p_arm.interval.start:  # p arm
                cn_cutoff_chrarm = cn_cutoff_chrarm + (arm_info.p_arm.ccn - 2.0)
            elif seed_intvs[0].start >= arm_info.q_arm.interval.start:  # q arm
                cn_cutoff_chrarm = cn_cutoff_chrarm + (arm_info.q_arm.ccn - 2.0)
            else:
                os.abort()
            for ci in range(len(seed_intvs))[::-1]:
                if seed_intvs[ci].cn < cn_cutoff_chrarm:
                    del seed_intvs[ci]
            if len(seed_intvs) > 0:
                lastseg: Interval | None = None
                sum_seed_len = 0
                for cns in seed_intvs:
                    if lastseg and cns.start - lastseg.end <= max_seg_gap:
                        sum_seed_len += cns.end - cns.start
                        lastseg.end = cns.end
                    elif lastseg is None:
                        lastseg = copy.deepcopy(cns)
                        sum_seed_len += len(cns)
                    elif sum_seed_len >= min_seed_size:
                        fp.write(
                            f"{lastseg.chr}\t{lastseg.start}\t{lastseg.end-1}\n"
                        )
                        sum_seed_len = 0
                        lastseg = cns
                if sum_seed_len >= min_seed_size:
                    fp.write(
                        f"{lastseg.chr}\t{lastseg.start}\t{lastseg.end-1}\n"  # type: ignore
                    )

    if not cnv_seeds:
        print(f"No seed intervals found with CN>={gain}.")
    print("Created " + output_filename)
