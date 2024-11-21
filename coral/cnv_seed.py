from __future__ import annotations

import logging
import os

import typer

from coral.constants import CHR_SIZES, CNSIZE_MAX

logger = logging.getLogger(__name__)


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
    blocked_intervals = []
    chr_arms: dict[str, list[list]] = {}

    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__))
    )
    with open(
        os.path.join(__location__, "annotations", "GRCh38_centromere.bed")
    ) as fp:
        for line in fp:
            s = line.strip().split()
            if len(s[0]) <= 5:  # chr1 - chrM
                blocked_intervals.append([s[0], int(s[1]), int(s[2])])
                if "p" in s[3]:
                    chr_arms[s[0]] = [[int(s[1])], [[], []], [int(s[1])]]
                if "q" in s[3]:
                    chr_arms[s[0]][0].append(int(s[2]))
                    chr_arms[s[0]][2].append(CHR_SIZES[s[0]] - int(s[2]))

    cnv_seeds = []
    cur_seed: list[tuple] = []
    for line in cn_seg_file:
        s = line.strip().split()
        if s[0] != "chromosome":
            if cn_seg_file.name.endswith(".cns"):
                cn = 2 * (2 ** float(s[4]))
            elif cn_seg_file.name.endswith(".bed"):
                cn = float(s[3])
            else:
                logger.error(cn_seg_file.name + "\n")
                logger.error("Invalid cn_seg file format!\n")
            # Require absolute CN >= max(gain, cn_cutoff_chrarm)
            if cn >= gain and (
                int(s[2]) <= chr_arms[s[0]][0][0]
                or int(s[1]) >= chr_arms[s[0]][0][1]
            ):  # type: ignore[possibly-undefined]
                # assume input CN segments sorted by chr and pos
                if (
                    len(cur_seed) > 0
                    and s[0] == cur_seed[-1][0]
                    and int(s[1]) - cur_seed[-1][2] <= max_seg_gap
                ):
                    cur_seed.append((s[0], int(s[1]), int(s[2]), cn))
                elif len(cur_seed) == 0:
                    cur_seed = [(s[0], int(s[1]), int(s[2]), cn)]
                else:
                    cnv_seeds.append(cur_seed)
                    cur_seed = [(s[0], int(s[1]), int(s[2]), cn)]
            if int(s[2]) <= chr_arms[s[0]][0][0]:
                chr_arms[s[0]][1][0].append((s[0], int(s[1]), int(s[2]), cn))
            if int(s[1]) >= chr_arms[s[0]][0][1]:
                chr_arms[s[0]][1][1].append((s[0], int(s[1]), int(s[2]), cn))
    cnv_seeds.append(cur_seed)

    for chr in chr_arms:
        sum_cns_len_parm = sum([cns[2] - cns[1] for cns in chr_arms[chr][1][0]])
        sum_cns_len_qarm = sum([cns[2] - cns[1] for cns in chr_arms[chr][1][1]])
        ccn_p, ccn_q = 2.0, 2.0
        if sum_cns_len_parm >= 0.5 * chr_arms[chr][2][0]:
            scns = sorted(chr_arms[chr][1][0], key=lambda cns: cns[3])
            sum_cns_len_ = 0
            for cns in scns:
                ccn_p = cns[3]
                sum_cns_len_ += cns[2] - cns[1]
                if sum_cns_len_ >= 0.49 * sum_cns_len_parm:
                    break
        if sum_cns_len_qarm >= 0.5 * chr_arms[chr][2][1]:
            scns = sorted(chr_arms[chr][1][1], key=lambda cns: cns[3])
            sum_cns_len_ = 0
            for cns in scns:
                ccn_q = cns[3]
                sum_cns_len_ += cns[2] - cns[1]
                if sum_cns_len_ >= 0.49 * sum_cns_len_qarm:
                    break
        chr_arms[chr].append([ccn_p, ccn_q])

    if output_prefix:
        output_filename = f"{output_prefix}_CNV_SEEDS.bed"
    else:
        output_filename = cn_seg_file.name.replace(".cns", "CNV_SEEDS.bed")
    with open(output_filename, "w") as fp:
        for seed in cnv_seeds:
            sum_seed_len = sum([cns[2] - cns[1] for cns in seed])
            cn_cutoff_chrarm = gain
            if sum_seed_len > CNSIZE_MAX:
                cn_cutoff_chrarm = 1.2 * gain
            if seed[-1][2] <= chr_arms[seed[-1][0]][0][0]:  # p arm
                cn_cutoff_chrarm = (
                    cn_cutoff_chrarm + chr_arms[seed[-1][0]][-1][0] - 2.0
                )
            elif seed[0][1] >= chr_arms[seed[-1][0]][0][1]:  # q arm
                cn_cutoff_chrarm = (
                    cn_cutoff_chrarm + chr_arms[seed[-1][0]][-1][1] - 2.0
                )
            else:
                os.abort()
            for ci in range(len(seed))[::-1]:
                if seed[ci][3] < cn_cutoff_chrarm:
                    del seed[ci]
            if len(seed) > 0:
                lastseg: list = []
                sum_seed_len = 0
                for cns in seed:
                    if len(lastseg) > 0 and cns[1] - lastseg[2] <= max_seg_gap:
                        sum_seed_len += cns[2] - cns[1]
                        lastseg[2] = cns[2]
                    elif len(lastseg) == 0:
                        lastseg = list(cns)
                        sum_seed_len += cns[2] - cns[1]
                    elif sum_seed_len >= min_seed_size:
                        fp.write(
                            "%s\t%d\t%d\n"
                            % (lastseg[0], lastseg[1], lastseg[2] - 1)
                        )
                        sum_seed_len = 0
                        lastseg = list(cns)
                if sum_seed_len >= min_seed_size:
                    fp.write(
                        "%s\t%d\t%d\n"
                        % (lastseg[0], lastseg[1], lastseg[2] - 1)
                    )

    print("Created " + output_filename)
