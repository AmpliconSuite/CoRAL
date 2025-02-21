from __future__ import annotations

import logging
import pathlib
import sys
from typing import Any

import matplotlib as mpl
import pysam
import typer

from coral.breakpoint.breakpoint_utilities import (
    bpc2bp,
    cluster_bp_list,
    interval2bp,
    interval_include,
    interval_overlap,
    interval_overlap_l,
)
from coral.constants import CHR_SIZES
from coral.datatypes import Interval

mpl.use("Agg")
import matplotlib.pyplot as plt
from pylab import rcParams  # type: ignore[import-untyped]

rcParams["figure.figsize"] = [20, 8]
rcParams["pdf.fonttype"] = 42
mpl.rc("xtick", labelsize=25)
mpl.rc("ytick", labelsize=25)

from coral import cigar_parsing, cycle2bed

logger = logging.getLogger(__name__)


def fetch(lr_bamfh):
    read_length = dict()
    chimeric_alignments = dict()
    for read in lr_bamfh.fetch():
        rn = read.query_name
        if read.flag < 256 and rn not in read_length:
            read_length[rn] = read.query_length
        try:
            sa_list = read.get_tag("SA:Z:")[:-1].split(";")
            for sa in sa_list:
                try:
                    if sa not in chimeric_alignments[rn]:
                        chimeric_alignments[rn].append(sa)
                except:
                    chimeric_alignments[rn] = [sa]
        except:
            pass
    reads_wo_primary_alignment = []
    for r in chimeric_alignments:
        if r not in read_length:
            reads_wo_primary_alignment.append(r)
            continue
        rl = read_length[r]
        chimeric_alignments[r] = cigar_parsing.alignment_from_satags(
            chimeric_alignments[r], rl
        )
    for r in reads_wo_primary_alignment:
        del chimeric_alignments[r]
    return read_length, chimeric_alignments


def locate_hsrs(
    lr_bam: pathlib.Path,
    cycle_file: typer.FileText,
    cn_seg_file: typer.FileText,
    output_prefix: str,
    normal_cov: float,
    bp_match_cutoff: int,
    bp_match_cutoff_clustering: int,
):
    ecdna_intervals: list[Interval] = []
    ecdna_intervals_ext: list[Interval] = []

    cycle_filename = cycle_file.name
    if cycle_file.name.endswith("_cycles.txt"):
        # convert it to a bed
        init_char = "" if output_prefix.endswith("/") else "_"
        conv_cycle_fn = output_prefix + init_char + "converted_"
        conv_cycle_fn += "cycles.bed"
        cycle2bed.convert_cycles_to_bed(cycle_file, conv_cycle_fn)
        cycle_filename = conv_cycle_fn

    elif not cycle_filename.endswith(".bed"):
        logger.error(cycle_filename + "\n")
        logger.error(
            "Cycles file must be either a valid *_cycles.txt file or a converted .bed file!\n"
        )
        sys.exit(1)

    with open(cycle_filename) as fp:
        for line in fp:
            if line.startswith("#"):
                continue

            s = line.strip().split()
            ecdna_intervals.append(Interval(s[0], int(s[1]), int(s[2])))
            ecdna_intervals_ext.append(
                Interval(
                    s[0],
                    int(s[1]) - bp_match_cutoff,
                    int(s[2]) + bp_match_cutoff,
                )
            )
    print("ecDNA intervals:")
    for ival in ecdna_intervals:
        print(ival)

    cns_dict: dict = {}
    for line in cn_seg_file:
        s = line.strip().split()
        if line.startswith("chromosome"):
            continue

        if s[0] not in cns_dict:
            cns_dict[s[0]] = {}

        if cn_seg_file.name.endswith(".cns"):
            cn = 2 * (2 ** float(s[4]))
        elif cn_seg_file.name.endswith(".bed"):
            cn = float(s[3])
        else:
            logger.error(cn_seg_file.name + "\n")
            logger.error("Invalid cn_seg file format!\n")

        try:
            cns_dict[s[0]].append([int(s[1]), int(s[2]), cn])  # type: ignore[possibly-undefined]
        except:
            cns_dict[s[0]] = [[int(s[1]), int(s[2]), cn]]  # type: ignore[possibly-undefined]

    lr_bamfh = pysam.AlignmentFile(str(lr_bam), "rb")
    read_length, chimeric_alignments = fetch(lr_bamfh)
    print("Fetched %d chimeric alignments." % len(chimeric_alignments))

    bp_list = []
    for r in chimeric_alignments:
        cycle_flag = False
        r_int = chimeric_alignments[r][0]
        rr_int = chimeric_alignments[r][1]
        q_ = chimeric_alignments[r][2]
        for interval in ecdna_intervals:
            i = interval_overlap_l(interval, rr_int)
            if i >= 0 and interval_include(rr_int[i], interval):
                cycle_flag = True
                break
        if cycle_flag:
            bassigned = [0 for i in range(len(rr_int) - 1)]
            """
			Breakpoint from local alignment i and i + 1
			"""
            for ri in range(len(rr_int) - 1):
                if (
                    q_[ri] >= 20
                    and q_[ri + 1] >= 20
                    and interval_overlap_l(rr_int[ri], ecdna_intervals) == -1
                    and interval_overlap_l(rr_int[ri + 1], ecdna_intervals) >= 0
                ) or (
                    q_[ri] >= 20
                    and q_[ri + 1] >= 20
                    and interval_overlap_l(rr_int[ri], ecdna_intervals) >= 0
                    and interval_overlap_l(rr_int[ri + 1], ecdna_intervals)
                    == -1
                ):
                    bp_list.append(
                        interval2bp(
                            rr_int[ri],
                            rr_int[ri + 1],
                            (r, ri, ri + 1),
                            int(r_int[ri + 1][0]) - int(r_int[ri][1]),
                        )
                        + [q_[ri], q_[ri + 1]]
                    )
                    bassigned[ri] = 1
            """
			Breakpoint from local alignment i - 1 and i + 1
			"""
            for ri in range(1, len(rr_int) - 1):
                if (
                    bassigned[ri - 1] == 0
                    and bassigned[ri] == 0
                    and q_[ri] < 10
                    and q_[ri - 1] >= 20
                    and q_[ri + 1] >= 20
                    and interval_overlap_l(rr_int[ri - 1], ecdna_intervals)
                    == -1
                    and interval_overlap(rr_int[ri + 1], ecdna_intervals) >= 0  # type: ignore[arg-type]
                ) or (
                    bassigned[ri - 1] == 0
                    and bassigned[ri] == 0
                    and q_[ri] < 10
                    and q_[ri - 1] >= 20
                    and q_[ri + 1] >= 20
                    and interval_overlap_l(rr_int[ri - 1], ecdna_intervals)
                    == -1
                    and interval_overlap(rr_int[ri + 1], ecdna_intervals) >= 0  # type: ignore[arg-type]
                ):
                    bp_list.append(
                        interval2bp(
                            rr_int[ri - 1],
                            rr_int[ri + 1],
                            (r, ri - 1, ri + 1),
                            int(r_int[ri + 1][0]) - int(r_int[ri - 1][1]),
                        )
                        + [q_[ri - 1], q_[ri + 1]]
                    )

    bp_clusters = cluster_bp_list(
        bp_list, float(normal_cov) * 0.5, bp_match_cutoff_clustering
    )
    bp_refined: list[Any] = []
    bp_stats = []
    for c in bp_clusters:
        if len(c) >= float(normal_cov) * 0.5:
            bp_cluster_r = c
            while len(bp_cluster_r) >= float(normal_cov) * 0.5:
                bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(
                    bp_cluster_r, bp_match_cutoff
                )
                # print (bp[:6])
                if len(set(bpr)) >= float(normal_cov) * 0.5:
                    bpi_ = -1
                    for bpi in range(len(bp_refined)):
                        bp_ = bp_refined[bpi]
                        if (
                            bp[0] == bp_[0]
                            and bp[3] == bp_[3]
                            and bp[2] == bp_[2]
                            and bp[5] == bp_[5]
                            and abs(bp[1] - bp_[1]) <= bp_match_cutoff
                            and abs(bp[4] - bp_[4]) < bp_match_cutoff
                        ):
                            bp_refined[bpi][-1] |= set(bpr)
                            bpi_ = bpi
                            break
                    if bpi_ < 0:
                        bp_refined.append(bp + [bpr])
                        bp_stats.append(bp_stats_)
    print(
        f"Found {len(bp_refined)} breakpoints connecting ecDNA and chromosomes."
    )
    lr_bamfh.close()
    sum_sizes = sum(CHR_SIZES.values())
    agg_size = 0
    xtick_pos = []
    starting_pos = dict()
    for chr in CHR_SIZES.keys():
        agg_size += CHR_SIZES[chr]
        if agg_size < sum_sizes:
            plt.plot(
                [agg_size * 100.0 / sum_sizes, agg_size * 100.0 / sum_sizes],
                [-1, 1000000],
                "k--",
                linewidth=2,
            )
        xtick_pos.append((agg_size - 0.5 * CHR_SIZES[chr]) * 100.0 / sum_sizes)
        starting_pos[chr] = (agg_size - CHR_SIZES[chr]) * 100.0 / sum_sizes

    for bp in bp_refined:
        if (
            interval_overlap_l(
                Interval(bp[0], bp[1], bp[1]), ecdna_intervals_ext
            )
            >= 0
            and interval_overlap_l(
                Interval(bp[3], bp[4], bp[4]), ecdna_intervals_ext
            )
            < 0
        ):
            if bp[3] in starting_pos:
                cn = 0.0
                for seg in cns_dict[bp[3]]:
                    if bp[4] > seg[0] and bp[4] < seg[1]:
                        cn = seg[2]
                        break
                if cn <= 5.0 and len(bp[-1]) <= float(normal_cov) * 2.5:
                    print("Breakpoint", bp[:6], "Support = ", len(bp[-1]))
                    xpos = starting_pos[bp[3]] + bp[4] * 100.0 / sum_sizes
                    ypos = len(bp[-1])
                    plt.plot(xpos, ypos, "bo")
        elif (
            interval_overlap_l(
                Interval(bp[0], bp[1], bp[1]), ecdna_intervals_ext
            )
            < 0
            and interval_overlap_l(
                Interval(bp[3], bp[4], bp[4]), ecdna_intervals_ext
            )
            >= 0
        ):
            if bp[0] in starting_pos:
                cn = 0.0
                for seg in cns_dict[bp[0]]:
                    if bp[1] > seg[0] and bp[1] < seg[1]:
                        cn = seg[2]
                        break
                if cn <= 5.0 and len(bp[-1]) <= float(normal_cov) * 2.5:
                    print("Breakpoint", bp[:6], "Support = ", len(bp[-1]))
                    xpos = starting_pos[bp[0]] + bp[1] * 100.0 / sum_sizes
                    ypos = len(bp[-1])
                    plt.plot(xpos, ypos, "bo")

    plt.xlim([0, 100])
    plt.ylim([1, 500])
    plt.yscale("log")
    plt.xticks(xtick_pos, list(range(1, 23)) + ["X", "Y"])  # type: ignore[arg-type]
    plt.title(output_prefix + " integration loci", fontsize=25)
    plt.ylabel("Long read support", fontsize=25)
    plt.tight_layout()
    out_img_name = "integration_sites_" + output_prefix
    plt.savefig(out_img_name + ".png")
    print("\nCreated " + out_img_name + ".png")
