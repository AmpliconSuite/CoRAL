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
    interval_overlap_l,
)
from coral.constants import CHR_SIZES, CHR_TAG_TO_IDX
from coral.datatypes import BPAlignments, Interval

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
        chimeric_alignments[r] = cigar_parsing.alignment_from_satags(
            chimeric_alignments[r], r
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
        cas = chimeric_alignments[r]
        cas = [ca for ca in cas if ca.ref_interval.chr in CHR_TAG_TO_IDX]
        ref_intvs = [ca.ref_interval for ca in cas]
        for interval in ecdna_intervals:
            i = interval_overlap_l(interval, ref_intvs)
            if i is not None and ref_intvs[i].does_overlap(interval):
                cycle_flag = True
                break
        if cycle_flag:
            bassigned = [0 for i in range(len(cas) - 1)]
            """
			Breakpoint from local alignment i and i + 1
			"""
            for ri in range(len(cas) - 1):
                if (
                    cas[ri].mapq >= 20
                    and cas[ri + 1].mapq >= 20
                    and interval_overlap_l(cas[ri].ref_interval, ecdna_intervals) is None
                    and interval_overlap_l(cas[ri + 1].ref_interval, ecdna_intervals) is not None
                ) or (
                    cas[ri].mapq >= 20
                    and cas[ri + 1].mapq >= 20
                    and interval_overlap_l(cas[ri].ref_interval, ecdna_intervals) is not None
                    and interval_overlap_l(cas[ri + 1].ref_interval, ecdna_intervals) is None
                ):
                    bp_list.append(
                        interval2bp(
                            cas[ri],
                            cas[ri + 1],
                            BPAlignments(r, ri, ri + 1),
                            cas[ri + 1].query_bounds.start - cas[ri].query_bounds.end,
                        )
                    )
                    bassigned[ri] = 1
            """
			Breakpoint from local alignment i - 1 and i + 1
			"""
            for ri in range(1, len(cas) - 1):
                if (
                    bassigned[ri - 1] == 0
                    and bassigned[ri] == 0
                    and cas[ri].mapq < 10
                    and cas[ri - 1].mapq >= 20
                    and cas[ri + 1].mapq >= 20
                    and interval_overlap_l(cas[ri - 1].ref_interval, ecdna_intervals) is None
                    and interval_overlap_l(cas[ri + 1].ref_interval, ecdna_intervals) is not None
                ) or (
                    bassigned[ri - 1] == 0
                    and bassigned[ri] == 0
                    and cas[ri].mapq < 10
                    and cas[ri - 1].mapq >= 20
                    and cas[ri + 1].mapq >= 20
                    and interval_overlap_l(cas[ri - 1].ref_interval, ecdna_intervals) is not None
                    and interval_overlap_l(cas[ri + 1].ref_interval, ecdna_intervals) is None
                ):
                    bp_list.append(
                        interval2bp(
                            cas[ri - 1],
                            cas[ri + 1],
                            BPAlignments(r, ri - 1, ri + 1),
                            cas[ri + 1].query_bounds.start - cas[ri - 1].query_bounds.end,
                        )
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
                        bp_, reads_ = bp_refined[bpi]
                        if (
                            bp.chr1 == bp_.chr1
                            and bp.chr2 == bp_.chr2
                            and bp.strand1 == bp_.strand1
                            and bp.strand2 == bp_.strand2
                            and abs(bp.pos1 - bp_.pos1) <= bp_match_cutoff
                            and abs(bp.pos2 - bp_.pos2) < bp_match_cutoff
                        ):
                            bp_refined[bpi][1] |= set(bpr)
                            bpi_ = bpi
                            break
                    if bpi_ < 0:
                        bp_refined.append([bp, set(bpr)])
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

    for bp_entry in bp_refined:
        bp_obj, bp_reads = bp_entry
        if (
            interval_overlap_l(
                Interval(bp_obj.chr1, bp_obj.pos1, bp_obj.pos1), ecdna_intervals_ext
            )
            is not None
            and interval_overlap_l(
                Interval(bp_obj.chr2, bp_obj.pos2, bp_obj.pos2), ecdna_intervals_ext
            )
            is None
        ):
            if bp_obj.chr2 in starting_pos:
                cn = 0.0
                for seg in cns_dict[bp_obj.chr2]:
                    if bp_obj.pos2 > seg[0] and bp_obj.pos2 < seg[1]:
                        cn = seg[2]
                        break
                if cn <= 5.0 and len(bp_reads) <= float(normal_cov) * 2.5:
                    print("Breakpoint", bp_obj, "Support = ", len(bp_reads))
                    xpos = starting_pos[bp_obj.chr2] + bp_obj.pos2 * 100.0 / sum_sizes
                    ypos = len(bp_reads)
                    plt.plot(xpos, ypos, "bo")
        elif (
            interval_overlap_l(
                Interval(bp_obj.chr1, bp_obj.pos1, bp_obj.pos1), ecdna_intervals_ext
            )
            is None
            and interval_overlap_l(
                Interval(bp_obj.chr2, bp_obj.pos2, bp_obj.pos2), ecdna_intervals_ext
            )
            is not None
        ):
            if bp_obj.chr1 in starting_pos:
                cn = 0.0
                for seg in cns_dict[bp_obj.chr1]:
                    if bp_obj.pos1 > seg[0] and bp_obj.pos1 < seg[1]:
                        cn = seg[2]
                        break
                if cn <= 5.0 and len(bp_reads) <= float(normal_cov) * 2.5:
                    print("Breakpoint", bp_obj, "Support = ", len(bp_reads))
                    xpos = starting_pos[bp_obj.chr1] + bp_obj.pos1 * 100.0 / sum_sizes
                    ypos = len(bp_reads)
                    plt.plot(xpos, ypos, "bo")

    plt.xlim([0, 100])
    plt.ylim([1, 500])
    plt.yscale("log")
    plt.xticks(xtick_pos, list(range(1, 23)) + ["X", "Y"])  # type: ignore[arg-type]
    plt.title(output_prefix + " integration loci", fontsize=25)
    plt.ylabel("Long read support", fontsize=25)
    plt.tight_layout()
    out_img_name = output_prefix + "_integration_sites"
    plt.savefig(out_img_name + ".png")
    print("\nCreated " + out_img_name + ".png")
