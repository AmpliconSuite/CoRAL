#!/usr/bin/env python3
from __future__ import annotations

import functools
import importlib.resources
import io
import logging
import os
import pathlib
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import DefaultDict, Optional

import colorama
import intervaltree
import matplotlib as mpl
import numpy as np
import typer

from coral import datatypes
from coral.breakpoint import (
    breakpoint_utilities,  # type: ignore[import-untyped]
)
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.datatypes import Interval
from coral.summary.parsing import parse_cycle_file

mpl.use("Agg")
import matplotlib.pyplot as plt
import pysam
from matplotlib import gridspec, ticker
from matplotlib.patches import Arc, Rectangle
from pylab import rcParams  # type: ignore[import-untyped]

from coral import (
    bam_types,
    core_types,
    core_utils,
    supplemental_data,
)

rcParams["pdf.fonttype"] = 42


logger = logging.getLogger(__name__)


# makes a gene object from parsed refGene data
# this stores global properties for the gene
class Gene:
    def __init__(self, gchrom, gstart, gend, gdata):
        self.gchrom = gchrom
        self.gstart = gstart
        self.gend = gend
        self.gname = gdata[-4]
        self.strand = gdata[3]
        self.height = 0.5
        # self.highlight_name = highlight_name
        estarts = [int(x) for x in gdata[9].rsplit(",") if x]
        eends = [int(x) for x in gdata[10].rsplit(",") if x]
        self.eposns = list(zip(estarts, eends))

    def __str__(self):
        return f"Gene Name: {self.gname}, Chromosome: {self.gchrom}, Start: {self.gstart}, End: {self.gend}, Strand: {self.strand}"


@dataclass
class GraphViz:
    """ """

    lr_bamfh: Optional[pysam.AlignmentFile] = None
    bam: Optional[bam_types.BAMWrapper] = None
    graph: datatypes.BreakpointGraph | None = None

    graph_amplified_intervals: dict[core_types.ChrTag, list[datatypes.Interval]] = field(
        default_factory=lambda: defaultdict(list)
    )
    num_amplified_intervals: int = 0
    cycle_amplified_intervals: dict[
        core_types.ChrTag, list[datatypes.Interval]
    ] = field(default_factory=lambda: defaultdict(list))
    discordant_edges: list[datatypes.BreakpointEdge] = field(
        default_factory=list
    )

    cycles: dict[int, datatypes.ReconstructedCycle] = field(
        default_factory=dict
    )
    genes: DefaultDict[str, intervaltree.IntervalTree] = field(
        default_factory=lambda: defaultdict(intervaltree.IntervalTree)
    )
    plot_bounds: tuple[str, int, int] | None = None

    def open_bam(self, bam_fn: str) -> None:
        self.lr_bamfh = pysam.AlignmentFile(bam_fn, "rb")
        self.bam = bam_types.BAMWrapper(bam_fn, "rb")

    @property
    def sequence_edges_by_chr(self) -> dict[str, list[datatypes.SequenceEdge]]:
        seq_edges_by_chr: dict[str, list[datatypes.SequenceEdge]] = defaultdict(
            list
        )
        for seq_edge in self.graph.sequence_edges:  # type: ignore[union-attr]
            seq_edges_by_chr[seq_edge.chr].append(seq_edge)
        return seq_edges_by_chr

    def parse_genes(
        self,
        ref_genome: core_types.ReferenceGenome,
        gene_subset_list=None,
        restrict_to_bushman=False,
    ) -> None:
        bushman_set = set()
        if restrict_to_bushman:
            bushman_filepath = (
                importlib.resources.files(supplemental_data)
                / "Bushman_group_allOnco_May2018.tsv"
            )
            with bushman_filepath.open("r") as infile:
                _ = next(infile)
                for line in infile:
                    if not (bushman_fields := line.rstrip().rsplit()):
                        continue
                    bushman_set.add(bushman_fields[-1].strip('"'))

        seen_names = set()
        ref_gene_filepath = (
            importlib.resources.files(supplemental_data)
            / f"refGene_{ref_genome}.txt"
        )
        with ref_gene_filepath.open("r") as infile:
            for line in infile:
                if not (fields := line.rstrip().rsplit()):
                    continue
                curr_chrom: str = fields[2]
                if ref_genome in {
                    core_types.ReferenceGenome.hg19,
                    core_types.ReferenceGenome.hg38,
                } and not curr_chrom.startswith("chr"):
                    curr_chrom = "chr" + curr_chrom

                tstart = int(fields[4])
                tend = int(fields[5])
                gname = fields[-4]
                is_other_feature = gname.startswith(("LOC", "LINC", "MIR"))
                if (
                    (restrict_to_bushman
                    and gname not in bushman_set)
                    or (gene_subset_list
                    and gname not in gene_subset_list)
                ):
                    continue

                if gname not in seen_names and not is_other_feature:
                    seen_names.add(gname)
                    curr_gene = Gene(curr_chrom, tstart, tend, fields)
                    self.genes[curr_chrom][tstart:tend] = curr_gene

    def update_graph_intervals(self) -> None:
        for chrom in self.sequence_edges_by_chr:
            lstart, lend = -2, -2
            if chrom not in self.graph_amplified_intervals:
                self.graph_amplified_intervals[chrom] = []
            for seq_edge in self.sequence_edges_by_chr[chrom]:
                start = seq_edge.start
                end = seq_edge.end
                if start != lend + 1:
                    if lstart >= 0:
                        self.graph_amplified_intervals[chrom].append(
                            datatypes.Interval(chrom, lstart, lend)
                        )
                        self.num_amplified_intervals += 1
                    lstart = start
                    lend = end
                else:
                    lend = end
            self.graph_amplified_intervals[chrom].append(
                datatypes.Interval(chrom, lstart, lend)
            )
            self.num_amplified_intervals += 1

    def merge_intervals(
        self, interval_list: list[tuple[int, int]], padding: float = 0.0
    ) -> list[tuple[int, int]]:
        # takes a list of interval tuples (pos1, pos2). Assumes from same chrom
        # return a list of interval tuples that are merged if overlapping or directly adjacent, or within padding distance
        sorted_intervals = sorted(interval_list)
        merged = [sorted_intervals[0]]
        for current in sorted_intervals[1:]:
            prev = merged[-1]
            if current[0] <= prev[1] + padding:
                merged[-1] = (prev[0], max(prev[1], current[1]))
            else:
                merged.append(current)

        return merged

    def update_cycle_amplified_intervals(
        self,
        cycle_ids: list[int] | None = None,
        *,
        cycle_only: bool = False,
        graph_given: bool = False,
    ) -> None:
        """Derive amplified intervals from (selected) cycles"""
        self.num_amplified_intervals = 0
        if cycle_ids is None:
            cycle_ids = list(self.cycles.keys())
        if cycle_only:
            cycle_ids = [
                cycle_id
                for cycle_id in self.cycles
                if self.cycles[cycle_id].is_cyclic
            ]
        if graph_given:  # if the graph file is given, use this to set the amplified intervals
            for cycle_id in cycle_ids:
                for cycle_seg in self.cycles[cycle_id].segments:
                    for graph_intv in self.graph_amplified_intervals[cycle_seg.chr]:
                        if not graph_intv.encompasses(cycle_seg):
                            continue

                        if (
                            graph_intv
                            not in self.cycle_amplified_intervals[cycle_seg.chr]
                        ):
                            self.cycle_amplified_intervals[
                                cycle_seg.chr
                            ].append(graph_intv)
                        break

        else:  # if the graph file is not given extract from the cycles file
            # collect a list of intervals for each chrom
            cycle_ivald = defaultdict(list)
            for cycle_id in self.cycles:
                for cycle_seg in self.cycles[cycle_id].segments:
                    cycle_ivald[cycle_seg.chr].append(
                        (cycle_seg.start, cycle_seg.end)
                    )

                # merge
                for chrom, ival_list in cycle_ivald.items():
                    merged = self.merge_intervals(ival_list, padding=10000)
                    self.cycle_amplified_intervals[chrom] = [
                        datatypes.Interval(chrom, start, end)
                        for start, end in merged
                    ]

        for chr_tag in self.cycle_amplified_intervals:
            self.cycle_amplified_intervals[chr_tag] = sorted(
                self.cycle_amplified_intervals[chr_tag]
            )
            self.num_amplified_intervals += len(
                self.cycle_amplified_intervals[chr_tag]
            )

    def set_gene_heights(
        self, rel_genes: list[Gene], padding: float = 0.0
    ) -> None:
        if not rel_genes:
            return
        gname_to_gobj = {x.gname: x for x in rel_genes}
        # merge intervals
        intervals = [(x.gstart, x.gend) for x in rel_genes]
        merged = self.merge_intervals(intervals, padding=padding)

        gene_ival_t = intervaltree.IntervalTree[str]()
        for x in rel_genes:
            gene_ival_t.addi(x.gstart, x.gend, x.gname)

        for mi in merged:
            ghits: list[intervaltree.Interval] = gene_ival_t[mi[0] : mi[1]]
            gene_heights = np.linspace(0.15, 0.75, len(ghits))
            for g, h in zip(ghits, gene_heights):
                gname_to_gobj[g.data].height = h

    def plot_graph(
        self,
        title: str,
        output_fn: str,
        margin_between_intervals: float = 2,
        height: float = 7.5,
        fontsize: float = 18,
        dpi: int = 300,
        max_cov_cutoff: float = float("inf"),
        quality_threshold: float = 0,
        gene_font_size: float = 12,
        *,
        hide_genes: bool = False,
    ) -> None:
        """Plot discordant edges and coverage on sequence edges in breakpoint
        graph."""
        if not self.plot_bounds:
            width: int = max(15, 2 * self.num_amplified_intervals)
        else:
            width = 15
        fig = plt.figure(figsize=(width, height))
        if not hide_genes:
            gs = gridspec.GridSpec(2, 1, height_ratios=[8, 2])
        else:
            gs = gridspec.GridSpec(2, 1, height_ratios=[8, 0.000001])
        ax = fig.add_subplot(gs[0, 0])
        plt.subplots_adjust(
            left=73 / 1000.0,
            right=1 - 73 / 1000.0,
            bottom=1 / 4.0,
            top=1 - 1 / 20.0,
        )
        ax.set_title(title, fontsize=fontsize)
        ax2 = ax.twinx()
        ax3 = fig.add_subplot(gs[1, 0], sharex=ax)
        # ax.yaxis.set_label_coords(-0.05, 0.25)
        # ax2.yaxis.set_label_coords(1.05, 0.33)
        ax.xaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.spines["left"].set_visible(False)
        ax3.spines["right"].set_visible(False)
        ax3.spines["top"].set_visible(False)

        # Draw sequence edges
        total_len_amp = 0  # Total length of amplified intervals
        for chrom in self.graph_amplified_intervals:
            total_len_amp += sum(
                [len(graph_intv) for graph_intv in self.graph_amplified_intervals[chrom]],
            )
        # sorted_chrs = sorted(self.intervals_from_graph.keys(), key = lambda chr: CHR_TAG_TO_IDX[chr])
        zoom_factor = 1.0
        if self.plot_bounds:
            zoom_factor = (
                float(self.plot_bounds[2] - self.plot_bounds[1]) / total_len_amp
            )
        sorted_chrs = breakpoint_utilities.sort_chrom_names(
            self.graph_amplified_intervals.keys()
        )
        amplified_intervals_start = {}
        ymax = 0
        ylim_params = [0, 0]
        x = margin_between_intervals
        for chrom in sorted_chrs:
            interval_idx = 0
            amplified_intervals_start[chrom] = [x]
            for seq in self.sequence_edges_by_chr[chrom]:
                if chrom not in self.graph_amplified_intervals or (
                    interval_idx >= len(self.graph_amplified_intervals[chrom])
                    or seq.start > self.graph_amplified_intervals[chrom][interval_idx].end
                ):
                    # int_ = self.intervals_from_graph[chrom][interval_idx]
                    x += margin_between_intervals
                    amplified_intervals_start[chrom].append(x)
                    interval_idx += 1
                x1 = x
                x += (seq.end - seq.start) * 100.0 / total_len_amp
                x2 = x
                if self.plot_bounds:
                    if chrom != self.plot_bounds[0]:
                        continue  # Skip if chromosome doesn't match plot bounds

                    if not (
                        seq.end >= self.plot_bounds[1]
                        and seq.start <= self.plot_bounds[2]
                    ):
                        continue  # Skip if interval doesn't overlap with plot bounds

                y = seq.cn
                ylim_params[0] += (seq.cn ** 2)
                ylim_params[1] += (seq.cn * seq.lr_nc / 1.25)
                ymax = max(y, ymax)

                ax2.hlines(y, x1, x2, color="black", lw=6, zorder=2)

            x += margin_between_intervals

        # Draw amplified interval separators
        if not self.plot_bounds:
            for chrom in amplified_intervals_start:
                if chrom != sorted_chrs[0]:
                    ax.axvline(
                        x=amplified_intervals_start[chrom][0]
                        - margin_between_intervals * 0.5,
                        linestyle="--",
                        lw=2,
                        zorder=2,
                    )
                    ax3.axvline(
                        x=amplified_intervals_start[chrom][0]
                        - margin_between_intervals * 0.5,
                        linestyle="--",
                        lw=2,
                        zorder=2,
                    )
                for i in range(1, len(amplified_intervals_start[chrom])):
                    ax.axvline(
                        x=amplified_intervals_start[chrom][i]
                        - margin_between_intervals * 0.5,
                        linestyle=":",
                        lw=2,
                        zorder=2,
                    )

        # Draw discordant edges
        colorcode = {
            "+-": "red",
            "++": "magenta",
            "-+": (139 / 256.0, 69 / 256.0, 19 / 256.0),
            "--": "teal",
            "interchromosomal": "blue",
        }
        assert self.graph
        avg_bp_rc = (
            sum([bp.lr_count for bp in self.graph.discordant_edges])
            * 1.0
            / max(len(self.graph.discordant_edges), 1)
        )
        for bp in self.graph.discordant_edges:
            chr1 = bp.node1.chr
            pos1 = bp.node1.pos
            chr2 = bp.node2.chr
            pos2 = bp.node2.pos
            int1 = 0
            int2 = 0
            ort = f"{bp.node1.strand}{bp.node2.strand}"
            if chr1 in self.graph_amplified_intervals and chr2 in self.graph_amplified_intervals:
                while pos1 > self.graph_amplified_intervals[chr1][int1].end:
                    int1 += 1
                bp_x1 = (
                    amplified_intervals_start[chr1][int1]
                    + (pos1 - self.graph_amplified_intervals[chr1][int1].start)
                    * 100.0
                    / total_len_amp
                )
                while pos2 > self.graph_amplified_intervals[chr2][int2].end:
                    int2 += 1
                bp_x2 = (
                    amplified_intervals_start[chr2][int2]
                    + (pos2 - self.graph_amplified_intervals[chr2][int2].start)
                    * 100.0
                    / total_len_amp
                )
                # check if either bp overlaps before plotting
                if self.plot_bounds:
                    # Check if both breakpoints belong to the same chromosome
                    # as in plot bounds
                    hit1 = (
                        chr1 == self.plot_bounds[0]
                        and self.plot_bounds[1] <= pos1 <= self.plot_bounds[2]
                    )
                    hit2 = (
                        chr2 == self.plot_bounds[0]
                        and self.plot_bounds[1] <= pos2 <= self.plot_bounds[2]
                    )
                    if not hit1 and not hit2:
                        continue

                arc = Arc(
                    ((bp_x1 + bp_x2) * 0.5, 0),
                    bp_x1 - bp_x2,
                    2 * ymax,
                    theta1=0,
                    theta2=180,
                    color=colorcode[ort],
                    lw=min(3 * (bp.lr_count / avg_bp_rc), 3),
                    zorder=3,
                )
                ax2.add_patch(arc)

            else:
                print("Could not place " + str(bp))
                continue

        #ax2.set_ylim(0, 1.4 * ymax)
        ax2.set_ylabel("CN", fontsize=fontsize)
        ax2.tick_params(axis="y", labelsize=fontsize)

        # Draw coverage within amplified intervals
        max_cov = 0
        for chrom in sorted_chrs:
            for inti in range(len(self.graph_amplified_intervals[chrom])):
                graph_intv = self.graph_amplified_intervals[chrom][inti]
                if self.plot_bounds:
                    if chrom != self.plot_bounds[0]:
                        continue  # Skip if chromosome doesn't match plot bounds

                    if not (
                        graph_intv.end >= self.plot_bounds[1]
                        and graph_intv.start <= self.plot_bounds[2]
                    ):
                        continue  # Skip if interval doesn't overlap with plot bounds

                window_size = 150
                ival_len = graph_intv.end - graph_intv.start
                if self.plot_bounds:
                    ival_len = self.plot_bounds[2] - self.plot_bounds[1]

                if ival_len >= 1_000_000:
                    window_size = 10_000
                elif ival_len >= 100_000:
                    window_size = 1_000

                for w in range(graph_intv.start, graph_intv.end, window_size):
                    intv = Interval(chrom, w, w + window_size)
                    cov = (
                        self.bam.count_raw_coverage(  # type: ignore
                            intv,
                            quality_threshold=0,
                            read_callback_type="nofilter",
                        )
                        / window_size
                    )
                    max_cov = max(cov, max_cov)
                    x = (
                        amplified_intervals_start[chrom][inti]
                        + (w - graph_intv.start) * 100.0 / total_len_amp
                    )
                    rect = Rectangle(
                        (x, 0),
                        window_size * 100.0 / total_len_amp,
                        cov,
                        color="silver",
                        zorder=1,
                    )
                    ax.add_patch(rect)
                w = graph_intv.end - (
                    (graph_intv.end - graph_intv.start + 1) % window_size
                )
                if w < graph_intv.end:
                    cov = (
                        self.bam.count_raw_coverage(
                            Interval(chrom, w, w + window_size)
                        )
                        * 1.0
                        / window_size
                    )
                    max_cov = max(cov, max_cov)
                    x = (
                        amplified_intervals_start[chrom][inti]
                        + (w - graph_intv.start) * 100.0 / total_len_amp
                    )
                    rect = Rectangle(
                        (x, 0),
                        window_size * 100.0 / total_len_amp,
                        cov,
                        color="silver",
                        zorder=1,
                    )
                    ax.add_patch(rect)
        ylim_params[1] /= max_cov
        ax2.set_ylim(0, ylim_params[0] / ylim_params[1])
        ax.set_ylabel("Coverage", fontsize=fontsize)
        ax.set_ylim(0, min(1.25 * max_cov, max_cov_cutoff))
        ax.tick_params(axis="y", labelsize=fontsize)

        # draw genes below plot
        if not hide_genes:
            for chrom in sorted_chrs:
                for inti in range(len(self.graph_amplified_intervals[chrom])):
                    graph_intv = self.graph_amplified_intervals[chrom][inti]
                    if self.plot_bounds:
                        if chrom != self.plot_bounds[0]:
                            continue  # Skip if chromosome doesn't match plot bounds

                        if not (
                            graph_intv.end >= self.plot_bounds[1]
                            and graph_intv.start <= self.plot_bounds[2]
                        ):
                            continue  # Skip if interval doesn't overlap with plot bounds

                    rel_genes = [
                        x.data
                        for x in self.genes[chrom][
                            graph_intv.start : graph_intv.end
                        ]
                    ]
                    gene_padding = total_len_amp * 0.02
                    self.set_gene_heights(rel_genes, gene_padding)

                    for gene_obj in rel_genes:  # plot line for the gene
                        height = gene_obj.height
                        cut_gs = max(graph_intv.start, gene_obj.gstart)
                        cut_ge = min(graph_intv.end, gene_obj.gend)
                        # if self.plot_bounds:
                        #     cut_gs = max(self.plot_bounds[1], cut_gs)
                        #     cut_ge = min(self.plot_bounds[2], cut_ge)

                        gene_start = (
                            amplified_intervals_start[chrom][inti]
                            + (cut_gs - graph_intv.start)
                            * 100.0
                            / total_len_amp
                        )
                        gene_end = (
                            amplified_intervals_start[chrom][inti]
                            + (cut_ge - graph_intv.start)
                            * 100.0
                            / total_len_amp
                        )
                        ax3.hlines(
                            height,
                            gene_start,
                            gene_end,
                            color="cornflowerblue",
                            lw=4.5,
                        )  # Draw horizontal bars for genes
                        if self.plot_bounds:
                            if (
                                cut_ge < self.plot_bounds[1]
                                or cut_gs > self.plot_bounds[2]
                            ):
                                continue

                            cut_gs = max(self.plot_bounds[1], cut_gs)
                            cut_ge = min(self.plot_bounds[2], cut_ge)
                            gene_start = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_gs - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )
                            gene_end = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_ge - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )

                        ax3.text(
                            (gene_start + gene_end) / 2,
                            height + 0.05,
                            gene_obj.gname,
                            ha="center",
                            va="bottom",
                            fontsize=gene_font_size,
                            style="italic",
                        )

                        if gene_obj.strand == "+":
                            ax3.plot(
                                gene_start,
                                height,
                                marker=">",
                                color="black",
                                markersize=7,
                            )
                        elif gene_obj.strand == "-":
                            ax3.plot(
                                gene_end,
                                height,
                                marker="<",
                                color="black",
                                markersize=7,
                            )

                        for (
                            exon_start,
                            exon_end,
                        ) in gene_obj.eposns:  # plot exon bars
                            if (
                                not exon_end > graph_intv.start
                                or not exon_start < graph_intv.end
                            ):
                                continue

                            cut_es = max(graph_intv.start, exon_start)
                            cut_ee = min(graph_intv.end, exon_end)
                            # if self.plot_bounds:
                            #     cut_es = max(self.plot_bounds[1], cut_es)
                            #     cut_ee = min(self.plot_bounds[2], cut_ee)
                            #
                            exon_start_pos = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_es - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )
                            exon_end_pos = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_ee - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )

                            exon_min_width = (
                                0.2 * zoom_factor
                            )  # Adjust the minimum width as needed
                            exon_width = exon_end_pos - exon_start_pos
                            if exon_width < exon_min_width:
                                diff = (exon_min_width - exon_width) / 2
                                exon_start_pos -= diff
                                exon_end_pos += diff

                            ax3.hlines(
                                height,
                                exon_start_pos,
                                exon_end_pos,
                                color="black",
                                lw=7.5,
                            )

        # Ticks and labels
        xtickpos = []
        xticklabels = []
        if not self.plot_bounds:
            ax.set_xlim(
                0,
                100
                + (self.num_amplified_intervals + 1) * margin_between_intervals,
            )
            ax2.set_xlim(
                0,
                100
                + (self.num_amplified_intervals + 1) * margin_between_intervals,
            )
            ax3.set_xlim(
                0,
                100
                + (self.num_amplified_intervals + 1) * margin_between_intervals,
            )
            for chrom in sorted_chrs:
                nint_chr = len(self.graph_amplified_intervals[chrom])
                for inti in range(len(amplified_intervals_start[chrom])):
                    if inti > 0:
                        xtickpos.append(
                            amplified_intervals_start[chrom][inti]
                            - margin_between_intervals,
                        )
                        # Add chr label x-axis tick in the middle of all intervals
                        if (
                            nint_chr % 2 == 0
                            and inti == (nint_chr - 2) // 2 + 1
                        ):
                            xtickpos.append(
                                amplified_intervals_start[chrom][inti]
                                - margin_between_intervals * 0.5,
                            )
                        xtickpos.append(amplified_intervals_start[chrom][inti])
                        # Add chr label x-axis tick in the middle of all intervals
                        if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                            xtickpos.append(
                                (
                                    amplified_intervals_start[chrom][inti]
                                    + amplified_intervals_start[chrom][inti + 1]
                                    - margin_between_intervals
                                )
                                * 0.5,
                            )
                    else:
                        if chrom != sorted_chrs[0]:
                            xtickpos.append(
                                amplified_intervals_start[chrom][0]
                                - margin_between_intervals,
                            )
                        xtickpos.append(amplified_intervals_start[chrom][0])
                        if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                            chri = sorted_chrs.index(chrom)
                            if chri == len(sorted_chrs) - 1:
                                amplified_intervals_end = (
                                    100
                                    + self.num_amplified_intervals
                                    * margin_between_intervals
                                )
                            else:
                                amplified_intervals_end = (
                                    amplified_intervals_start[
                                        sorted_chrs[chri + 1]
                                    ][0]
                                    - margin_between_intervals
                                )
                            xtickpos.append(
                                (
                                    amplified_intervals_start[chrom][inti]
                                    + amplified_intervals_end
                                )
                                * 0.5,
                            )
            xtickpos.append(
                100 + self.num_amplified_intervals * margin_between_intervals
            )
            for chrom in sorted_chrs:
                nint_chr = len(self.graph_amplified_intervals[chrom])
                for inti in range(nint_chr):
                    graph_intv = self.graph_amplified_intervals[chrom][inti]
                    xticklabels.append(f"{graph_intv.start:,}   ")
                    if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                        xticklabels.append(chrom)
                    xticklabels.append(f"{graph_intv.end:,}   ")
                    if nint_chr % 2 == 0 and inti == (nint_chr - 2) // 2:
                        xticklabels.append(chrom)
            ax3.set_xticks(xtickpos)
            ax3.set_xticklabels(xticklabels, size=fontsize)
            ticks_labels = ax3.get_xticklabels()
            for ti in range(len(xticklabels)):
                if xticklabels[ti] not in sorted_chrs:
                    ticks_labels[ti].set_rotation(90)
                else:
                    ax3.xaxis.get_major_ticks()[ti].tick1line.set_visible(False)

        else:  # self.plot_bounds are given
            # look up the segment
            pchrom, pstart, pend = self.plot_bounds
            nint_chr = len(self.graph_amplified_intervals[pchrom])
            relint = None
            rint_ = None
            for inti in range(nint_chr):
                istart, iend = (
                    self.graph_amplified_intervals[pchrom][inti].start,
                    self.graph_amplified_intervals[pchrom][inti].end,
                )
                if istart <= pstart <= iend:
                    relint = inti
                    rint_ = self.graph_amplified_intervals[pchrom][inti]
                    break

            if relint is None:
                print(
                    f"Could not identify region {pchrom}:{pstart}-{pend} in graph regions. Region should be fully contained in graph.",
                )

            else:
                plot_start = (
                    amplified_intervals_start[pchrom][relint]
                    + (pstart - rint_.start) * 100.0 / total_len_amp
                )
                plot_end = (
                    amplified_intervals_start[pchrom][relint]
                    + (pend - rint_.start) * 100.0 / total_len_amp
                )
                xtickpos.append(plot_start)
                xtickpos.append(plot_end)
                xticklabels.append(pchrom + ":" + str(pstart))
                xticklabels.append(pchrom + ":" + str(pend))
                ax3.set_xticks(xtickpos)
                ax3.set_xticklabels(xticklabels, size=fontsize - 4)

                ax.set_xlim(plot_start, plot_end)
                ax2.set_xlim(plot_start, plot_end)
                ax3.set_xlim(plot_start, plot_end)

        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
        ax3.set_ylim(0, 1)
        fig.subplots_adjust(hspace=0)
        plt.savefig(output_fn + ".png", dpi=dpi)
        plt.savefig(output_fn + ".pdf")

    def close_bam(self):
        try:
            self.lr_bamfh.close()
        except AttributeError:
            pass

    def plot_cycles(
        self,
        title: str,
        output_fn: str,
        num_cycles: int | None = None,
        margin_between_intervals: int = 2,
        fontsize: int = 18,
        dpi: int = 300,
        gene_font_size: int = 12,
        *,
        cycle_only: bool = False,
        hide_genes: bool = False,
    ) -> None:
        """Plot cycles & paths returned from cycle decomposition"""
        width = max(15, 2 * self.num_amplified_intervals)
        cycles_to_plot = list(self.cycles)  # Get cycle ID keys
        if num_cycles is not None:
            cycles_to_plot = [
                cycle_id
                for cycle_id in self.cycles
                if int(cycle_id) <= num_cycles
            ]
        if cycle_only:
            cycles_to_plot = [
                cycle_id
                for cycle_id in cycles_to_plot
                if self.cycles[cycle_id].is_cyclic
            ]
        cycles_to_plot = sorted(cycles_to_plot)
        height = sum(
            [
                3 * len(self.cycles[cycle_id].segments) - 1
                for cycle_id in cycles_to_plot
            ]
        ) + 9 * (len(cycles_to_plot) - 1)
        fig = plt.figure(figsize=(width, max(4, height * 0.25)))
        if not hide_genes:
            vrat = 50 / height
            gs = gridspec.GridSpec(2, 1, height_ratios=[8, vrat])
        else:
            gs = gridspec.GridSpec(2, 1, height_ratios=[8, 0.000001])
        ax = fig.add_subplot(gs[0, 0])
        ax.set_title(title, fontsize=fontsize)
        ax.xaxis.set_visible(False)
        ax3 = fig.add_subplot(gs[1, 0], sharex=ax)
        ax3.yaxis.set_visible(False)
        ax3.spines["left"].set_visible(False)
        ax3.spines["right"].set_visible(False)
        ax3.spines["top"].set_visible(False)

        # Compute the x coordinates for each amplified interval
        total_len_amp = 0  # Total length of amplified intervals
        for chrom in self.cycle_amplified_intervals:
            total_len_amp += sum(
                [
                    int_.end - int_.start + 1
                    for int_ in self.cycle_amplified_intervals[chrom]
                ],
            )
        # sorted_chrs = sorted(self.intervals_from_cycle.keys(), key = lambda chr: CHR_TAG_TO_IDX[chr])
        sorted_chrs = breakpoint_utilities.sort_chrom_names(
            self.cycle_amplified_intervals.keys()
        )
        amplified_intervals_start = {}
        x: float = margin_between_intervals
        for chrom in sorted_chrs:
            amplified_intervals_start[chrom] = [x]
            for interval_idx in range(
                len(self.cycle_amplified_intervals[chrom])
            ):
                intv = self.cycle_amplified_intervals[chrom][interval_idx]
                x += (intv.end - intv.start) * 100.0 / total_len_amp
                x += margin_between_intervals
                if (
                    interval_idx
                    < len(self.cycle_amplified_intervals[chrom]) - 1
                ):
                    amplified_intervals_start[chrom].append(x)

        # Draw amplified interval separators
        for chrom in amplified_intervals_start:
            if chrom != sorted_chrs[0]:
                ax.axvline(
                    x=amplified_intervals_start[chrom][0]
                    - margin_between_intervals * 0.5,
                    linestyle="--",
                    lw=2,
                )
                ax3.axvline(
                    x=amplified_intervals_start[chrom][0]
                    - margin_between_intervals * 0.5,
                    linestyle="--",
                    lw=2,
                )
            for i in range(1, len(amplified_intervals_start[chrom])):
                ax.axvline(
                    x=amplified_intervals_start[chrom][i]
                    - margin_between_intervals * 0.5,
                    linestyle=":",
                    lw=2,
                )

        # Draw cycles
        y_cur = -2
        extension = 1.5
        cycleticks = []
        cycleticklabels = []
        for cycle_id in cycles_to_plot:
            ystart_cycle_id = y_cur
            cycle_min_x = float("inf")
            cycle_max_x = 0.0
            for i, segment in enumerate(self.cycles[cycle_id].segments):
                # Segment i
                interval_idx = 0
                while (
                    segment.start
                    > self.cycle_amplified_intervals[segment.chr][
                        interval_idx
                    ].end
                ):
                    interval_idx += 1
                x1 = (
                    amplified_intervals_start[segment.chr][interval_idx]
                    + (
                        segment.start
                        - self.cycle_amplified_intervals[segment.chr][
                            interval_idx
                        ].start
                    )
                    * 100.0
                    / total_len_amp
                )
                cycle_min_x = min(x1, cycle_min_x)
                xlen = (segment.end - segment.start) * 100.0 / total_len_amp
                cycle_max_x = max(x1 + xlen, cycle_max_x)
                rect = Rectangle(
                    (x1, y_cur),
                    xlen,
                    1,
                    facecolor="antiquewhite",
                    linewidth=2,
                    edgecolor="dimgrey",
                )
                ax.add_patch(rect)

                # Connections between segment i and i + 1
                if i < len(self.cycles[cycle_id].segments) - 1:
                    nseg = self.cycles[cycle_id].segments[i + 1]
                    interval_idx_n = 0
                    while (
                        nseg.start
                        > self.cycle_amplified_intervals[nseg.chr][
                            interval_idx_n
                        ].end
                    ):
                        interval_idx_n += 1
                    if segment.strand == "+" and nseg.strand == "-":
                        x2 = x1 + xlen
                        x2n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x2n += (
                            (
                                nseg.end
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        ax.vlines(
                            x=max(x2, x2n) + extension,
                            ymin=y_cur + 0.5,
                            ymax=y_cur - 1.5,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur + 0.5,
                            xmin=x2,
                            xmax=max(x2, x2n) + extension,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur - 1.5,
                            xmin=x2n,
                            xmax=max(x2, x2n) + extension,
                            colors="b",
                            lw=2,
                        )
                        y_cur -= 2
                    elif segment.strand == "-" and nseg.strand == "+":
                        x1n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x1n += (
                            (
                                nseg.start
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        ax.vlines(
                            x=min(x1, x1n) - extension,
                            ymin=y_cur + 0.5,
                            ymax=y_cur - 1.5,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur + 0.5,
                            xmin=min(x1, x1n) - extension,
                            xmax=x1,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur - 1.5,
                            xmin=min(x1, x1n) - extension,
                            xmax=x1n,
                            colors="b",
                            lw=2,
                        )
                        y_cur -= 2
                    elif segment.strand == "+" and nseg.strand == "+":
                        x2 = x1 + xlen
                        x1n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x1n += (
                            (
                                nseg.start
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        if x2 <= x1n:
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x2,
                                xmax=x1n,
                                colors="b",
                                lw=2,
                            )
                        else:
                            ax.vlines(
                                x=x2 + extension,
                                ymin=y_cur - 0.5,
                                ymax=y_cur + 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.vlines(
                                x=x1n - extension,
                                ymin=y_cur - 1.5,
                                ymax=y_cur - 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x2,
                                xmax=x2 + extension,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 0.5,
                                xmin=x1n - extension,
                                xmax=x2 + extension,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 1.5,
                                xmin=x1n - extension,
                                xmax=x1n,
                                colors="b",
                                lw=2,
                            )
                            y_cur -= 2
                    else:  # seg[3] == '-' and nseg[3] == '-'
                        x2n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x2n += (
                            (
                                nseg.end
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        if x1 >= x2n:
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x2n,
                                xmax=x1,
                                colors="b",
                                lw=2,
                            )
                        else:
                            ax.vlines(
                                x=x1 - extension,
                                ymin=y_cur - 0.5,
                                ymax=y_cur + 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.vlines(
                                x=x2n + extension,
                                ymin=y_cur - 1.5,
                                ymax=y_cur - 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x1 - extension,
                                xmax=x1,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 0.5,
                                xmin=x1 - extension,
                                xmax=x2n + extension,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 1.5,
                                xmin=x2n,
                                xmax=x2n + extension,
                                colors="b",
                                lw=2,
                            )
                            y_cur -= 2

            # First and last segments
            if not self.cycles[cycle_id].is_cyclic:  # Paths
                seg = self.cycles[cycle_id].segments[0]
                interval_idx = 0
                while (
                    seg.start
                    > self.cycle_amplified_intervals[seg.chr][interval_idx].end
                ):
                    interval_idx += 1
                if seg.strand == "+":
                    x1: float = (
                        amplified_intervals_start[seg.chr][interval_idx]
                        + (
                            seg.start
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x1 - 2 * extension,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
                else:
                    x2: float = (
                        amplified_intervals_start[seg.chr][interval_idx]
                        + (
                            seg.end
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x2,
                        xmax=x2 + 2 * extension,
                        colors="b",
                        lw=2,
                    )
                seg = self.cycles[cycle_id].segments[-1]
                interval_idx = 0
                while (
                    seg.start
                    > self.cycle_amplified_intervals[seg.chr][interval_idx].end
                ):
                    interval_idx += 1
                if seg.strand == "+":
                    x2: float = amplified_intervals_start[seg.chr][interval_idx]
                    x2 += (
                        (
                            seg.end
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x2,
                        xmax=x2 + 2 * extension,
                        colors="b",
                        lw=2,
                    )
                else:
                    x1: float = amplified_intervals_start[seg.chr][interval_idx]
                    x1 += (
                        (
                            seg.start
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x1 - 2 * extension,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
            else:  # Cycles
                # xmin_ = 0.5
                xmin_ = cycle_min_x - extension
                xmax_ = cycle_max_x + extension

                if len(self.cycles[cycle_id].segments) > 1:
                    xmin_ -= extension
                    xmax_ += extension

                # xmax_ = 99.5 + (self.num_amplified_intervals + 1) * margin_between_intervals
                seg1 = self.cycles[cycle_id].segments[0]
                interval_idx1 = 0
                while (
                    seg1.start
                    > self.cycle_amplified_intervals[seg1.chr][
                        interval_idx1
                    ].end
                ):
                    interval_idx1 += 1
                seg2 = self.cycles[cycle_id].segments[-1]
                interval_idx2 = 0
                while (
                    seg2.start
                    > self.cycle_amplified_intervals[seg2.chr][
                        interval_idx2
                    ].end
                ):
                    interval_idx2 += 1
                if seg1.strand == "-" and seg2.strand == "+":
                    x2: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x2 += (
                        (
                            seg1.end
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x2n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x2n += (
                        (
                            seg2.end
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmax_,
                        ymin=y_cur + 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x2,
                        xmax=xmax_,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5, xmin=x2n, xmax=xmax_, colors="b", lw=2
                    )
                elif seg1.strand == "+" and seg2.strand == "-":
                    x1: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x1 += (
                        (
                            seg1.start
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x1n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x1n += (
                        (
                            seg2.start
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmin_,
                        ymin=y_cur + 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=xmin_,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5, xmin=xmin_, xmax=x1n, colors="b", lw=2
                    )
                elif seg1.strand == "-" and seg2.strand == "-":
                    x2: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x2 += (
                        (
                            seg1.end
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x1n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x1n += (
                        (
                            seg2.end
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmax_,
                        ymin=y_cur - 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.vlines(
                        x=x1n - extension,
                        ymin=y_cur - 0.5,
                        ymax=y_cur + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x2,
                        xmax=xmax_,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x1n - extension,
                        xmax=x1n,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur - 0.5,
                        xmin=x1n - extension,
                        xmax=xmax_,
                        colors="b",
                        lw=2,
                    )
                else:
                    x1: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x1 += (
                        (
                            seg1.start
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x2n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x2n += (
                        (
                            seg2.end
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmin_,
                        ymin=y_cur - 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.vlines(
                        x=x2n + extension,
                        ymin=y_cur - 0.5,
                        ymax=y_cur + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=xmin_,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x2n,
                        xmax=x2n + extension,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur - 0.5,
                        xmin=xmin_,
                        xmax=x2n + extension,
                        colors="b",
                        lw=2,
                    )

            # Separators between cycles; ticks
            ax.hlines(
                y=y_cur - 2,
                xmin=-1,
                xmax=101
                + (self.num_amplified_intervals + 1) * margin_between_intervals,
                colors="k",
            )
            cycleticks.append((y_cur + ystart_cycle_id) * 0.5)
            if self.cycles[cycle_id].is_cyclic:
                cycleticklabels.append(
                    f"cycle {cycle_id}:\n"
                    f"CN = {self.cycles[cycle_id].overall_cn:.2f}\n"
                    f"weighted-CN = {self.cycles[cycle_id].total_cn_weighted_length:,.2e}"
                )
            else:
                cycleticklabels.append(
                    f"path {cycle_id}:\n"
                    f"CN = {self.cycles[cycle_id].overall_cn:.2f}\n"
                    f"weighted-CN = {self.cycles[cycle_id].total_cn_weighted_length:,.2e}"
                )
            y_cur -= 4

        if not hide_genes:
            for chrom in sorted_chrs:
                for int_idx in range(
                    len(self.cycle_amplified_intervals[chrom])
                ):
                    intv = self.cycle_amplified_intervals[chrom][int_idx]
                    rel_genes = [
                        x.data for x in self.genes[chrom][intv.start : intv.end]
                    ]
                    gene_padding = total_len_amp * 0.02
                    self.set_gene_heights(rel_genes, gene_padding)

                    for gene_obj in rel_genes:  # plot gene lines
                        height = gene_obj.height
                        cut_gs = max(intv.start, gene_obj.gstart)
                        cut_ge = min(intv.end, gene_obj.gend)
                        gene_start = (
                            amplified_intervals_start[chrom][int_idx]
                            + (cut_gs - intv.start) * 100.0 / total_len_amp
                        )
                        gene_end = (
                            amplified_intervals_start[chrom][int_idx]
                            + (cut_ge - intv.start) * 100.0 / total_len_amp
                        )
                        ax3.hlines(
                            height,
                            gene_start,
                            gene_end,
                            color="cornflowerblue",
                            lw=4.5,
                        )  # Draw horizontal bars for genes
                        ax3.text(
                            (gene_start + gene_end) / 2,
                            height + 0.05,
                            gene_obj.gname,
                            ha="center",
                            va="bottom",
                            fontsize=gene_font_size,
                            style="italic",
                        )

                        if gene_obj.strand == "+":
                            ax3.plot(
                                gene_start,
                                height,
                                marker=">",
                                color="black",
                                markersize=7,
                            )
                        elif gene_obj.strand == "-":
                            ax3.plot(
                                gene_end,
                                height,
                                marker="<",
                                color="black",
                                markersize=7,
                            )

                        for (
                            exon_start,
                            exon_end,
                        ) in gene_obj.eposns:  # plot exon bars
                            if (
                                not exon_end > intv.start
                                or not exon_start < intv.end
                            ):
                                continue

                            cut_es = max(intv.start, exon_start)
                            cut_ee = min(intv.end, exon_end)
                            exon_start_pos = (
                                amplified_intervals_start[chrom][int_idx]
                                + (cut_es - intv.start) * 100.0 / total_len_amp
                            )
                            exon_end_pos = (
                                amplified_intervals_start[chrom][int_idx]
                                + (cut_ee - intv.start) * 100.0 / total_len_amp
                            )

                            exon_min_width = (
                                0.2  # Adjust the minimum width as needed
                            )
                            exon_width = exon_end_pos - exon_start_pos
                            if exon_width < exon_min_width:
                                diff = (exon_min_width - exon_width) / 2
                                exon_start_pos -= diff
                                exon_end_pos += diff

                            ax3.hlines(
                                height,
                                exon_start_pos,
                                exon_end_pos,
                                color="black",
                                lw=7.5,
                            )

        # Ticks an labels
        ax.set_xlim(
            -1,
            101 + (self.num_amplified_intervals + 1) * margin_between_intervals,
        )
        ax.set_ylim(y_cur + 2, 0)
        xtickpos: list[float | int] = []
        for chrom in sorted_chrs:
            nint_chr = len(self.cycle_amplified_intervals[chrom])
            for int_idx in range(len(amplified_intervals_start[chrom])):
                if int_idx > 0:
                    xtickpos.append(
                        amplified_intervals_start[chrom][int_idx]
                        - margin_between_intervals,
                    )
                    if nint_chr % 2 == 0 and int_idx == (nint_chr - 2) // 2 + 1:
                        xtickpos.append(
                            amplified_intervals_start[chrom][int_idx]
                            - margin_between_intervals * 0.5,
                        )
                    xtickpos.append(amplified_intervals_start[chrom][int_idx])
                    if nint_chr % 2 == 1 and int_idx == (nint_chr - 1) // 2:
                        xtickpos.append(
                            (
                                amplified_intervals_start[chrom][int_idx]
                                + amplified_intervals_start[chrom][int_idx + 1]
                                - margin_between_intervals
                            )
                            * 0.5,
                        )
                else:
                    if chrom != sorted_chrs[0]:
                        xtickpos.append(
                            amplified_intervals_start[chrom][0]
                            - margin_between_intervals,
                        )
                    xtickpos.append(amplified_intervals_start[chrom][0])
                    if nint_chr % 2 == 1 and int_idx == (nint_chr - 1) // 2:
                        chri = sorted_chrs.index(chrom)
                        if chri == len(sorted_chrs) - 1:
                            amplified_intervals_end = (
                                100
                                + self.num_amplified_intervals
                                * margin_between_intervals
                            )
                        else:
                            amplified_intervals_end = (
                                amplified_intervals_start[
                                    sorted_chrs[chri + 1]
                                ][0]
                                - margin_between_intervals
                            )
                        xtickpos.append(
                            (
                                amplified_intervals_start[chrom][int_idx]
                                + amplified_intervals_end
                            )
                            * 0.5,
                        )
        xtickpos.append(
            100 + self.num_amplified_intervals * margin_between_intervals
        )
        xticklabels = []
        for chrom in sorted_chrs:
            nint_chr = len(self.cycle_amplified_intervals[chrom])
            for int_idx in range(nint_chr):
                intv = self.cycle_amplified_intervals[chrom][int_idx]
                xticklabels.append(f"{intv.start:}   ")
                if nint_chr % 2 == 1 and int_idx == (nint_chr - 1) // 2:
                    xticklabels.append(chrom)
                xticklabels.append(f"{intv.end:}   ")
                if nint_chr % 2 == 0 and int_idx == (nint_chr - 2) // 2:
                    xticklabels.append(chrom)
        ax3.set_xticks(xtickpos)
        ax3.set_xticklabels(xticklabels, size=fontsize)
        ticks_labels = ax3.get_xticklabels()
        for ti in range(len(xticklabels)):
            if xticklabels[ti] not in sorted_chrs:
                ticks_labels[ti].set_rotation(90)
            else:
                ax3.xaxis.get_major_ticks()[ti].tick1line.set_visible(False)

        ax.set_yticks(cycleticks)
        ax.set_yticklabels(cycleticklabels, fontsize=fontsize)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        plt.tight_layout()
        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
        ax3.set_ylim(0, 1)
        fig.subplots_adjust(hspace=0)
        plt.savefig(output_fn + ".png", dpi=dpi)
        plt.savefig(output_fn + ".pdf")


@core_utils.profile_fn_with_call_counter
def plot_amplicon(
    ref: core_types.ReferenceGenome,
    bam_path: pathlib.Path | None,
    graph_file: io.TextIOWrapper | None,
    cycle_file: io.TextIOWrapper | None,
    output_prefix: str,
    num_cycles: int | None,
    max_coverage: float,
    min_mapq: float,
    gene_subset_list: list[str],
    gene_fontsize: int,
    region: str | None,
    *,
    should_plot_graph: bool,
    should_plot_cycles: bool,
    should_hide_genes: bool,
    should_restrict_to_bushman_genes: bool,
    should_plot_only_cyclic_walks: bool,
) -> None:
    if should_plot_graph:
        if not graph_file:
            print("Please specify the breakpoint graph file to plot.")
            sys.exit(1)
        if not bam_path:
            print("Please specify the bam file to plot.")
            sys.exit(1)

    if should_plot_cycles and not cycle_file:
        print("Please specify the cycle file, in *.bed format, to plot.")
        sys.exit(1)

    g = GraphViz()
    g.parse_genes(ref, set(gene_subset_list), should_restrict_to_bushman_genes)
    if should_plot_graph:
        bp_graph = parse_breakpoint_graph(graph_file)  # type: ignore[arg-type]
        g.open_bam(bam_path)
        g.graph = bp_graph
        if region:
            pchrom = region.split(":")[0]
            pb1, pb2 = region.split(":")[1].rsplit("-")
            g.plot_bounds = (pchrom, int(pb1), int(pb2))
        g.update_graph_intervals()
        gtitle = output_prefix
        if "/" in output_prefix:
            gtitle = output_prefix.split("/")[-1]
        g.plot_graph(
            gtitle,
            output_prefix + "_graph",
            max_cov_cutoff=max_coverage,
            quality_threshold=min_mapq,
            hide_genes=should_hide_genes,
            gene_font_size=gene_fontsize,
        )

    if should_plot_cycles:
        recon_cycles = parse_cycle_file(cycle_file, output_prefix, num_cycles)
        for cycle_id, cycle in recon_cycles.items():
            g.cycles[cycle_id] = cycle
        cycle_ids_ = None
        if num_cycles:
            cycle_ids_ = [i + 1 for i in range(num_cycles)]

        graph_given_ = graph_file is not None
        if graph_given_:
            g.graph = bp_graph
            g.update_graph_intervals()
        g.update_cycle_amplified_intervals(
            cycle_ids=cycle_ids_,
            cycle_only=should_plot_only_cyclic_walks,
            graph_given=graph_given_,
        )
        gtitle = output_prefix
        if "/" in output_prefix:
            gtitle = output_prefix.split("/")[-1]
        g.plot_cycles(
            gtitle,
            output_prefix + "_cycles",
            num_cycles=num_cycles,
            cycle_only=should_plot_only_cyclic_walks,
            hide_genes=should_hide_genes,
            gene_font_size=gene_fontsize,
        )
    g.close_bam()
    if graph_file:
        print(
            f"Visualization completed for {colorama.Fore.LIGHTCYAN_EX}"
            f"{graph_file.name} ({cycle_file.name if cycle_file else 'no cycles'}"  # type: ignore[union-attr]
            f"{colorama.Style.RESET_ALL})"
        )
    elif cycle_file:
        print(
            f"Visualization completed for {colorama.Fore.LIGHTCYAN_EX}{cycle_file.name}"
        )


def draw_cycle(
    g: BreakpointGraph,
    cycle: list[int],
    cycle_id: int,
    cycle_only: bool,
    graph_given: bool,
) -> None:
    pass
