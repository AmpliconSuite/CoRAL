#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os

from intervaltree import IntervalTree
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc
import matplotlib.ticker as ticker
from pylab import rcParams
import pysam

import global_names
from breakpoint_utilities import *

rcParams['pdf.fonttype'] = 42


# makes a gene object from parsed refGene data
# this stores global properties for the gene
class gene(object):
    def __init__(self, gchrom, gstart, gend, gdata, highlight_name):
        self.gchrom = gchrom
        self.gstart = gstart
        self.gend = gend
        self.gname = gdata[-4]
        self.strand = gdata[3]
        self.height = 0.5
        self.highlight_name = highlight_name
        estarts = [int(x) for x in gdata[9].rsplit(",") if x]
        eends = [int(x) for x in gdata[10].rsplit(",") if x]
        self.eposns = list(zip(estarts, eends))

    def __str__(self):
        return f"Gene Name: {self.gname}, Chromosome: {self.gchrom}, Start: {self.gstart}, End: {self.gend}, Strand: {self.strand}"


class graph_vis:
    """
    """
    lr_bamfh = ""
    maxCN = 0.0
    sequence_edges_by_chr = dict()
    amplified_intervals_from_graph = dict()
    num_amplified_intervals = 0
    amplified_intervals_from_cycle = dict()
    discordant_edges = []
    cycles = dict()
    cycle_flags = dict()
    genes = defaultdict(IntervalTree)

    def __init__(self, bam_fn):
        self.lr_bamfh = pysam.AlignmentFile(bam_fn, "rb")

    def parse_genes(self, ref, gene_highlight_list, restrict_to_bushman=False):
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        if ref == "GRCh37" or ref == "hg19":
            refGene_name = "refGene_hg19.txt"
        elif ref == "GRCm38" or ref == "mm10":
            refGene_name = "refGene_mm10.txt"
        else:
            refGene_name = "refGene_" + ref + ".txt"

        seenNames = set()
        bushman_set = set()
        if restrict_to_bushman:
            with open(os.path.join(__location__, "annotations", "Bushman_group_allOnco_May2018.tsv")) as infile:
                _ = next(infile)
                for line in infile:
                    fields = line.rstrip().rsplit()
                    if not fields:
                        continue
                    bushman_set.add(fields[-1].strip("\""))

        with open(os.path.join(__location__, "annotations", refGene_name)) as infile:
            for line in infile:
                fields = line.rsplit("\t")
                currChrom = fields[2]
                if (ref == "GRCh37" or ref == "GRCm38") and not currChrom.startswith("hpv"):
                    currChrom = currChrom[3:]

                tstart = int(fields[4])
                tend = int(fields[5])
                gname = fields[-4]
                is_other_feature = (gname.startswith("LOC") or gname.startswith("LINC") or gname.startswith("MIR"))
                if restrict_to_bushman and gname not in bushman_set:
                    continue

                if gname not in seenNames and not is_other_feature:
                    seenNames.add(gname)
                    currGene = gene(currChrom, tstart, tend, fields, gname in gene_highlight_list)
                    self.genes[currChrom][tstart:tend] = currGene


    def parse_graph_file(self, graph_fn):
        with open(graph_fn, 'r') as fp:
            for line in fp:
                s = line.strip().split("\t")
                if s[0] == "sequence":
                    schr = s[1].split(":")[0]
                    start = int(s[1].split(":")[1][:-1])
                    end = int(s[2].split(":")[1][:-1])
                    try:
                        #altered to use AA graph (not CoRAL)
                        self.sequence_edges_by_chr[schr].append([schr, start, end, float(s[3]), int(s[6]), int(s[5])])
                    except:
                        self.sequence_edges_by_chr[schr] = [[schr, start, end, float(s[3]), int(s[6]), int(s[5])]]
                    if float(s[3]) > self.maxCN:
                        self.maxCN = float(s[3])
                elif s[0] == "discordant":
                    b1 = s[1].split("->")[0]
                    b2 = s[1].split("->")[1]
                    chr1 = b1.split(":")[0]
                    pos1 = int(b1.split(":")[1][:-1])
                    o1 = b1.split(":")[1][-1]
                    chr2 = b2.split(":")[0]
                    pos2 = int(b2.split(":")[1][:-1])
                    o2 = b2.split(":")[1][-1]
                    self.discordant_edges.append([chr1, pos1, o1, chr2, pos2, o2, float(s[2]), int(s[3])])


    def graph_amplified_intervals(self):
        for chr in self.sequence_edges_by_chr.keys():
            lstart, lend = -2, -2
            if chr not in self.amplified_intervals_from_graph:
                self.amplified_intervals_from_graph[chr] = []
            for se in self.sequence_edges_by_chr[chr]:
                start = se[1]
                end = se[2]
                if start != lend + 1:
                    if lstart >= 0:
                        self.amplified_intervals_from_graph[chr].append([lstart, lend])
                        self.num_amplified_intervals += 1
                    lstart = start
                    lend = end
                else:
                    lend = end
            self.amplified_intervals_from_graph[chr].append([lstart, lend])
            self.num_amplified_intervals += 1


    def parse_cycle_file(self, cycle_fn):
        with open(cycle_fn, 'r') as fp:
            for line in fp:
                s = line.strip().split("\t")
                if s[0][0] == '#':
                    continue
                if s[4] not in self.cycles:
                    self.cycles[s[4]] = [[s[0], int(s[1]), int(s[2]), s[3]]]
                    if s[5] == "True":
                        self.cycle_flags[s[4]] = [True, float(s[6])]
                    else:
                        self.cycle_flags[s[4]] = [False, float(s[6])]
                else:
                    self.cycles[s[4]].append([s[0], int(s[1]), int(s[2]), s[3]])


    def cycle_amplified_intervals(self, cycle_ids = None, cycle_only = False):
        """
        Derive amplified intervals from (selected) cycles
        """
        self.num_amplified_intervals = 0
        if cycle_ids == None:
            cycle_ids = [cycle_id for cycle_id in self.cycle_flags.keys()]
        if cycle_only:
            cycle_ids = [cycle_id for cycle_id in self.cycle_flags.keys() if self.cycle_flags[cycle_id][0]]
        for cycle_id in cycle_ids:
            for segment in self.cycles[cycle_id]:
                for int_ in self.amplified_intervals_from_graph[segment[0]]:
                    if interval_include(segment, [segment[0], int_[0], int_[1]]):
                        if segment[0] not in self.amplified_intervals_from_cycle:
                            self.amplified_intervals_from_cycle[segment[0]] = []
                        if int_ not in self.amplified_intervals_from_cycle[segment[0]]:
                            self.amplified_intervals_from_cycle[segment[0]].append(int_)
                        break
        for chr in self.amplified_intervals_from_cycle.keys():
            self.amplified_intervals_from_cycle[chr] = sorted(self.amplified_intervals_from_cycle[chr])
            self.num_amplified_intervals += len(self.amplified_intervals_from_cycle[chr])


    def set_gene_heights(self, rel_genes, padding=0):
        if not rel_genes:
            return

        gname_to_gobj = {x.gname: x for x in rel_genes}
        # merge intervals
        intervals = [(x.gstart, x.gend) for x in rel_genes]
        sorted_intervals = sorted(intervals)
        merged = [sorted_intervals[0]]
        for current in sorted_intervals[1:]:
            prev = merged[-1]
            if current[0] <= prev[1] + padding:
                merged[-1] = (prev[0], max(prev[1], current[1]))
            else:
                merged.append(current)

        gene_ival_t = IntervalTree()
        for x in rel_genes:
            gene_ival_t.addi(x.gstart, x.gend, x.gname)

        for mi in merged:
            ghits = gene_ival_t[mi[0]:mi[1]]
            gene_heights = np.linspace(0.15, 0.75, len(ghits))
            for g, h in zip(ghits, gene_heights):
                gname_to_gobj[g.data].height = h


    def plot_graph(self, title, output_fn, margin_between_intervals = 2, height = 7.5, fontsize = 18, dpi = 300,
                   max_cov_cutoff = float('inf'), quality_threshold=0, hide_genes=False, gene_font_size=12):
        """
        Plot discordant edges and coverage on sequence edges in breakpoint graph
        """
        width = max(15, 2 * self.num_amplified_intervals)
        fig = plt.figure(figsize = (width, height))
        if not hide_genes:
            gs = gridspec.GridSpec(2, 1, height_ratios=[8,2])
        else:
            gs = gridspec.GridSpec(2, 1, height_ratios=[8, 0])
        ax = fig.add_subplot(gs[0,0])
        plt.subplots_adjust(left=73/1000.0, right=1-73/1000.0, bottom=1/4.0, top=1-1/20.0)
        ax.set_title(title, fontsize = fontsize)
        ax2 = ax.twinx()
        ax3 = fig.add_subplot(gs[1, 0], sharex=ax)
        # ax.yaxis.set_label_coords(-0.05, 0.25)
        # ax2.yaxis.set_label_coords(1.05, 0.33)
        ax.xaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['top'].set_visible(False)

        # Draw sequence edges
        total_len_amp = 0 # Total length of amplified intervals
        for chrom in self.amplified_intervals_from_graph.keys():
            total_len_amp += sum([int_[1] - int_[0] + 1 for int_ in self.amplified_intervals_from_graph[chrom]])
        sorted_chrs = sorted(self.amplified_intervals_from_graph.keys(), key = lambda chr: global_names.chr_idx[chr])
        amplified_intervals_start = dict()
        ymax = 0
        x = margin_between_intervals
        for chrom in sorted_chrs:
            interval_idx = 0
            amplified_intervals_start[chrom] = [x]
            for seq in self.sequence_edges_by_chr[chrom]:
                if seq[1] > self.amplified_intervals_from_graph[chrom][interval_idx][1]:
                    int_ = self.amplified_intervals_from_graph[chrom][interval_idx]
                    x += margin_between_intervals
                    amplified_intervals_start[chrom].append(x)
                    interval_idx += 1
                x1 = x
                x += (seq[2] - seq[1]) * 100.0 / total_len_amp
                x2 = x
                y = seq[3]
                if y > ymax:
                    ymax = y
                ax2.hlines(y, x1, x2, color = "black", lw = 6, zorder = 2)
            x += margin_between_intervals

        # Draw amplified interval separators
        for chrom in amplified_intervals_start.keys():
            if chrom != sorted_chrs[0]:
                ax.axvline(x =amplified_intervals_start[chrom][0] - margin_between_intervals * 0.5, linestyle ="--", lw = 2, zorder = 2)
                ax3.axvline(x =amplified_intervals_start[chrom][0] - margin_between_intervals * 0.5, linestyle ="--", lw = 2, zorder = 2)
            for i in range(1, len(amplified_intervals_start[chrom])):
                ax.axvline(x =amplified_intervals_start[chrom][i] - margin_between_intervals * 0.5, linestyle =":", lw = 2, zorder = 2)

        # Draw discordant edges
        colorcode = {"+-": "red", "++": "magenta", "-+": (139/256.0, 69/256.0, 19/256.0), "--": "teal", "interchromosomal": "blue"}
        avg_bp_rc = sum([bp[7] for bp in self.discordant_edges]) * 1.0 / len(self.discordant_edges)
        for bp in self.discordant_edges:
            chr1 = bp[0]
            pos1 = bp[1]
            chr2 = bp[3]
            pos2 = bp[4]
            int1 = 0
            try:
                while pos1 > self.amplified_intervals_from_graph[chr1][int1][1]:
                    int1 += 1
                bp_x1 = amplified_intervals_start[chr1][int1] + (pos1 - self.amplified_intervals_from_graph[chr1][int1][0]) * 100.0 / total_len_amp
                int2 = 0
                while pos2 > self.amplified_intervals_from_graph[chr2][int2][1]:
                    int2 += 1
            except:
                print(str(bp), "not placeable")
                continue

            bp_x2 = amplified_intervals_start[chr2][int2] + (pos2 - self.amplified_intervals_from_graph[chr2][int2][0]) * 100.0 / total_len_amp
            ort = bp[2] + bp[5]
            if chr1 != chr2:
                ort = "interchromosomal"
            arc = Arc(((bp_x1 + bp_x2) * 0.5, 0), bp_x1 - bp_x2, 2 * ymax, theta1 = 0, theta2 = 180,
                    color = colorcode[ort], lw = min(3 * (bp[7] / avg_bp_rc), 3), zorder = 3)
            ax2.add_patch(arc)

        ax2.set_ylim(0, 1.4 * ymax)
        ax2.set_ylabel('CN', fontsize = fontsize)
        ax2.tick_params(axis = 'y', labelsize = fontsize)

        # Draw coverage within amplified intervals
        max_cov = 0
        for chrom in sorted_chrs:
            for inti in range(len(self.amplified_intervals_from_graph[chrom])):
                int_ = self.amplified_intervals_from_graph[chrom][inti]
                window_size = 150
                if int_[1] - int_[0] >= 1000000:
                    window_size = 10000
                elif int_[1] - int_[0] >= 100000:
                    window_size = 1000
                for w in range(int_[0], int_[1], window_size):
                    cov = sum([sum(nc) for nc in self.lr_bamfh.count_coverage(chrom, w, w + window_size,
                                                                              quality_threshold = quality_threshold, read_callback = 'nofilter')]) * 1.0 / window_size
                    if cov > max_cov:
                        max_cov = cov
                    x = amplified_intervals_start[chrom][inti] + (w - int_[0]) * 100.0 / total_len_amp
                    rect = Rectangle((x, 0), window_size * 100.0 / total_len_amp, cov, color ='silver', zorder = 1)
                    ax.add_patch(rect)
                w = int_[1] - ((int_[1] - int_[0] + 1) % window_size)
                if w < int_[1]:
                    cov = sum([sum(nc) for nc in self.lr_bamfh.count_coverage(chrom, w, w + window_size,
                                                                              quality_threshold = quality_threshold, read_callback = 'nofilter')]) * 1.0 / window_size
                    if cov > max_cov:
                        max_cov = cov
                    x = amplified_intervals_start[chrom][inti] + (w - int_[0]) * 100.0 / total_len_amp
                    rect = Rectangle((x, 0), window_size * 100.0 / total_len_amp, cov, color ='silver', zorder = 1)
                    ax.add_patch(rect)
        ax.set_ylabel('Coverage', fontsize = fontsize)
        ax.set_ylim(0, min(1.4 * max_cov, max_cov_cutoff))
        ax.tick_params(axis = 'y', labelsize = fontsize)

        # draw genes below plot
        if not hide_genes:
            for chrom in sorted_chrs:
                for inti in range(len(self.amplified_intervals_from_graph[chrom])):
                    int_ = self.amplified_intervals_from_graph[chrom][inti]
                    rel_genes = [x.data for x in self.genes[chrom][int_[0]:int_[1]]]
                    gene_padding = total_len_amp * 0.02
                    self.set_gene_heights(rel_genes, gene_padding)

                    for gene_obj in rel_genes:
                        height = gene_obj.height
                        print(gene_obj)
                        gene_start = amplified_intervals_start[chrom][inti] + (
                                    gene_obj.gstart - int_[0]) * 100.0 / total_len_amp
                        gene_end = amplified_intervals_start[chrom][inti] + (
                                    gene_obj.gend - int_[0]) * 100.0 / total_len_amp
                        ax3.hlines(height, gene_start, gene_end, color='cornflowerblue', lw=4.5)  # Draw horizontal bars for genes

                        ax3.text((gene_start + gene_end) / 2, height + 0.05, gene_obj.gname, ha='center', va='bottom',
                                 fontsize=gene_font_size)

                        if gene_obj.strand == '+':
                            ax3.plot(gene_start, height, marker='>', color='black', markersize=7)
                        elif gene_obj.strand == '-':
                            ax3.plot(gene_end, height, marker='<', color='black', markersize=7)

                        for exon_start, exon_end in gene_obj.eposns:
                            exon_start_pos = amplified_intervals_start[chrom][inti] + (
                                        exon_start - int_[0]) * 100.0 / total_len_amp
                            exon_end_pos = amplified_intervals_start[chrom][inti] + (
                                        exon_end - int_[0]) * 100.0 / total_len_amp

                            exon_min_width = 0.2  # Adjust the minimum width as needed
                            exon_width = exon_end_pos - exon_start_pos
                            if exon_width < exon_min_width:
                                diff = (exon_min_width - exon_width) / 2
                                exon_start_pos -= diff
                                exon_end_pos += diff

                            ax3.hlines(height, exon_start_pos, exon_end_pos, color='black', lw=7.5)

        # Ticks and labels
        ax.set_xlim(0, 100 + (self.num_amplified_intervals + 1) * margin_between_intervals)
        ax2.set_xlim(0, 100 + (self.num_amplified_intervals + 1) * margin_between_intervals)
        ax3.set_xlim(0, 100 + (self.num_amplified_intervals + 1) * margin_between_intervals)
        xtickpos = []
        for chrom in sorted_chrs:
            nint_chr = len(self.amplified_intervals_from_graph[chrom])
            for inti in range(len(amplified_intervals_start[chrom])):
                if inti > 0:
                    xtickpos.append(amplified_intervals_start[chrom][inti] - margin_between_intervals)
                    if nint_chr % 2 == 0 and inti == (nint_chr - 2) // 2 + 1:
                        xtickpos.append(amplified_intervals_start[chrom][inti] - margin_between_intervals * 0.5)
                    xtickpos.append(amplified_intervals_start[chrom][inti])
                    if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                        xtickpos.append((amplified_intervals_start[chrom][inti] + amplified_intervals_start[chrom][inti + 1] -
                                         margin_between_intervals) * 0.5)
                else:
                    if chrom != sorted_chrs[0]:
                        xtickpos.append(amplified_intervals_start[chrom][0] - margin_between_intervals)
                    xtickpos.append(amplified_intervals_start[chrom][0])
                    if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                        chri = sorted_chrs.index(chrom)
                        if chri == len(sorted_chrs) - 1:
                            amplified_intervals_end = 100 + self.num_amplified_intervals * margin_between_intervals
                        else:
                            amplified_intervals_end = amplified_intervals_start[sorted_chrs[chri + 1]][0] - margin_between_intervals
                        xtickpos.append((amplified_intervals_start[chrom][inti] + amplified_intervals_end) * 0.5)
        xtickpos.append(100 + self.num_amplified_intervals * margin_between_intervals)
        xticklabels = []
        for chrom in sorted_chrs:
            nint_chr = len(self.amplified_intervals_from_graph[chrom])
            for inti in range(nint_chr):
                int_ = self.amplified_intervals_from_graph[chrom][inti]
                xticklabels.append(str(int_[0]) + "   ")
                if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                    xticklabels.append(chrom)
                xticklabels.append(str(int_[1]) + "   ")
                if nint_chr % 2 == 0 and inti == (nint_chr - 2) // 2:
                    xticklabels.append(chrom)
        ax3.set_xticks(xtickpos)
        ax3.set_xticklabels(xticklabels, size = fontsize)
        ticks_labels = ax3.get_xticklabels()
        for ti in range(len(xticklabels)):
            if not xticklabels[ti] in sorted_chrs:
                ticks_labels[ti].set_rotation(90)
            else:
                ax3.xaxis.get_major_ticks()[ti].tick1line.set_visible(False)

        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
        ax3.set_ylim(0,1)
        fig.subplots_adjust(hspace=0)
        plt.savefig(output_fn + ".png", dpi = dpi)
        plt.savefig(output_fn + ".pdf")


    def close_bam(self):
        self.lr_bamfh.close()


    def plotcycle(self, title, output_fn, num_cycles = -1, cycle_only = False, margin_between_intervals = 2,
            fontsize = 18, dpi = 150):
        """
        Plot cycles & paths returned from cycle decomposition
        """
        width = max(15, 2 * self.num_amplified_intervals)
        cycles_to_plot = [cycle_id for cycle_id in self.cycles.keys()]
        if num_cycles > 0:
            cycles_to_plot = [cycle_id for cycle_id in cycles_to_plot if int(cycle_id) <= num_cycles]
        if cycle_only:
            cycles_to_plot = [cycle_id for cycle_id in cycles_to_plot if self.cycle_flags[cycle_id][0]]
        cycles_to_plot = sorted(cycles_to_plot)
        height = sum([2 * len(self.cycles[cycle_id]) - 1 for cycle_id in cycles_to_plot]) + 6 * (len(cycles_to_plot) - 1)
        fig, ax = plt.subplots(figsize = (width, max(4, height * 0.25)))
        ax.set_title(title, fontsize = fontsize)

        # Compute the x coordinates for each amplified interval
        total_len_amp = 0 # Total length of amplified intervals
        for chr in self.amplified_intervals_from_cycle.keys():
            total_len_amp += sum([int_[1] - int_[0] + 1 for int_ in self.amplified_intervals_from_cycle[chr]])
        sorted_chrs = sorted(self.amplified_intervals_from_cycle.keys(), key = lambda chr: global_names.chr_idx[chr])
        amplified_intervals_start = dict()
        x = margin_between_intervals
        for chr in sorted_chrs:
            amplified_intervals_start[chr] = [x]
            for interval_idx in range(len(self.amplified_intervals_from_cycle[chr])):
                int_ = self.amplified_intervals_from_cycle[chr][interval_idx]
                x += (int_[1] - int_[0]) * 100.0 / total_len_amp
                x += margin_between_intervals
                if interval_idx < len(self.amplified_intervals_from_cycle[chr]) - 1:
                    amplified_intervals_start[chr].append(x)

        # Draw amplified interval separators
        for chr in amplified_intervals_start.keys():
            if chr != sorted_chrs[0]:
                plt.axvline(x = amplified_intervals_start[chr][0] - margin_between_intervals * 0.5, linestyle = "--", lw = 2)
            for i in range(1, len(amplified_intervals_start[chr])):
                plt.axvline(x = amplified_intervals_start[chr][i] - margin_between_intervals * 0.5, linestyle = ":", lw = 2)

        # Draw cycles
        y_cur = -2
        extension = 1.5
        cycleticks = []
        cycleticklabels = []
        for cycle_id in cycles_to_plot:
            ystart_cycle_id = y_cur
            for i in range(len(self.cycles[cycle_id])):
                # Segment i
                seg = self.cycles[cycle_id][i]
                interval_idx = 0
                while seg[1] > self.amplified_intervals_from_cycle[seg[0]][interval_idx][1]:
                    interval_idx += 1
                x1 = amplified_intervals_start[seg[0]][interval_idx] + (seg[1] - self.amplified_intervals_from_cycle[seg[0]][interval_idx][0]) * 100.0 / total_len_amp
                xlen = (seg[2] - seg[1]) * 100.0 / total_len_amp
                rect = Rectangle((x1, y_cur), xlen, 1, facecolor='#FFFEF7', linewidth = 2, edgecolor = 'k')
                ax.add_patch(rect)

                # Connections between segment i and i + 1
                if i < len(self.cycles[cycle_id]) - 1:
                    nseg = self.cycles[cycle_id][i + 1]
                    interval_idx_n = 0
                    while nseg[1] > self.amplified_intervals_from_cycle[nseg[0]][interval_idx_n][1]:
                        interval_idx_n += 1
                    if seg[3] == '+' and nseg[3] == '-':
                        x2 = x1 + xlen
                        x2n = amplified_intervals_start[nseg[0]][interval_idx_n]
                        x2n += (nseg[2] - self.amplified_intervals_from_cycle[nseg[0]][interval_idx_n][0]) * 100.0 / total_len_amp
                        ax.vlines(x = max(x2, x2n) + extension, ymin = y_cur + 0.5, ymax = y_cur - 1.5, colors = 'b', lw = 2)
                        ax.hlines(y = y_cur + 0.5, xmin = x2, xmax = max(x2, x2n) + extension, colors = 'b', lw = 2)
                        ax.hlines(y = y_cur - 1.5, xmin = x2n, xmax = max(x2, x2n) + extension, colors = 'b', lw = 2)
                        y_cur -= 2
                    elif seg[3] == '-' and nseg[3] == '+':
                        x1n = amplified_intervals_start[nseg[0]][interval_idx_n]
                        x1n += (nseg[1] - self.amplified_intervals_from_cycle[nseg[0]][interval_idx_n][0]) * 100.0 / total_len_amp
                        ax.vlines(x = min(x1, x1n) - extension, ymin = y_cur + 0.5, ymax = y_cur - 1.5, colors = 'b', lw = 2)
                        ax.hlines(y = y_cur + 0.5, xmin = min(x1, x1n) - extension, xmax = x1, colors = 'b', lw = 2)
                        ax.hlines(y = y_cur - 1.5, xmin = min(x1, x1n) - extension, xmax = x1n, colors = 'b', lw = 2)
                        y_cur -= 2
                    elif seg[3] == '+' and nseg[3] == '+':
                        x2 = x1 + xlen
                        x1n = amplified_intervals_start[nseg[0]][interval_idx_n]
                        x1n += (nseg[1] - self.amplified_intervals_from_cycle[nseg[0]][interval_idx_n][0]) * 100.0 / total_len_amp
                        if x2 <= x1n:
                            ax.hlines(y = y_cur + 0.5, xmin = x2, xmax = x1n, colors = 'b', lw = 2)
                        else:
                            ax.vlines(x = x2 + extension, ymin = y_cur - 0.5, ymax = y_cur + 0.5, colors = 'b', lw = 2)
                            ax.vlines(x = x1n - extension, ymin = y_cur - 1.5, ymax = y_cur - 0.5, colors = 'b', lw = 2)
                            ax.hlines(y = y_cur + 0.5, xmin = x2, xmax = x2 + extension, colors = 'b', lw = 2)
                            ax.hlines(y = y_cur - 0.5, xmin = x1n - extension, xmax = x2 + extension, colors = 'b', lw = 2)
                            ax.hlines(y = y_cur - 1.5, xmin = x1n - extension, xmax = x1n, colors = 'b', lw = 2)
                            y_cur -= 2
                    else: # seg[3] == '-' and nseg[3] == '-'
                        x2n = amplified_intervals_start[nseg[0]][interval_idx_n]
                        x2n += (nseg[2] - self.amplified_intervals_from_cycle[nseg[0]][interval_idx_n][0]) * 100.0 / total_len_amp
                        if x1 >= x2n:
                            ax.hlines(y = y_cur + 0.5, xmin = x2n, xmax = x1, colors = 'b', lw = 2)
                        else:
                            ax.vlines(x = x1 - extension, ymin = y_cur - 0.5, ymax = y_cur + 0.5, colors = 'b', lw = 2)
                            ax.vlines(x = x2n + extension, ymin = y_cur - 1.5, ymax = y_cur - 0.5, colors = 'b', lw = 2)
                            ax.hlines(y = y_cur + 0.5, xmin = x1 - extension, xmax = x1, colors = 'b', lw = 2)
                            ax.hlines(y = y_cur - 0.5, xmin = x1 - extension, xmax = x2n + extension, colors = 'b', lw = 2)
                            ax.hlines(y = y_cur - 1.5, xmin = x2n, xmax = x2n + extension, colors = 'b', lw = 2)
                            y_cur -= 2

            # First and last segments
            if not self.cycle_flags[cycle_id][0]: # Paths
                seg = self.cycles[cycle_id][0]
                interval_idx = 0
                while seg[1] > self.amplified_intervals_from_cycle[seg[0]][interval_idx][1]:
                    interval_idx += 1
                if seg[3] == '+':
                    x1 = amplified_intervals_start[seg[0]][interval_idx] + \
                        (seg[1] - self.amplified_intervals_from_cycle[seg[0]][interval_idx][0]) * 100.0 / total_len_amp
                    ax.hlines(y = ystart_cycle_id + 0.5, xmin = x1 - 2 * extension, xmax = x1, colors = 'b', lw = 2)
                else:
                    x2 = amplified_intervals_start[seg[0]][interval_idx] + \
                        (seg[2] - self.amplified_intervals_from_cycle[seg[0]][interval_idx][0]) * 100.0 / total_len_amp
                    ax.hlines(y = ystart_cycle_id + 0.5, xmin = x2, xmax = x2 + 2 * extension, colors = 'b', lw = 2)
                seg = self.cycles[cycle_id][-1]
                interval_idx = 0
                while seg[1] > self.amplified_intervals_from_cycle[seg[0]][interval_idx][1]:
                    interval_idx += 1
                if seg[3] == '+':
                    x2 = amplified_intervals_start[seg[0]][interval_idx]
                    x2 += (seg[2] - self.amplified_intervals_from_cycle[seg[0]][interval_idx][0]) * 100.0 / total_len_amp
                    ax.hlines(y = y_cur + 0.5, xmin = x2, xmax = x2 + 2 * extension, colors = 'b', lw = 2)
                else:
                    x1 = amplified_intervals_start[seg[0]][interval_idx]
                    x1 += (seg[1] - self.amplified_intervals_from_cycle[seg[0]][interval_idx][0]) * 100.0 / total_len_amp
                    ax.hlines(y = y_cur + 0.5, xmin = x1 - 2 * extension, xmax = x1, colors = 'b', lw = 2)
            else: # Cycles
                xmin_ = 0.5
                xmax_ = 99.5 + (self.num_amplified_intervals + 1) * margin_between_intervals
                seg1 = self.cycles[cycle_id][0]
                interval_idx1 = 0
                while seg1[1] > self.amplified_intervals_from_cycle[seg1[0]][interval_idx1][1]:
                    interval_idx1 += 1
                seg2 = self.cycles[cycle_id][-1]
                interval_idx2 = 0
                while seg2[1] > self.amplified_intervals_from_cycle[seg2[0]][interval_idx2][1]:
                    interval_idx2 += 1
                if seg1[3] == '-' and seg2[3] == '+':
                    x2 = amplified_intervals_start[seg1[0]][interval_idx1]
                    x2 += (seg1[2] - self.amplified_intervals_from_cycle[seg1[0]][interval_idx1][0]) * 100.0 / total_len_amp
                    x2n = amplified_intervals_start[seg2[0]][interval_idx2]
                    x2n += (seg2[2] - self.amplified_intervals_from_cycle[seg2[0]][interval_idx2][0]) * 100.0 / total_len_amp
                    ax.vlines(x = xmax_, ymin = y_cur + 0.5, ymax = ystart_cycle_id + 0.5, colors = 'b', lw = 2)
                    ax.hlines(y = ystart_cycle_id + 0.5, xmin = x2, xmax = xmax_, colors = 'b', lw = 2)
                    ax.hlines(y = y_cur + 0.5, xmin = x2n, xmax = xmax_, colors = 'b', lw = 2)
                elif seg1[3] == '+' and seg2[3] == '-':
                    x1 = amplified_intervals_start[seg1[0]][interval_idx1]
                    x1 += (seg1[1] - self.amplified_intervals_from_cycle[seg1[0]][interval_idx1][0]) * 100.0 / total_len_amp
                    x1n = amplified_intervals_start[seg2[0]][interval_idx2]
                    x1n += (seg2[1] - self.amplified_intervals_from_cycle[seg2[0]][interval_idx2][0]) * 100.0 / total_len_amp
                    ax.vlines(x = xmin_, ymin = y_cur + 0.5, ymax = ystart_cycle_id + 0.5, colors = 'b', lw = 2)
                    ax.hlines(y = ystart_cycle_id + 0.5, xmin = xmin_, xmax = x1, colors = 'b', lw = 2)
                    ax.hlines(y = y_cur + 0.5, xmin = xmin_, xmax = x1n, colors = 'b', lw = 2)
                elif seg1[3] == '-' and seg2[3] == '-':
                    x2 = amplified_intervals_start[seg1[0]][interval_idx1]
                    x2 += (seg1[2] - self.amplified_intervals_from_cycle[seg1[0]][interval_idx1][0]) * 100.0 / total_len_amp
                    x1n = amplified_intervals_start[seg2[0]][interval_idx2]
                    x1n += (seg2[1] - self.amplified_intervals_from_cycle[seg2[0]][interval_idx2][0]) * 100.0 / total_len_amp
                    ax.vlines(x = xmax_, ymin = y_cur - 0.5, ymax = ystart_cycle_id + 0.5, colors = 'b', lw = 2)
                    ax.vlines(x = x1n - extension, ymin = y_cur - 0.5, ymax = y_cur + 0.5, colors = 'b', lw = 2)
                    ax.hlines(y = ystart_cycle_id + 0.5, xmin = x2, xmax = xmax_, colors = 'b', lw = 2)
                    ax.hlines(y = y_cur + 0.5, xmin = x1n - extension, xmax = x1n, colors = 'b', lw = 2)
                    ax.hlines(y = y_cur - 0.5, xmin = x1n - extension, xmax = xmax_, colors = 'b', lw = 2)
                else:
                    x1 = amplified_intervals_start[seg1[0]][interval_idx1]
                    x1 += (seg1[1] - self.amplified_intervals_from_cycle[seg1[0]][interval_idx1][0]) * 100.0 / total_len_amp
                    x2n = amplified_intervals_start[seg2[0]][interval_idx2]
                    x2n += (seg2[2] - self.amplified_intervals_from_cycle[seg2[0]][interval_idx2][0]) * 100.0 / total_len_amp
                    ax.vlines(x = xmin_, ymin = y_cur - 0.5, ymax = ystart_cycle_id + 0.5, colors = 'b', lw = 2)
                    ax.vlines(x = x2n + extension, ymin = y_cur - 0.5, ymax = y_cur + 0.5, colors = 'b', lw = 2)
                    ax.hlines(y = ystart_cycle_id + 0.5, xmin = xmin_, xmax = x1, colors = 'b', lw = 2)
                    ax.hlines(y = y_cur + 0.5, xmin = x2n, xmax = x2n + extension, colors = 'b', lw = 2)
                    ax.hlines(y = y_cur - 0.5, xmin = xmin_, xmax = x2n + extension, colors = 'b', lw = 2)

            # Separators between cycles; ticks
            ax.hlines(y = y_cur - 2, xmin = -1, xmax = 101 + (self.num_amplified_intervals + 1) * margin_between_intervals, colors = 'k')
            cycleticks.append((y_cur + ystart_cycle_id) * 0.5)
            if self.cycle_flags[cycle_id][0]:
                cycleticklabels.append("cycle " + cycle_id + ":\nCN = " + str(round(self.cycle_flags[cycle_id][1], 2)))
            else:
                cycleticklabels.append("path " + cycle_id + ":\nCN = " + str(round(self.cycle_flags[cycle_id][1], 2)))
            y_cur -= 4

        # Ticks an labels
        plt.xlim(-1, 101 + (self.num_amplified_intervals + 1) * margin_between_intervals)
        plt.ylim(y_cur + 2, 0)
        xtickpos = []
        for chr in sorted_chrs:
            nint_chr = len(self.amplified_intervals_from_cycle[chr])
            for inti in range(len(amplified_intervals_start[chr])):
                if inti > 0:
                    xtickpos.append(amplified_intervals_start[chr][inti] - margin_between_intervals)
                    if nint_chr % 2 == 0 and inti == (nint_chr - 2) // 2 + 1:
                        xtickpos.append(amplified_intervals_start[chr][inti] - margin_between_intervals * 0.5)
                    xtickpos.append(amplified_intervals_start[chr][inti])
                    if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                        xtickpos.append((amplified_intervals_start[chr][inti] + amplified_intervals_start[chr][inti + 1] -
                                margin_between_intervals) * 0.5)
                else:
                    if chr != sorted_chrs[0]:
                        xtickpos.append(amplified_intervals_start[chr][0] - margin_between_intervals)
                    xtickpos.append(amplified_intervals_start[chr][0])
                    if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                        chri = sorted_chrs.index(chr)
                        if chri == len(sorted_chrs) - 1:
                            amplified_intervals_end = 100 + self.num_amplified_intervals * margin_between_intervals
                        else:
                            amplified_intervals_end = amplified_intervals_start[sorted_chrs[chri + 1]][0] - margin_between_intervals
                        xtickpos.append((amplified_intervals_start[chr][inti] + amplified_intervals_end) * 0.5)
        xtickpos.append(100 + self.num_amplified_intervals * margin_between_intervals)
        xticklabels = []
        for chr in sorted_chrs:
            nint_chr = len(self.amplified_intervals_from_cycle[chr])
            for inti in range(nint_chr):
                int_ = self.amplified_intervals_from_cycle[chr][inti]
                xticklabels.append(str(int_[0]) + "   ")
                if nint_chr % 2 == 1 and inti == (nint_chr - 1) // 2:
                    xticklabels.append(chr)
                xticklabels.append(str(int_[1]) + "   ")
                if nint_chr % 2 == 0 and inti == (nint_chr - 2) // 2:
                    xticklabels.append(chr)
        ax.set_xticks(xtickpos)
        ax.set_xticklabels(xticklabels, size = fontsize)
        ticks_labels = ax.get_xticklabels()
        for ti in range(len(xticklabels)):
            if xticklabels[ti][:3] != 'chr':
                ticks_labels[ti].set_rotation(90)
        ax.set_yticks(cycleticks)
        ax.set_yticklabels(cycleticklabels, fontsize = fontsize)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.tight_layout()
        plt.savefig(output_fn + '.png', dpi = dpi)
        plt.savefig(output_fn + '.pdf')




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Long read only amplicon reconstruction pipeline.")
    parser.add_argument("--ref", help = "Name of reference genome used", choices=["hg19", "hg38", "mm10"], required=True)
    parser.add_argument("--bam", help = "Sorted & indexed bam file.", required = True)
    parser.add_argument("--graph", help = "AmpliconSuite-formatted *.graph file.", required=True)
    parser.add_argument("--cycle", help = "AmpliconSuite-formatted cycles file, in *.bed format.")
    parser.add_argument("--output_prefix", help = "Prefix of output files.", required = True)
    parser.add_argument("--plot_graph", help = "Visualize breakpoint graph.",  action = 'store_true')
    parser.add_argument("--plot_cycles", help = "Visualize (selected) cycles.",  action = 'store_true')
    parser.add_argument("--cycles_only", help = "Only plot cycles.",  action = 'store_true')
    parser.add_argument("--num_cycles", help = "Only plot the first NUM_CYCLES cycles.",  type = int)
    parser.add_argument("--max_coverage", help = "Limit the maximum visualized coverage in the graph", type = float, default = float('inf'))
    parser.add_argument("--min_mapq", help = "Minimum mapping quality to count read in coverage plotting", type = float, default = 0)
    parser.add_argument("--gene_subset_list", help = "List of genes to visualize (will show all by default)", nargs='+', default=[])
    parser.add_argument("--hide_genes", help="Do not show gene track", action='store_true', default=False)
    parser.add_argument("--gene_fontsize", help="Change size of gene font", type=float, default=12)
    parser.add_argument("--bushman_genes", help="Reduce gene set to the Bushman cancer-related gene set", action='store_true', default=False)
    args = parser.parse_args()

    if args.plot_graph and args.graph == None:
        print ("Please specify the breakpoint graph file to plot.")
        os.abort()
    if args.plot_cycles and args.cycle == None:
        print ("Please specify the cycle file, in *.bed format, to plot.")
        os.abort()

    g = graph_vis(args.bam)
    g.parse_genes(args.ref, set(args.gene_subset_list), args.bushman_genes)
    g.parse_graph_file(args.graph)
    if args.plot_graph:
        g.graph_amplified_intervals()
        gtitle = args.output_prefix
        if '/' in args.output_prefix:
            gtitle = args.output_prefix.split('/')[-1]
        g.plot_graph(gtitle, args.output_prefix + "_graph", max_cov_cutoff=args.max_coverage, quality_threshold=args.min_mapq,
                     hide_genes=args.hide_genes, gene_font_size=args.gene_fontsize)
    if args.plot_cycles:
        g.parse_cycle_file(args.cycle)
        cycle_ids_ = None
        cycle_only_ = False
        if args.num_cycles:
            cycle_ids_ = [str(i + 1) for i in range(args.num_cycles)]
        if args.cycles_only:
            cycle_only_ = True
        g.cycle_amplified_intervals(cycle_ids = cycle_ids_, cycle_only = cycle_only_)
        gtitle = args.output_prefix
        if '/' in args.output_prefix:
            gtitle = args.output_prefix.split('/')[-1]
        if args.num_cycles:
            g.plotcycle(gtitle, args.output_prefix + "_cycles", num_cycles = args.num_cycles, cycle_only = cycle_only_)
        else:
            g.plotcycle(gtitle, args.output_prefix + "_cycles", cycle_only = cycle_only_)
    g.close_bam()
    print ("Visualization completed.")


