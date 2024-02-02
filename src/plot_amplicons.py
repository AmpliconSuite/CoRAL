import os
import argparse
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc
from pylab import rcParams
rcParams['pdf.fonttype'] = 42

import global_names
from breakpoint_utilities import *


class graph_vis:
	"""
	Class for visualizing breakpoint graph and cycles 
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


	def __init__(self, bam_fn):
		self.lr_bamfh = pysam.AlignmentFile(bam_fn, "rb")

    
	def parse_graph_file(self, graph_fn):
		with open(graph_fn, 'r') as fp:
			for line in fp:
				s = line.strip().split("\t")
				if s[0] == "sequence":
					schr = s[1].split(":")[0]
					start = int(s[1].split(":")[1][:-1])
					end = int(s[2].split(":")[1][:-1])
					try:
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


	def plot_graph(self, title, output_fn, margin_between_intervals = 2, height = 7.5, fontsize = 18, dpi = 150):
		"""
		Plot discordant edges and coverage on sequence edges in breakpoint graph
		"""
		width = max(15, 2 * self.num_amplified_intervals)
		fig, ax = plt.subplots(figsize = (width, height))
		ax.set_title(title, fontsize = fontsize)
		ax2 = ax.twinx()

		# Draw sequence edges
		total_len_amp = 0 # Total length of amplified intervals
		for chr in self.amplified_intervals_from_graph.keys():
			total_len_amp += sum([int_[1] - int_[0] + 1 for int_ in self.amplified_intervals_from_graph[chr]])
		sorted_chrs = sorted(self.amplified_intervals_from_graph.keys(), key = lambda chr: global_names.chr_idx[chr])
		amplified_intervals_start = dict()
		ymax = 0
		x = margin_between_intervals
		for chr in sorted_chrs:
			interval_idx = 0
			amplified_intervals_start[chr] = [x]
			for seq in self.sequence_edges_by_chr[chr]:
				if seq[1] > self.amplified_intervals_from_graph[chr][interval_idx][1]:
					int_ = self.amplified_intervals_from_graph[chr][interval_idx]
					x += margin_between_intervals
					amplified_intervals_start[chr].append(x)
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
		for chr in amplified_intervals_start.keys():
			if chr != sorted_chrs[0]:
				plt.axvline(x = amplified_intervals_start[chr][0] - margin_between_intervals * 0.5, linestyle = "--", lw = 2, zorder = 2)
			for i in range(1, len(amplified_intervals_start[chr])):
				plt.axvline(x = amplified_intervals_start[chr][i] - margin_between_intervals * 0.5, linestyle = ":", lw = 2, zorder = 2)
		
		# Draw discordant edges
		colorcode = {"+-": "red", "++": "magenta", "-+": "brown", "--": "teal", "interchromosomal": "blue"}
		avg_bp_rc = sum([bp[7] for bp in self.discordant_edges]) * 1.0 / len(self.discordant_edges)
		for bp in self.discordant_edges:
			chr1 = bp[0]
			pos1 = bp[1]
			chr2 = bp[3]
			pos2 = bp[4]
			int1 = 0
			while pos1 > self.amplified_intervals_from_graph[chr1][int1][1]:
				int1 += 1
			bp_x1 = amplified_intervals_start[chr1][int1] + (pos1 - self.amplified_intervals_from_graph[chr1][int1][0]) * 100.0 / total_len_amp
			int2 = 0
			while pos2 > self.amplified_intervals_from_graph[chr2][int2][1]:
				int2 += 1
			bp_x2 = amplified_intervals_start[chr2][int2] + (pos2 - self.amplified_intervals_from_graph[chr2][int2][0]) * 100.0 / total_len_amp
			ort = bp[2] + bp[5]
			if chr1 != chr2:
				ort = "interchromosomal"
			arc = Arc(((bp_x1 + bp_x2) * 0.5, 0), bp_x1 - bp_x2, 2 * ymax, theta1 = 0, theta2 = 180, \
					color = colorcode[ort], lw = min(3 * (bp[7] / avg_bp_rc), 3), zorder = 3)
			ax2.add_patch(arc)
		ax2.set_ylim(0, 1.4 * ymax)
		ax2.set_ylabel('CN', fontsize = fontsize)
		ax2.tick_params(axis = 'y', labelsize = fontsize)

		# Draw coverage within amplified intervals
		max_cov = 0
		for chr in sorted_chrs:
			for inti in range(len(self.amplified_intervals_from_graph[chr])):
				int_ = self.amplified_intervals_from_graph[chr][inti]
				window_size = 150
				if int_[1] - int_[0] >= 1000000:
					window_size = 10000
				elif int_[1] - int_[0] >= 100000:
					window_size = 1000
				for w in range(int_[0], int_[1], window_size):
					cov = sum([sum(nc) for nc in self.lr_bamfh.count_coverage(chr, w, w + window_size, 
							quality_threshold = 0, read_callback = 'nofilter')]) * 1.0 / window_size
					if cov > max_cov:
						max_cov = cov
					x = amplified_intervals_start[chr][inti] + (w - int_[0]) * 100.0 / total_len_amp
					rect = Rectangle((x, 0), window_size * 100.0 / total_len_amp, cov, color ='silver', zorder = 1)
					ax.add_patch(rect)
				w = int_[1] - ((int_[1] - int_[0] + 1) % window_size)
				if w < int_[1]:
					cov = sum([sum(nc) for nc in self.lr_bamfh.count_coverage(chr, w, w + window_size, 
							quality_threshold = 0, read_callback = 'nofilter')]) * 1.0 / window_size
					if cov > max_cov:
						max_cov = cov
					x = amplified_intervals_start[chr][inti] + (w - int_[0]) * 100.0 / total_len_amp
					rect = Rectangle((x, 0), window_size * 100.0 / total_len_amp, cov, color ='silver', zorder = 1)
					ax.add_patch(rect)
		ax.set_ylabel('Coverage', fontsize = fontsize)
		ax.set_ylim(0, 1.4 * max_cov)
		ax.tick_params(axis = 'y', labelsize = fontsize)

		# Ticks an labels
		ax.set_xlim(0, 100 + (self.num_amplified_intervals + 1) * margin_between_intervals)
		ax2.set_xlim(0, 100 + (self.num_amplified_intervals + 1) * margin_between_intervals)
		xtickpos = []
		for chr in sorted_chrs:
			nint_chr = len(self.amplified_intervals_from_graph[chr])
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
			nint_chr = len(self.amplified_intervals_from_graph[chr])
			for inti in range(nint_chr):
				int_ = self.amplified_intervals_from_graph[chr][inti]
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
		
		plt.tight_layout()
		plt.savefig(output_fn, dpi = dpi)


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
		plt.savefig(output_fn, dpi = dpi)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Long read only amplicon reconstruction pipeline.")
	parser.add_argument("--lr_bam", help = "Sorted indexed (long read) bam file.", required = True)
	parser.add_argument("--graph", help = "Long read *.graph file.", required = True)
	parser.add_argument("--cycle", help = "Long read cycles file, in *.bed format.", required = True)
	parser.add_argument("--output_prefix", help = "Prefix of output files.", required = True)
	parser.add_argument("--plot_graph", help = "Visualize breakpoint graph.",  action = 'store_true')
	parser.add_argument("--plot_cycles", help = "Visualize (selected) cycles.",  action = 'store_true')
	parser.add_argument("--cycles_only", help = "Only plot cycles.",  action = 'store_true')
	parser.add_argument("--num_cycles", help = "Only plot the first NUM_CYCLES cycles.",  type = int)
	args = parser.parse_args()

	if args.plot_graph and args.graph == None:
		print ("Please specify the breakpoint graph file to plot.")
		os.abort()
	if args.plot_cycles and args.cycle == None:
		print ("Please specify the cycle file, in *.bed format, to plot.")
		os.abort()

	g = graph_vis(args.lr_bam)
	g.parse_graph_file(args.graph)
	g.parse_cycle_file(args.cycle)
	g.graph_amplified_intervals()
	if args.plot_graph:
		gtitle = args.output_prefix
		if '/' in args.output_prefix:
			gtitle = args.output_prefix.split('/')[-1]
		g.plot_graph(gtitle, args.output_prefix + "_graph.png")
	if args.plot_cycles:
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
			g.plotcycle(gtitle, args.output_prefix + "_cycles.png", num_cycles = args.num_cycles, cycle_only = cycle_only_)
		else:
			g.plotcycle(gtitle, args.output_prefix + "_cycles.png", cycle_only = cycle_only_)
	g.close_bam()
	print ("Visualization completed.")


