import os
import sys

import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = [20, 8]
rcParams['pdf.fonttype'] = 42
mpl.rc('xtick', labelsize = 25) 
mpl.rc('ytick', labelsize = 25)

from breakpoint_utilities import *
import cigar_parsing
import cycle2bed
import global_names


def fetch(lr_bamfh):
	read_length = dict()
	chimeric_alignments = dict()
	for read in lr_bamfh.fetch():
		rn = read.query_name
		#if rn == "0d54d7a9-8c1b-45d9-ad33-d66bef4f63ff":
		#	print (read.query_sequence)
		if read.flag < 256 and rn not in read_length:
			read_length[rn] = read.query_length
		try:
			sa_list = read.get_tag('SA:Z:')[:-1].split(';')
			for sa in sa_list:
				try:
					if sa not in chimeric_alignments[rn]:
						chimeric_alignments[rn].append(sa) 
				except:
					chimeric_alignments[rn] = [sa]
		except:
			pass
	reads_wo_primary_alignment = []
	for r in chimeric_alignments.keys():
		if r not in read_length:
			reads_wo_primary_alignment.append(r)
			continue
		rl = read_length[r]
		chimeric_alignments[r] = cigar_parsing.alignment_from_satags(chimeric_alignments[r], rl)
	for r in reads_wo_primary_alignment:
		del chimeric_alignments[r]
	return read_length, chimeric_alignments


def locate_hsrs(args):
	chr_sizes = global_names.chr_sizes
	ecdna_intervals = []
	ecdna_intervals_ext = []

	cycle_fn = args.cycles
	if args.cycles.endswith("_cycles.txt"):
		# convert it to a bed
		init_char = "" if args.output_prefix.endswith("/") else "_"
		conv_cycle_fn = args.output_prefix + init_char + "converted_"
		conv_cycle_fn += "cycles.bed"
		cycle2bed.convert_cycles_to_bed(args.cycles, conv_cycle_fn)
		cycle_fn = conv_cycle_fn

	elif not args.cycles.endswith(".bed"):
		sys.stderr.write(args.cycles + "\n")
		sys.stderr.write("Cycles file must be either a valid *_cycles.txt file or a converted .bed file!\n")
		sys.exit(1)

	with open(cycle_fn, 'r') as fp:
		for line in fp:
			if line.startswith("#"):
				continue

			s = line.strip().split()
			ecdna_intervals.append([s[0], int(s[1]), int(s[2])])
			ecdna_intervals_ext.append([s[0], int(s[1]) - args.bp_match_cutoff, int(s[2]) + args.bp_match_cutoff])
	print("ecDNA intervals:")
	for ival in ecdna_intervals:
		print(ival)

	cns_dict = dict()
	with open(args.cn_seg, 'r') as fp:
		for line in fp:
			s = line.strip().split()
			if line.startswith('chromosome'):
				continue

			if s[0] not in cns_dict:
				cns_dict[s[0]] = dict()

			if args.cn_seg.endswith(".cns"):
				cn = 2 * (2 ** float(s[4]))
			elif args.cn_seg.endswith(".bed"):
				cn = float(s[3])
			else:
				sys.stderr.write(args.cn_seg + "\n")
				sys.stderr.write("Invalid cn_seg file format!\n")

			try:
				cns_dict[s[0]].append([int(s[1]), int(s[2]), cn])
			except:
				cns_dict[s[0]] = [[int(s[1]), int(s[2]), cn]]


	lr_bamfh = pysam.AlignmentFile(args.lr_bam, 'rb')
	read_length, chimeric_alignments = fetch(lr_bamfh)
	print ("Fetched %d chimeric alignments." %len(chimeric_alignments))

	bp_list = []
	for r in chimeric_alignments.keys():
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
				if q_[ri] >= 20 and q_[ri + 1] >= 20 and interval_overlap_l(rr_int[ri], ecdna_intervals) == -1 and \
					interval_overlap_l(rr_int[ri + 1], ecdna_intervals) >= 0:
					bp_list.append(interval2bp(rr_int[ri], rr_int[ri + 1], (r, ri, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri][1])) + [q_[ri], q_[ri + 1]])
					bassigned[ri] = 1
				elif q_[ri] >= 20 and q_[ri + 1] >= 20 and interval_overlap_l(rr_int[ri], ecdna_intervals) >= 0 and \
					interval_overlap_l(rr_int[ri + 1], ecdna_intervals) == -1:
					bp_list.append(interval2bp(rr_int[ri], rr_int[ri + 1], (r, ri, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri][1])) + [q_[ri], q_[ri + 1]])
					bassigned[ri] = 1
			"""
			Breakpoint from local alignment i - 1 and i + 1
			"""
			for ri in range(1, len(rr_int) - 1):
				if bassigned[ri - 1] == 0 and bassigned[ri] == 0 and \
					q_[ri] < 10 and q_[ri - 1] >= 20 and q_[ri + 1] >= 20 and \
					interval_overlap_l(rr_int[ri - 1], ecdna_intervals) == -1 and \
					interval_overlap(rr_int[ri + 1], ecdna_intervals) >= 0:
					bp_list.append(interval2bp(rr_int[ri - 1], rr_int[ri + 1], (r, ri - 1, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri - 1][1])) + [q_[ri - 1], q_[ri + 1]])
				elif bassigned[ri - 1] == 0 and bassigned[ri] == 0 and \
					q_[ri] < 10 and q_[ri - 1] >= 20 and q_[ri + 1] >= 20 and \
					interval_overlap_l(rr_int[ri - 1], ecdna_intervals) == -1 and \
					interval_overlap(rr_int[ri + 1], ecdna_intervals) >= 0:
					bp_list.append(interval2bp(rr_int[ri - 1], rr_int[ri + 1], (r, ri - 1, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri - 1][1])) + [q_[ri - 1], q_[ri + 1]])
	
	bp_clusters = cluster_bp_list(bp_list, float(args.normal_cov) * 0.5, args.bp_match_cutoff_clustering)
	bp_refined = []
	bp_stats = []
	for c in bp_clusters:
		if len(c) >= float(args.normal_cov) * 0.5:
			bp_cluster_r = c
			while len(bp_cluster_r) >= float(args.normal_cov) * 0.5:
				bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(bp_cluster_r, args.bp_match_cutoff)
				#print (bp[:6])
				if len(set(bpr)) >= float(args.normal_cov) * 0.5:
					bpi_ = -1
					for bpi in range(len(bp_refined)):
						bp_ = bp_refined[bpi]
						if bp[0] == bp_[0] and bp[3] == bp_[3] and bp[2] == bp_[2] and bp[5] == bp_[5] and \
							abs(bp[1] - bp_[1]) <= args.bp_match_cutoff and abs(bp[4] - bp_[4]) < args.bp_match_cutoff:
							bp_refined[bpi][-1] |= set(bpr)
							bpi_ = bpi
							break
					if bpi_ < 0:
						bp_refined.append(bp + [bpr])
						bp_stats.append(bp_stats_)
	print ("Found %d breakpoints connecting ecDNA and chromosomes." %len(bp_refined))
	lr_bamfh.close()
	sum_sizes = sum(chr_sizes.values())
	agg_size = 0
	xtick_pos = []
	starting_pos = dict()
	for chr in chr_sizes.keys():
		agg_size += chr_sizes[chr]
		if agg_size < sum_sizes:
			plt.plot([agg_size * 100.0 / sum_sizes, agg_size * 100.0 / sum_sizes], [-1, 1000000], 'k--', linewidth = 2)
		xtick_pos.append((agg_size - 0.5 * chr_sizes[chr]) * 100.0 / sum_sizes)
		starting_pos[chr] = (agg_size - chr_sizes[chr]) * 100.0 / sum_sizes

	for bp in bp_refined:
		if interval_overlap_l([bp[0], bp[1], bp[1]], ecdna_intervals_ext) >= 0 and interval_overlap_l([bp[3], bp[4], bp[4]], ecdna_intervals_ext) < 0:
			if bp[3] in starting_pos.keys():
				cn = 0.0
				for seg in cns_dict[bp[3]]:
					if bp[4] > seg[0] and bp[4] < seg[1]:
						cn = seg[2]
						break
				if cn <= 5.0 and len(bp[-1]) <= float(args.normal_cov) * 2.5:
					print ("Breakpoint", bp[:6], "Support = ", len(bp[-1]))
					xpos = starting_pos[bp[3]] + bp[4] * 100.0 / sum_sizes
					ypos = len(bp[-1])
					plt.plot(xpos, ypos, 'bo')
		elif interval_overlap_l([bp[0], bp[1], bp[1]], ecdna_intervals_ext) < 0 and interval_overlap_l([bp[3], bp[4], bp[4]], ecdna_intervals_ext) >= 0:
			if bp[0] in starting_pos.keys():
				cn = 0.0
				for seg in cns_dict[bp[0]]:
					if bp[1] > seg[0] and bp[1] < seg[1]:
						cn = seg[2]
						break
				if cn <= 5.0 and len(bp[-1]) <= float(args.normal_cov) * 2.5:
					print ("Breakpoint", bp[:6], "Support = ", len(bp[-1]))
					xpos = starting_pos[bp[0]] + bp[1] * 100.0 / sum_sizes
					ypos = len(bp[-1])
					plt.plot(xpos, ypos, 'bo')

	plt.xlim([0, 100])
	plt.ylim([1, 500])
	plt.yscale('log')
	plt.xticks(xtick_pos, list(range(1, 23)) + ['X', 'Y'])
	plt.title(args.output_prefix + " integration loci", fontsize = 25)
	plt.ylabel('Long read support', fontsize = 25)
	plt.tight_layout()
	out_img_name = "integration_sites_" + args.output_prefix
	plt.savefig(out_img_name + '.png')
	print('\nCreated ' + out_img_name + '.png')





