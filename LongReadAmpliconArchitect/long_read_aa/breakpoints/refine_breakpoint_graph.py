import time
start_time = time.time()
import logging
import math
import random
import pysam
import argparse
import sys
import os
import numpy as np
from collections import Counter
import cvxopt
import cvxopt.modeling
import gurobipy as gp
from gurobipy import GRB

import cigar_parsing
from long_read_aa.breakpoints import global_names
global_names.TSTART = start_time


edge_type_to_index = {'s': 1, 'c': 2, 'd': 3}


def interval_overlap(int1, int2):
	"""
	Check if two intervals in the form of [chr, s, e] overlap
	"""
	return (int1[0] == int2[0] and int(int1[1]) <= int(int2[2]) and int(int2[1]) <= int(int1[2]))


def interval_overlap_l(int1, intl):
	"""
	Check if an interval in the form of [chr, s, e] overlaps with a list of intervals
	"""
	for int2i in range(len(intl)):
		if interval_overlap(int1, intl[int2i]):
			return int2i
	return -1


def bp_match(bp1, bp2, rgap, bp_distance_cutoff):
	"""
	Check if two breakpoints match
	A breakpoint (chr1, e1, chr2, s2) must either satisfy chr1 > chr2 or chr1 == chr2 and e1 >= s2
	"""
	if bp1[0] == bp2[0] and bp1[3] == bp2[3] and bp1[2] == bp2[2] and bp1[5] == bp2[5]:
		if rgap <= 0:
			return (abs(int(bp1[1]) - int(bp2[1])) < bp_distance_cutoff[0] and \
				abs(int(bp1[4]) - int(bp2[4])) < bp_distance_cutoff[1])
		rgap_ = rgap
		consume_rgap = [0, 0]
		if bp1[2] == '+' and int(bp1[1]) <= int(bp2[1]) - bp_distance_cutoff[0]:
			rgap_ -= (int(bp2[1]) - bp_distance_cutoff[0] - int(bp1[1]) + 1)
			consume_rgap[0] = 1
		if bp1[2] == '-' and int(bp1[1]) >= int(bp2[1]) + bp_distance_cutoff[0]:
			rgap_ -= (int(bp1[1]) - int(bp2[1]) - bp_distance_cutoff[0] + 1)
			consume_rgap[0] = 1
		if bp1[5] == '+' and int(bp1[4]) <= int(bp2[4]) - bp_distance_cutoff[1]:
			rgap_ -= (int(bp2[4]) - bp_distance_cutoff[1] - int(bp1[4]) + 1)
			consume_rgap[1] = 1
		if bp1[5] == '-' and int(bp1[4]) >= int(bp2[4]) + bp_distance_cutoff[1]:
			rgap_ -= (int(bp1[4]) - int(bp2[4]) - bp_distance_cutoff[1] + 1)
			consume_rgap[1] = 1
		return (((consume_rgap[0] == 1 and rgap_ >= 0) or (abs(int(bp1[1]) - int(bp2[1])) < bp_distance_cutoff[0])) and \
			((consume_rgap[1] == 1 and rgap_ >= 0) or (abs(int(bp1[4]) - int(bp2[4])) < bp_distance_cutoff[1])))
	return False


def interval2bp(R1, R2, r = (), rgap = 0):
	"""
	Convert split/chimeric alignment to breakpoint
	"""
	if (global_names.chr_idx[R2[0]] < global_names.chr_idx[R1[0]]) or (global_names.chr_idx[R2[0]] == global_names.chr_idx[R1[0]] and R2[1] < R1[2]):
		return [R1[0], R1[2], R1[3], R2[0], R2[1], global_names.neg_plus_minus[R2[3]], r, rgap, 0]
	return [R2[0], R2[1], global_names.neg_plus_minus[R2[3]], R1[0], R1[2], R1[3], (r[0], r[2], r[1]), rgap, 1]


def bpc2bp(bp_cluster, bp_distance_cutoff):
	"""
	Call exact breakpoint from a breakpoint cluster
	"""
	#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp_cluster = %s" %(bp_cluster))
	bp = bp_cluster[0][:-2]
	bp[1] = 0 if bp[2] == '+' else 1000000000
	bp[4] = 0 if bp[5] == '+' else 1000000000
	bpr = []
	bp_stats = [0, 0, 0, 0]
	bp_stats_ = [0, 0, 0, 0, 0, 0]
	for bp_ in bp_cluster:
		bp_stats[0] += bp_[1]
		bp_stats[2] += (bp_[1] * bp_[1])
		bp_stats[1] += bp_[4]
		bp_stats[3] += (bp_[4] * bp_[4])
	for i in range(4):
		bp_stats[i] /= (len(bp_cluster) * 1.0)
	try:
		bp_stats[2] = max(bp_distance_cutoff / 2.99, math.sqrt(bp_stats[2] - bp_stats[0] * bp_stats[0]))
	except:
		bp_stats[2] = bp_distance_cutoff / 2.99
	try:
		bp_stats[3] = max(bp_distance_cutoff / 2.99, math.sqrt(bp_stats[3] - bp_stats[1] * bp_stats[1]))
	except:
		bp_stats[3] = bp_distance_cutoff / 2.99
	#print ("bp_stats", bp_stats)
	bp1_list = []
	bp4_list = []
	for bp_ in bp_cluster:
		if bp_[1] <= bp_stats[0] + 3 * bp_stats[2] and bp_[1] >= bp_stats[0] - 3 * bp_stats[2] and \
			bp_[4] <= bp_stats[1] + 3 * bp_stats[3] and bp_[4] >= bp_stats[1] - 3 * bp_stats[3]:
			bp1_list.append(bp_[1])
			bp4_list.append(bp_[4])
			#if (bp_[2] == '+' and bp_[1] > bp[1]) or (bp_[2] == '-' and bp_[1] < bp[1]):
			#	bp[1] = bp_[1]
			#if (bp_[5] == '+' and bp_[4] > bp[4]) or (bp_[5] == '-' and bp_[4] < bp[4]):
			#	bp[4] = bp_[4]
	if len(bp1_list) > 0:
		bp1_counter = Counter(bp1_list)
		if len(bp1_counter.most_common(2)) == 1 or bp1_counter.most_common(2)[0][1] > bp1_counter.most_common(2)[1][1]:
			bp[1] = bp1_counter.most_common(2)[0][0]
		else:
			if len(bp1_list) % 2 == 1:
				bp[1] = int(np.median(bp1_list))
			elif bp_[2] == '+':
				bp[1] = int(math.ceil(np.median(bp1_list)))
			else:
				bp[1] = int(math.floor(np.median(bp1_list)))
	if len(bp4_list) > 0:
		bp4_counter = Counter(bp4_list)
		if len(bp4_counter.most_common(2)) == 1 or bp4_counter.most_common(2)[0][1] > bp4_counter.most_common(2)[1][1]:
			bp[4] = bp4_counter.most_common(2)[0][0]
		else:
			if len(bp4_list) % 2 == 1:
				bp[4] = int(np.median(bp4_list))
			elif bp_[5] == '+':
				bp[4] = int(math.ceil(np.median(bp4_list)))
			else:
				bp[4] = int(math.floor(np.median(bp4_list)))
	#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp = %s" %(bp))
	bp_cluster_r = []
	for bp_ in bp_cluster:
		if bp_match(bp_, bp, bp_[7] * 1.2, [bp_distance_cutoff, bp_distance_cutoff]):
			bpr.append(bp_[6])
			bp_stats_[0] += bp_[1]
			bp_stats_[2] += (bp_[1] * bp_[1])
			bp_stats_[1] += bp_[4]
			bp_stats_[3] += (bp_[4] * bp_[4])
			if bp_[-3] == 0:
				bp_stats_[4] += bp_[-2]
				bp_stats_[5] += bp_[-1]
			else:
				bp_stats_[4] += bp_[-1]
				bp_stats_[5] += bp_[-2]
		else:
			bp_cluster_r.append(bp_)
	if len(bpr) == 0:
		return bp, bpr, [0, 0, 0, 0, 0, 0], []
	for i in range(6):
		bp_stats_[i] /= (len(bpr) * 1.0)
	#print (bp_stats_)
	try:
		bp_stats_[2] = math.sqrt(bp_stats_[2] - bp_stats_[0] * bp_stats_[0])
	except:
		bp_stats_[2] = 0
	try:
		bp_stats_[3] = math.sqrt(bp_stats_[3] - bp_stats_[1] * bp_stats_[1])
	except:
		bp_stats_[3] = 0
	return bp, bpr, bp_stats_, bp_cluster_r


def cluster_bp_list(bp_list, min_cluster_size, bp_distance_cutoff):
	"""
	Clustering the breakpoints in bp_list
	"""
	bp_dict = dict()
	for bpi in range(len(bp_list)):
		bp = bp_list[bpi]
		try:
			bp_dict[(bp[0], bp[3], bp[2], bp[5])].append(bpi)
		except:
			bp_dict[(bp[0], bp[3], bp[2], bp[5])] = [bpi]
	
	bp_clusters = []
	for bp_chr_or in bp_dict.keys():
		if len(bp_dict[bp_chr_or]) >= min_cluster_size:
			bp_clusters_ = []
			for bpi in bp_dict[bp_chr_or]:
				bp = bp_list[bpi]
				bpcim = -1
				for bpci in range(len(bp_clusters_)):
					for lbp in bp_clusters_[bpci]:
						if abs(int(bp[1]) - int(lbp[1])) < bp_distance_cutoff and \
							abs(int(bp[4]) - int(lbp[4])) < bp_distance_cutoff:
							bpcim = bpci
							break
					if bpcim > 0:
						break
				if bpcim > 0:
					bp_clusters_[bpcim].append(bp)
				else:
					bp_clusters_.append([bp])
			bp_clusters += bp_clusters_
		else:
			bp_clusters.append([bp_list[bpi] for bpi in bp_dict[bp_chr_or]])
	return bp_clusters


class bam_to_breakpoint_hybrid():

	sr_bamfh = "" # Short read bam file
	lr_bamfh = "" # Long read bam file

	min_bp_match_cutoff_ = 100 # Breakpoint matching cutoffs
	min_bp_match_cutoff = [] # Each breakpoint has a unique matching cutoff
	"""
	min_clip_cutoff = 1000 # Not in use for now
	max_interval_cutoff = 2000000 # Not in use for now
	"""
	min_cluster_cutoff = 3 # Hard cutoff for considering a long read breakpoint cluster
	max_bp_distance_cutoff = 2000 # Used for breakpoint clustering - if the distance of two breakpoint positions are greater than this cutoff, then start a new cluster 
	small_del_cutoff = 10000 # +- breakpoints (small deletions) with the two ends less than this cutoff are treated specially
	min_del_len = 10000 # The minimum length of all +- (deletion) breakpoints returned by AA  

	read_length = dict() # Map read name -> read length
	chimeric_alignments = dict() # Map read name -> chimeric alignments (two or more records for one read)
	#large_clip_alignments = dict()
	large_indel_alignments = dict() # Map read name -> alignments with one record per read but large indels showing in CIGAR string
	#chimeric_alignments_bin = dict()
	
	amplicon_intervals = [] # Amplicon intervals
	seq_edges = [] # Sequence edges
	discordant_edges = [] # Discordant edges
	concordant_edges = [] # Concordant edges
	concordant_edges_reads = []
	path_constraints = [[], []]
	valid_path_constraints = [[], [], []]
	cycles = [[], []] # cycles, paths
	cycle_weights = [[], []] # cycles, paths
	path_constraints_satisfied = [[], []] # cycles, paths
	"""
	nodes: adjacent list - keys with format (chr, pos); 
	vals = [[flags], [sequence edges], [concordant edges], [discordant edges], [source edges]]  
	"""
	nodes = dict()
	endnodes = []
	"""
	discordant_edges_pos: For each concordant edge, which side has a discordant edge
	[[], []] = No discordant edges, i.e., CN boundary breakpoint
	[[], [edge indices]] = Only incoming breakpoint edge 
	[[edge indices], []] = Only outgoing breakpoint edge
	[[edge indices], [edge indices]] = Both incoming and outgoing breakpoint edges
	"""
	discordant_edges_pos = dict() 
	small_del_indices = [] # Indices of +- breakpoints on the same chr and < small_del_cutoff in discordant_edges
	small_del_indices_ = []
	source_edges = [] # Special AA discordant edges connected to a source node 
	new_bp_list_ = [] # Long read only breakpoints (discordant edges)
	new_bp_stats_ = [] # Statistics of long read only breakpoints (discordant edges)
	
	sr_length = 0.0 # Average short read length
	min_sr_alignment_length = 30.0 # Minimum alignment length of short reads - used for CN assignment
	max_sr_insert = 2000
	normal_cov_sr = 0 # Normal short read coverage - used for CN assignment
	normal_cov_lr = 0 # Normal long read coverage - used for CN assignment
	max_CN = 0.0
	

	def __init__(self, sr_bamfile, lr_bamfile, aacfile):
		self.sr_bamfh = pysam.AlignmentFile(sr_bamfile, 'rb')
		self.lr_bamfh = pysam.AlignmentFile(lr_bamfile, 'rb')


	def read_cycle(self, aacfile):
		"""
		Read in AA cycle file for amplicon intervals
		"""
		with open(aacfile, 'r') as fp:
			for line in fp:
				s = line.strip().split()
				if s[0] == 'Interval':
					self.amplicon_intervals.append([s[2], int(s[3]), int(s[4])])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Parsed %d amplicon intervals." %(len(self.amplicon_intervals)))
		for ai in self.amplicon_intervals:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Amplicon interval %s" %(ai))


	def read_cns(self, sr_cns, lr_cns):
		"""
		Read in (cnvkit) *.cns file and estimate the normal short read and long read coverage
		"""
		for cns in [sr_cns, lr_cns]:
			cns_intervals = []
			log2_cn = []
			with open(cns, 'r') as fp:
				for line in fp:
					s = line.strip().split()
					if s[0] != "chromosome":
						cns_intervals.append([s[0], int(s[1]), int(s[2])])
						log2_cn.append(float(s[4]))
			if cns == sr_cns:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total num SR copy number segments: %d." %(len(log2_cn)))
			else:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total num LR copy number segments: %d." %(len(log2_cn)))
			log2_cn_order = np.argsort(log2_cn)
			cns_intervals_median = []
			log2_cn_median = []
			im = int(len(log2_cn_order) / 2.4)
			ip = im + 1
			total_int_len = 0
			if cns_intervals[log2_cn_order[ip]][0] != 'chrY':
				cns_intervals_median.append(cns_intervals[log2_cn_order[ip]])
				log2_cn_median.append(log2_cn[log2_cn_order[ip]])
				total_int_len += (cns_intervals[log2_cn_order[ip]][2] - cns_intervals[log2_cn_order[ip]][1] + 1)
			if cns_intervals[log2_cn_order[im]][0] != 'chrY':
				cns_intervals_median.append(cns_intervals[log2_cn_order[im]])
				log2_cn_median.append(log2_cn[log2_cn_order[im]])
				total_int_len += (cns_intervals[log2_cn_order[im]][2] - cns_intervals[log2_cn_order[im]][1] + 1)
			i = 1
			while total_int_len < 10000000:
				if cns_intervals[log2_cn_order[ip + i]][0] != 'chrY':
					cns_intervals_median.append(cns_intervals[log2_cn_order[ip + i]])
					log2_cn_median.append(log2_cn[log2_cn_order[ip + i]])
					total_int_len += (cns_intervals[log2_cn_order[ip + i]][2] - cns_intervals[log2_cn_order[ip + i]][1] + 1)
				if cns_intervals[log2_cn_order[im - i]][0] != 'chrY':
					cns_intervals_median.append(cns_intervals[log2_cn_order[im - i]])
					log2_cn_median.append(log2_cn[log2_cn_order[im - i]])
					total_int_len += (cns_intervals[log2_cn_order[im - i]][2] - cns_intervals[log2_cn_order[im - i]][1] + 1)
				i += 1
			if cns == sr_cns:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Use %d SR copy number segments." %(len(cns_intervals_median)))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total length of SR copy number segments: %d." %(total_int_len))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Average SR copy number: %f." %(np.average(log2_cn_median)))
				nnc = 0
				sr_insert_sizes = []
				for i in range(len(cns_intervals_median)):
					sr_insert_sizes += [a.template_length for a in self.sr_bamfh.fetch(cns_intervals_median[i][0], \
												cns_intervals_median[i][1], \
												cns_intervals_median[i][2]) \
							if a.is_proper_pair and not a.is_reverse and a.template_length < 10000 and a.template_length > 0]
					nnc += sum([sum(nc) for nc in self.sr_bamfh.count_coverage(cns_intervals_median[i][0], \
												cns_intervals_median[i][1], \
												cns_intervals_median[i][2] + 1)])
					#print(nnc, cns_intervals_median[i])
				self.normal_cov_sr = nnc * 1.0 / total_int_len
				sr_insert_size = np.average(sr_insert_sizes)
				sr_insert_std = np.std(sr_insert_sizes)
				self.max_sr_insert = sr_insert_size + 9 * sr_insert_std
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Estimated short read insertion size: avg = %f; std = %f" \
						%(sr_insert_size, sr_insert_std))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset maximum short read insertion size to %f." %(self.max_sr_insert))
			else:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Use %d LR copy number segments." %(len(cns_intervals_median)))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total length of LR copy number segments: %d." %(total_int_len))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Average LR copy number: %f." %(np.average(log2_cn_median)))
				nnc = 0
				for i in range(len(cns_intervals_median)):
					nnc += sum([sum(nc) for nc in self.lr_bamfh.count_coverage(cns_intervals_median[i][0], \
												cns_intervals_median[i][1], \
												cns_intervals_median[i][2] + 1, \
												quality_threshold = 0, read_callback = 'nofilter')])
				self.normal_cov_lr = nnc * 1.0 / total_int_len
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "SR normal cov = %f." %(self.normal_cov_sr))
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "LR normal cov = %f." %(self.normal_cov_lr))
		self.min_cluster_cutoff = max(self.min_cluster_cutoff, 0.5 * self.normal_cov_lr)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset min_cluster_cutoff to %f." %(self.min_cluster_cutoff))

	
	def nextminus(self, chr, pos):
		"""
		Helper function to read_graph
		Return the distance to the next position towards 3' which has incoming breakpoint edges on chr 
		"""
		cr = -1
		pos_ = pos
		while (chr, pos_, '-') in self.nodes.keys():
			if pos_ != pos and (chr, pos_ - 1, pos_) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_ - 1, pos_)][1]) > 0:
				break
			if cr >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_, '-')][1][0]][-1]
			cr = max(cr, 0) + seglen
			pos_ = pos_ + seglen
		return cr


	def lastminus(self, chr, pos):
		"""
		Helper function to read_graph
		Return the distance to the next position towards 5' which has incoming breakpoint edges on chr 
		"""
		cl = -1
		pos_ = pos
		while (chr, pos_ - 1, '+') in self.nodes.keys():
			if pos_ != pos and (chr, pos_ - 1, pos_) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_ - 1, pos_)][1]) > 0:
				break
			if cl >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_ - 1, '+')][1][0]][-1]
			cl = max(cl, 0) + seglen
			pos_ = pos_ - seglen
		return cl


	def nextplus(self, chr, pos):
		"""
		Helper function to read_graph
		Return the next position towards 3' which has outgoing breakpoint edges on chr 
		"""
		cr = -1
		pos_ = pos
		while (chr, pos_ + 1, '-') in self.nodes.keys():
			if pos_ != pos and (chr, pos_, pos_ + 1) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_, pos_ + 1)][0]) > 0:
				break
			if cr >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_ + 1, '-')][1][0]][-1]
			cr = max(cr, 0) + seglen
			pos_ = pos_ + seglen
		return cr


	def lastplus(self, chr, pos):
		"""
		Helper function to read_graph
		Return the next position towards 5' which has outgoing breakpoint edges on chr 
		"""
		cl = -1
		pos_ = pos
		while (chr, pos_, '+') in self.nodes.keys():
			if pos_ != pos and (chr, pos_, pos_ + 1) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_, pos_ + 1)][0]) > 0:
				break
			if cl >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_, '+')][1][0]][-1]
			cl = max(cl, 0) + seglen
			pos_ = pos_ - seglen
		return cl


	def read_graph(self, aagfile):
		"""
		Read in AA breakpoint graph file, including sequence, concordant, discordant and source edges.
		"""
		with open(aagfile, 'r') as fp:
			for line in fp:
				s = line.strip().split()
				if s[0] == 'sequence':
					chr = s[1][: s[1].find(':')]
					l = int(s[1][s[1].find(':') + 1: -1])
					r = int(s[2][s[2].find(':') + 1: -1])
					self.nodes[(chr, l, '-')] = [[], [len(self.seq_edges)], [], [], []]
					self.nodes[(chr, r, '+')] = [[], [len(self.seq_edges)], [], [], []]
					self.seq_edges.append([chr, l, r, int(s[-1]), -1, r - l + 1])
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + 
							"New sequence edge %s." %(self.seq_edges[-1]))
				if s[0] == 'discordant':
					t = s[1].split('->')
					t0_ = t[0].split(':')
					t1_ = t[1].split(':')
					self.discordant_edges.append([t0_[0], int(t0_[1][:-1]), t0_[1][-1], t1_[0], int(t1_[1][:-1]), t1_[1][-1],
									[float(s[2]), int(s[3])], 0, s[-2], s[-1], set([])])
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + 
							"New discordant edge %s." %(self.discordant_edges[-1]))
					if t0_[1][-1] == '+':
						if (t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1)] = [[], []]
							self.nodes[(t0_[0], int(t0_[1][:-1]), '+')][0].append(len(self.discordant_edges) - 1)
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge connected to interval ends.")
						self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1)][0].append(len(self.discordant_edges) - 1)	
					else:
						if (t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1])) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1]))] = [[], []]
							self.nodes[(t0_[0], int(t0_[1][:-1]), '-')][0].append(len(self.discordant_edges) - 1)
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge connected to interval ends.")
						self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1]))][1].append(len(self.discordant_edges) - 1)
					if t1_[1][-1] == '+':
						if (t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1)] = [[], []]
							self.nodes[(t1_[0], int(t1_[1][:-1]), '+')][0].append(len(self.discordant_edges) - 1)
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge connected to interval ends.")
						self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1)][0].append(len(self.discordant_edges) - 1)
					else:
						if (t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1])) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1]))] = [[], []]
							self.nodes[(t1_[0], int(t1_[1][:-1]), '-')][0].append(len(self.discordant_edges) - 1)
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge connected to interval ends.")
						self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1]))][1].append(len(self.discordant_edges) - 1)
					if t0_[0] == t1_[0] and t0_[1][-1] == '-' and t1_[1][-1] == '+' and \
						abs(int(t0_[1][:-1]) - int(t1_[1][:-1])) < self.small_del_cutoff:
						self.small_del_indices.append(len(self.discordant_edges) - 1)
						self.min_del_len = min(self.min_del_len, abs(int(t0_[1][:-1]) - int(t1_[1][:-1])) - 2 * self.min_bp_match_cutoff_)
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge representing small del.")
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tReset min_del_len to %d." %(self.min_del_len))
				if s[0] == 'concordant':
					t = s[1].split('->')
					t0_ = t[0].split(':')
					t1_ = t[1].split(':')
					assert t0_[1][-1] == '+' and t1_[1][-1] == '-' and t0_[0] == t1_[0] and int(t1_[1][:-1]) == int(t0_[1][:-1]) + 1
					self.concordant_edges.append([t0_[0], int(t0_[1][:-1]), t0_[1][-1], t1_[0], int(t1_[1][:-1]), t1_[1][-1],
									int(s[3]), -1, s[-2], s[-1]])
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + 
							"New concordant edge %s." %(self.concordant_edges[-1]))
					self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1)] = [[], []]
				if s[0] == 'source':
					t = s[1].split('->')
					t0_ = t[0].split(':')
					t1_ = t[1].split(':')
					assert t0_[1][:-1] == '-1'
					endnode_ = 0
					for ai in self.amplicon_intervals:
						if t1_[0] == ai[0] and (int(t1_[1][:-1]) == ai[1] or int(t1_[1][:-1]) == ai[2]):
							endnode_ = 1
							break
					if endnode_ == 0:
						self.nodes[(t1_[0], int(t1_[1][:-1]), t1_[1][-1])][4].append(len(self.source_edges))
						self.source_edges.append(['source', -1, '-', t1_[0], int(t1_[1][:-1]), t1_[1][-1], float(s[2]), float(s[3]), s[-2], s[-1]])
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + 
								"New source edge %s." %(self.source_edges[-1]))
					else:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + 
								"Skip source edge %s connecting to interval end." %(str(['source', -1, '-', t1_[0], int(t1_[1][:-1]), t1_[1][-1], float(s[2]), float(s[3]), s[-2], s[-1]])))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Parsed %d sequence edges." %(len(self.seq_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Parsed %d concordant edges." %(len(self.concordant_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Parsed %d discordant edges." %(len(self.discordant_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Parsed %d source edges." %(len(self.source_edges)))
		# Set matching cutoff for each breakpoint
		for d in self.discordant_edges:
			c1l, c1r, c2l, c2r = -1, -1, -1, -1
			if d[2] == '-':
				c1r = self.nextminus(d[0], d[1])
				c1l = self.lastminus(d[0], d[1])
			else:
				c1r = self.nextplus(d[0], d[1])
				c1l = self.lastplus(d[0], d[1])
			if d[5] == '-':
				c2r = self.nextminus(d[3], d[4])
				c2l = self.lastminus(d[3], d[4])
			else:
				c2r = self.nextplus(d[3], d[4])
				c2l = self.lastplus(d[3], d[4])
			c1, c2 = self.min_bp_match_cutoff_, self.min_bp_match_cutoff_
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Breakpoint edge %s:" %(d))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tc1l = %d; c1r = %d." %(c1l, c1r))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tc2l = %d; c2r = %d." %(c2l, c2r))
			if c1l >= 0:
				c1 = min(c1, c1l // 2 + 1)
			if c1r >= 0:
				c1 = min(c1, c1r // 2 + 1)
			if c2l >= 0:
				c2 = min(c2, c2l // 2 + 1)
			if c2r >= 0:
				c2 = min(c2, c2r // 2 + 1)
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tc1 = %d; c2 = %d." %(c1, c2))
			self.min_bp_match_cutoff.append([c1, c2])
		self.min_del_len = max(200, self.min_del_len)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tReset min_del_len to %d." %(self.min_del_len))


	def fetch(self):
		"""
		Fetch all chimeric alignments and alignments showing large indels in CIGAR string in the long read bam
		"""
		for ai in self.amplicon_intervals:
			for read in self.lr_bamfh.fetch(ai[0], ai[1], ai[2] + 1):
				rn = read.query_name
				#rf = 0
				#q_[i] >= 20
				rq = read.mapping_quality
				if rq < 20:
					continue # Fixed on 04/09/23 - Also introduced quality control on small del breakpoints
				blocks = read.get_blocks()
				for bi in range(len(blocks) - 1):
					if abs(blocks[bi + 1][0] - blocks[bi][1]) > self.min_del_len:
						#rf = 1
						try:
							self.large_indel_alignments[rn].append([ai[0], blocks[bi + 1][0], blocks[bi][1], blocks[0][0], blocks[-1][1], rq])
						except:
							self.large_indel_alignments[rn] = [[ai[0], blocks[bi + 1][0], blocks[bi][1], blocks[0][0], blocks[-1][1], rq]]
				#if rf == 1:
				#	self.large_indel_alignments[rn].append([blocks[0][0], blocks[-1][1]])
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Fetched %d reads with large indels in CIGAR." %(len(self.large_indel_alignments)))
		for read in self.lr_bamfh.fetch():
			rn = read.query_name		
			if read.flag < 256 and rn not in self.read_length:
				self.read_length[rn] = read.query_length
			try:
				sa_list = read.get_tag('SA:Z:')[:-1].split(';')
				for sa in sa_list:
					try:
						if sa not in self.chimeric_alignments[rn]:
							self.chimeric_alignments[rn].append(sa) 
					except:
						self.chimeric_alignments[rn] = [sa]
			except:
				pass
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Fetched %d chimeric reads." %(len(self.chimeric_alignments)))
		reads_wo_primary_alignment = []
		for r in self.chimeric_alignments.keys():
			if r not in self.read_length:
				logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found chimeric read without primary alignment.")
				logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tRead name: %s; Read length: N/A." %r)
				logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAll CIGAR strings: %s." %(self.chimeric_alignments[r]))
				reads_wo_primary_alignment.append(r)
				continue
			rl = self.read_length[r]
			#if r == "dd37adb1-4cc5-4eec-ae64-c0842acb3df1":
			#	print (self.chimeric_alignments[r])
			self.chimeric_alignments[r] = cigar_parsing.alignment_from_satags(self.chimeric_alignments[r], rl)
			#if interval_overlap_l(['chr3', 111555000, 111556500], self.chimeric_alignments[r][1]) > 0: #and interval_overlap_l(['chr8', 127521500, 127521700], self.chimeric_alignments[r][1]) > 0:
			#	print (self.chimeric_alignments[r])
		for r in reads_wo_primary_alignment:
			del self.chimeric_alignments[r]
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Computed alignment intervals on all chimeric reads.")


	def find_breakpoints(self):
		"""
		For each chimeric alignment, first try to match the resulting breakpoint to the list of AA breakpoints
		Cluster the unmatched breakpoints
		"""
		new_bp_list = []
		for r in self.chimeric_alignments.keys():
			r_int = self.chimeric_alignments[r][0]
			rr_int = self.chimeric_alignments[r][1]
			q_ = self.chimeric_alignments[r][2]
			bassigned = [0 for i in range(len(rr_int) - 1)]
			"""
			Breakpoint from local alignment i and i + 1
			"""
			for i in range(len(rr_int) - 1):
				overlap_i = 0
				"""
				Match the breakpoint to the list of AA breakpoints
				"""
				for di in range(len(self.discordant_edges)):
					d = self.discordant_edges[di]
					if bp_match([rr_int[i][0], rr_int[i][2], rr_int[i][3], rr_int[i + 1][0], 
							rr_int[i + 1][1], global_names.neg_plus_minus[rr_int[i + 1][3]]], d, 
							int(r_int[i + 1][0]) - int(r_int[i][1]), self.min_bp_match_cutoff[di]):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add((r, i, i + 1)) 
						overlap_i = 1
						bassigned[i] = 1
					elif bp_match([rr_int[i + 1][0], rr_int[i + 1][1], global_names.neg_plus_minus[rr_int[i + 1][3]], 
							rr_int[i][0], rr_int[i][2], rr_int[i][3]], d, 
							int(r_int[i + 1][0]) - int(r_int[i][1]), self.min_bp_match_cutoff[di]):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add((r, i + 1, i))
						overlap_i = 1
						bassigned[i] = 1
				"""
				Add unmatched breakpoint to new_bp_list
				"""
				if overlap_i == 0 and interval_overlap_l(rr_int[i], self.amplicon_intervals) >= 0 and \
					interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) >= 0:
					if rr_int[i + 1][0] != rr_int[i][0] or rr_int[i + 1][3] != rr_int[i][3]:
						#print (r, r_int, rr_int)
						if q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list.append(interval2bp(rr_int[i], rr_int[i + 1], (r, i, i + 1), int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '+':
						#print (r, r_int, rr_int)
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list.append(interval2bp(rr_int[i], rr_int[i + 1], (r, i, i + 1), int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '-':
						#print (r, r_int, rr_int)
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list.append(interval2bp(rr_int[i], rr_int[i + 1], (r, i, i + 1), int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
			"""
			Breakpoint from local alignment i - 1 and i + 1
			"""		
			for i in range(1, len(rr_int) - 1):
				overlap_i = 0
				"""
				Match the breakpoint to the list of AA breakpoints
				"""
				if bassigned[i - 1] == 0 and bassigned[i] == 0 and q_[i] < 10 and q_[i - 1] >= 20 and q_[i + 1] >= 20 and \
					interval_overlap_l(rr_int[i - 1], self.amplicon_intervals) >= 0 and interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) >= 0:
					for di in range(len(self.discordant_edges)):
						d = self.discordant_edges[di]
						if bp_match([rr_int[i - 1][0], rr_int[i - 1][2], rr_int[i - 1][3], rr_int[i + 1][0], 
								rr_int[i + 1][1], global_names.neg_plus_minus[rr_int[i + 1][3]]], d, 
								int(r_int[i + 1][0]) - int(r_int[i - 1][1]), self.min_bp_match_cutoff[di]):
							self.discordant_edges[di][7] += 1
							self.discordant_edges[di][-1].add((r, i - 1, i + 1))
							overlap_i = 1
						elif bp_match([rr_int[i + 1][0], rr_int[i + 1][1], global_names.neg_plus_minus[rr_int[i + 1][3]], 
								rr_int[i - 1][0], rr_int[i - 1][2], rr_int[i - 1][3]], d, 
								int(r_int[i + 1][0]) - int(r_int[i - 1][1]), self.min_bp_match_cutoff[di]):
							self.discordant_edges[di][7] += 1
							self.discordant_edges[di][-1].add((r, i + 1, i - 1))
							overlap_i = 1
				"""
				Add unmatched breakpoint to new_bp_list
				"""
				if overlap_i == 0 and bassigned[i - 1] == 0 and bassigned[i] == 0 and q_[i] < 10 and q_[i - 1] >= 20 and q_[i + 1] >= 20 and \
					interval_overlap_l(rr_int[i - 1], self.amplicon_intervals) >= 0 and interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) >= 0:
					if rr_int[i + 1][0] != rr_int[i - 1][0] or rr_int[i + 1][3] != rr_int[i - 1][3]:
						new_bp_list.append(interval2bp(rr_int[i - 1], rr_int[i + 1], (r, i - 1, i + 1), int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '+':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i - 1][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list.append(interval2bp(rr_int[i - 1], rr_int[i + 1], (r, i - 1, i + 1), int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '-':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i - 1][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list.append(interval2bp(rr_int[i - 1], rr_int[i + 1], (r, i - 1, i + 1), int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found %d reads with new breakpoints." %(len(new_bp_list)))
	
		new_bp_clusters = cluster_bp_list(new_bp_list, self.min_cluster_cutoff, self.max_bp_distance_cutoff)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "These reads formed %d clusters." %(len(new_bp_clusters)))
		for c in new_bp_clusters:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "New cluster of size %d." %(len(c)))
			if len(c) >= self.min_cluster_cutoff:
				num_subcluster = 0
				bp_cluster_r = c
				#print (bp_cluster_r)
				while len(bp_cluster_r) >= self.min_cluster_cutoff:
					bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(bp_cluster_r, self.min_bp_match_cutoff_)
					#print (bp, bpr, bp_cluster_r)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSubcluster %d" %(num_subcluster))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp = %s" %(bp))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNum long read support = %d" %(len(set(bpr))))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp_stats = %s" %(bp_stats_))
					if len(set(bpr)) >= max(self.normal_cov_lr, 3.0):
						if interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals) >= 0 and \
							interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals) >= 0:
							existing = 0
							for bpi in range(len(self.discordant_edges)):
								bp_ = self.discordant_edges[bpi]
								if bp[0] == bp_[0] and bp[3] == bp_[3] and bp[2] == bp_[2] and bp[5] == bp_[5] and abs(bp[1] - bp_[1]) < 200 and abs(bp[4] - bp_[4]) < 200:
									#print (bp[:6], bp_[:6], bp_[-1])
									existing = 1
									self.discordant_edges[bpi][-1] |= set(bpr)
									self.discordant_edges[bpi][7] = len(self.discordant_edges[bpi][-1]) # A read may cover a breakpoint two times, now indicating which two local alignments covers the breakpoint
									break
							if existing == 0:
								if bp[0] == bp[3] and bp[2] == '-' and bp[5] == '+' and \
									abs(bp[1] - bp[4]) < self.small_del_cutoff:
									self.small_del_indices_.append(len(self.new_bp_list_))
								self.new_bp_list_.append(bp[:-3] + [set(bpr)]) # Todo (Done 04/09/23): if one read supports a breakpoint two times we should count as two reads
								self.new_bp_stats_.append(bp_stats_)
					else:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the subcluster %d." %(num_subcluster))
					#bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(, self.min_bp_match_cutoff_)
					num_subcluster += 1
			else:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the cluster.")


	def find_smalldel_breakpoints(self):
		"""
		For each alignment record with large indels, first try to match the resulting breakpoint to the list of AA breakpoints
		Cluster the unmatched breakpoints
		"""
		new_bp_list = []
		for r in self.large_indel_alignments.keys():
			for rr_gap_i in range(len(self.large_indel_alignments[r])):
				overlap_ = 0
				rr_gap = self.large_indel_alignments[r][rr_gap_i][:3]
				rr_gap_ = rr_gap
				if rr_gap[2] > rr_gap[1]:
					rr_gap_[2] = rr_gap[1]
					rr_gap_[1] = rr_gap[2]
				for di in range(len(self.discordant_edges)):
					d = self.discordant_edges[di]
					if bp_match([rr_gap_[0], rr_gap_[1], '-', rr_gap_[0], rr_gap_[2], '+'], d, 0, \
							self.min_bp_match_cutoff[di]):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add((r, rr_gap_i, rr_gap_i))
						overlap_ = 1
				if overlap_ == 0:
					new_bp_list.append([rr_gap_[0], rr_gap_[1], '-', rr_gap_[0], rr_gap_[2], '+', (r, rr_gap_i, rr_gap_i), 0, 0, -1, -1])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found %d reads with new small del breakpoints." %(len(new_bp_list)))

		new_bp_clusters = cluster_bp_list(new_bp_list, self.min_cluster_cutoff, self.max_bp_distance_cutoff)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "These reads formed %d clusters." %(len(new_bp_clusters)))
		for c in new_bp_clusters:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "New cluster of size %d." %(len(c)))
			if len(c) >= self.min_cluster_cutoff:
				num_subcluster = 0
				bp_cluster_r = c
				while len(bp_cluster_r) >= self.min_cluster_cutoff:
					bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(bp_cluster_r, self.min_bp_match_cutoff_)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSubcluster %d" %(num_subcluster))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp = %s" %(bp))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNum long read support = %d" %(len(set(bpr))))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp_stats = %s" %(bp_stats_))
					if len(set(bpr)) >= max(self.normal_cov_lr, 3.0):
						if interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals) >= 0 and \
							interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals) >= 0:
							existing = 0
							for bpi in range(len(self.discordant_edges)):
								bp_ = self.discordant_edges[bpi]
								if bp[0] == bp_[0] and bp[3] == bp_[3] and bp[2] == bp_[2] and bp[5] == bp_[5] and abs(bp[1] - bp_[1]) < 200 and abs(bp[4] - bp_[4]) < 200:
									#print (bp[:6], bp_[:6], bp_[-1])
									self.discordant_edges[bpi][-1] |= set(bpr)
									self.discordant_edges[bpi][7] = len(self.discordant_edges[bpi][-1])
									existing = 1
									break
							if existing == 0:
								self.small_del_indices_.append(len(self.new_bp_list_))
								self.new_bp_list_.append(bp[:-3] + [set(bpr)])  # Todo: if one read supports a breakpoint two times we should count as two reads
								self.new_bp_stats_.append(bp_stats_)
					else:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the subcluster %d." %(num_subcluster))
					num_subcluster += 1
			else:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the cluster.")


	def split_seg_bp(self):
		"""
		Split AA sequence edges where a long read breakpoint occurs in the middle
		"""
		split_seg = dict()
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Processing %d new LR only breakpoints." %(len(self.new_bp_list_)))
		for bpi in range(len(self.new_bp_list_)):
			bp = self.new_bp_list_[bpi]
			nsplit = [1, 1]
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "New breakpoint %s: nsplit = %s." %(bp[:6] + [len(bp[6])], nsplit))
			for seg in self.seq_edges:
				ep, em = -1, -1
				try:
					ep = 1 if len(self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0]) > 0 else 0
					em = 1 if len(self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1]) > 0 else 0
				except:
					pass
				if bp[0] == seg[0]:
					if bp[2] == '+':
						if seg[2] - self.min_bp_match_cutoff_ < bp[1] <= seg[2]:
							self.new_bp_list_[bpi][1] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
						elif seg[2] < bp[1] < seg[2] + self.min_bp_match_cutoff_ and em == 0:
							self.new_bp_list_[bpi][1] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
					if bp[2] == '-':
						if seg[2] < bp[1] <= seg[2] + self.min_bp_match_cutoff_ and ep + em >= 0:
							self.new_bp_list_[bpi][1] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
						elif seg[2] - self.min_bp_match_cutoff_ < bp[1] < seg[2] and ep == 0:
							self.new_bp_list_[bpi][1] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
				if bp[3] == seg[0]:
					if bp[5] == '+':
						if seg[2] - self.min_bp_match_cutoff_ < bp[4] <= seg[2]:
							self.new_bp_list_[bpi][4] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
						elif seg[2] < bp[4] < seg[2] + self.min_bp_match_cutoff_ and em == 0:
							self.new_bp_list_[bpi][4] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
					if bp[5] == '-':
						if seg[2] < bp[4] <= seg[2] + self.min_bp_match_cutoff_ and ep + em >= 0:
							self.new_bp_list_[bpi][4] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
						elif seg[2] - self.min_bp_match_cutoff_ < bp[4] < seg[2] and ep == 0:
							self.new_bp_list_[bpi][4] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "New breakpoint %s: nsplit = %s." %(bp[:6] + [len(bp[6])], nsplit))
			for segi in range(len(self.seq_edges)):	
				seg = self.seq_edges[segi]
				if nsplit[0] == 1 and bp[0] == seg[0] and seg[1] < int(bp[1]) < seg[2]:
					if bp[2] == '+':
						try:
							split_seg[segi].append((int(bp[1]), int(bp[1]) + 1, bpi, 1, '+'))
						except:
							split_seg[segi] = [(int(bp[1]), int(bp[1]) + 1, bpi, 1, '+')]
						try:
							self.discordant_edges_pos[(seg[0], int(bp[1]), int(bp[1]) + 1)][0].append(len(self.discordant_edges) + bpi)
						except:
							self.discordant_edges_pos[(seg[0], int(bp[1]), int(bp[1]) + 1)] = [[len(self.discordant_edges) + bpi], []]
					if bp[2] == '-':
						try:
							split_seg[segi].append((int(bp[1]) - 1, int(bp[1]), bpi, 1, '-'))
						except:
							split_seg[segi] = [(int(bp[1]) - 1, int(bp[1]), bpi, 1, '-')]
						try:
							self.discordant_edges_pos[(seg[0], int(bp[1]) - 1, int(bp[1]))][1].append(len(self.discordant_edges) + bpi)
						except:
							self.discordant_edges_pos[(seg[0], int(bp[1]) - 1, int(bp[1]))] = [[], [len(self.discordant_edges) + bpi]]
				if nsplit[1] == 1 and bp[3] == seg[0] and seg[1] < int(bp[4]) < seg[2]:
					if bp[5] == '+':
						try:
							split_seg[segi].append((int(bp[4]), int(bp[4]) + 1, bpi, 4, '+'))
						except:
							split_seg[segi] = [(int(bp[4]), int(bp[4]) + 1, bpi, 4, '+')]
						try:
							self.discordant_edges_pos[(seg[0], int(bp[4]), int(bp[4]) + 1)][0].append(len(self.discordant_edges) + bpi)
						except:
							self.discordant_edges_pos[(seg[0], int(bp[4]), int(bp[4]) + 1)] = [[len(self.discordant_edges) + bpi], []]
					if bp[5] == '-':
						try:
							split_seg[segi].append((int(bp[4]) - 1, int(bp[4]), bpi, 4, '-'))
						except:
							split_seg[segi] = [(int(bp[4]) - 1, int(bp[4]), bpi, 4, '-')]
						try:
							self.discordant_edges_pos[(seg[0], int(bp[4]) - 1, int(bp[4]))][1].append(len(self.discordant_edges) + bpi)
						except:
							self.discordant_edges_pos[(seg[0], int(bp[4]) - 1, int(bp[4]))] = [[], [len(self.discordant_edges) + bpi]]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will split the following %d sequence edges." %(len(split_seg)))
		for segi in split_seg.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will split the sequence edge at index %d." %(segi))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSplit at %s." %(split_seg[segi]))
		new_segs = dict()
		for segi in split_seg.keys():
			split_seg[segi].sort(key = lambda item: item[0])
			#print (split_seg[segi])
			sseg = self.seq_edges[segi]	
			new_segs[segi] = []
			lssi = 0
			del_list = []
			for ssi in range(len(split_seg[segi])):
				spp = split_seg[segi][ssi]
				if ssi > 0 and int(spp[0]) - int(split_seg[segi][lssi][1]) <= self.min_bp_match_cutoff_ and \
					spp[4] == split_seg[segi][lssi][4]:
					if self.new_bp_list_[spp[2]][spp[3] + 1] == '+':
						self.new_bp_list_[spp[2]][spp[3]] = split_seg[segi][lssi][0]
					else:
						self.new_bp_list_[spp[2]][spp[3]] = split_seg[segi][lssi][1]
					del_list.append(ssi)
				else:
					lssi = ssi
			for dssi in del_list[::-1]:
				del split_seg[segi][dssi]
			#print (split_seg[segi])
			for ssi in range(len(split_seg[segi])):
				if ssi == 0:
					new_segs[segi].append([sseg[0], sseg[1], split_seg[segi][ssi][0], -1, -1, 
								int(split_seg[segi][ssi][0]) - int(sseg[1]) + 1])
				else:
					new_segs[segi].append([sseg[0], split_seg[segi][ssi - 1][1], split_seg[segi][ssi][0], 
								-1, -1, int(split_seg[segi][ssi][0]) - \
								int(split_seg[segi][ssi - 1][1]) + 1])
			new_segs[segi].append([sseg[0], split_seg[segi][-1][1], sseg[2], -1, -1, 
					int(sseg[2]) - int(split_seg[segi][-1][1]) + 1])
		for segi in sorted(split_seg.keys())[::-1]:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Processing sequence edge at index %d." %(segi))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDeleted sequence edge %s." %(self.seq_edges[segi]))
			del self.seq_edges[segi]
			for nsi in range(len(new_segs[segi])):
				self.seq_edges.insert(segi + nsi, new_segs[segi][nsi]) # Add sequence edges
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNew sequence edge %s at index %d." %(new_segs[segi][nsi], segi + nsi))
				if nsi < len(new_segs[segi]) - 1:
					self.concordant_edges.append([new_segs[segi][nsi][0], new_segs[segi][nsi][2], '+', 
						new_segs[segi][nsi + 1][0], new_segs[segi][nsi + 1][1], '-', -1, -1, 'None', 'None'])
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNew concordant edge %s." %(self.concordant_edges[-1]))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tdiscordant_edges_pos = %s." 
							%(self.discordant_edges_pos[(new_segs[segi][nsi][0], new_segs[segi][nsi][2], new_segs[segi][nsi + 1][1])]))
					self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi][2], '+')] = [[], [], [], [], []]
					self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi + 1][1], '-')] = [[], [], [], [], []]
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNew node %s; adjacent list = %s." 
							%(str((new_segs[segi][nsi][0], new_segs[segi][nsi][2], '+')), self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi][2], '+')]))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNew node %s; adjacent list = %s." 
							%(str((new_segs[segi][nsi][0], new_segs[segi][nsi + 1][1], '-')), self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi + 1][1], '-')]))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num sequence edges in breakpoint graph = %d." %(len(self.seq_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num concordant edges in breakpoint graph = %d." %(len(self.concordant_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num discordant edges in breakpoint graph = %d." %(len(self.discordant_edges) + len(self.new_bp_list_)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num source edges in breakpoint graph = %d." %(len(self.source_edges)))
				

	def del_sr_bp(self, aa_downsampled = False):
		"""
		Delete short read breakpoints without long read coverage;
		Merge the corresponding sequence edges where no breakpoints left at the end after breakpoint edge deletion.
		"""
		del_list = [] # The list of short read ONLY breakpoints to be deleted
		source_del_list = [] # The list of source edges to be deleted
		bpi_map = dict()
		bpi_ = 0
		if aa_downsampled:
			for bpi in range(len(self.discordant_edges)):
				if self.discordant_edges[bpi][7] < 0.5 * self.normal_cov_lr and (self.discordant_edges[bpi][6][0] < 1.0 or self.discordant_edges[bpi][6][1] < min(10.0, 0.5 * self.normal_cov_sr)):
					del_list.append(bpi)
				else:
					bpi_map[bpi] = bpi_
					bpi_ += 1
		else:
			for bpi in range(len(self.discordant_edges)):
				if self.discordant_edges[bpi][7] < 0.5 * self.normal_cov_lr and (self.discordant_edges[bpi][6][0] < 1.0 or self.discordant_edges[bpi][6][1] < 0.5 * self.normal_cov_sr):
					del_list.append(bpi)
				else:
					bpi_map[bpi] = bpi_
					bpi_ += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will delete the discordant edges at index %s." %(del_list))
		srci_map = dict()
		srci_ = 0
		if aa_downsampled:
			for srci in range(len(self.source_edges)):
				if self.source_edges[srci][7] > 0.5 and (self.source_edges[srci][6] < 1.0 or self.source_edges[srci][7] < min(10.0, 0.5 * self.normal_cov_sr)):
					source_del_list.append(srci)
				else:
					srci_map[srci] = srci_
					srci_ += 1
		else:
			for srci in range(len(self.source_edges)):
				if self.source_edges[srci][7] > 0.5 and (self.source_edges[srci][6] < 1.0 or self.source_edges[srci][7] < 0.5 * self.normal_cov_sr):
					source_del_list.append(srci)
				else:
					srci_map[srci] = srci_
					srci_ += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will delete the source edge at index %s." %(source_del_list))
		"""
		Delete breakpoint and source edges
		"""
		for bpi in del_list[::-1]:
			del self.discordant_edges[bpi]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "After deletion, there are %d remaining discordant edges." %(len(self.discordant_edges)))
		for srci in source_del_list[::-1]:
			del self.source_edges[srci]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "After deletion, there are %d remaining source edges." %(len(self.source_edges)))
		"""
		Modify the indices in node[0] and node[4] for each node in the adjacent list
		"""
		for node in self.nodes.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Modifying the adjacent list for node %s." %(str(node)))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBefore modification: %s." %(self.nodes[node]))
			if len(self.nodes[node][0]) > 0:
				assert len(self.nodes[node][0]) == 1
				if self.nodes[node][0][0] in del_list:
					del self.nodes[node][0][0]
				else:
					self.nodes[node][0][0] = bpi_map[self.nodes[node][0][0]]
			if len(self.nodes[node][4]) > 0:
				#assert len(self.nodes[node][4]) == 1
				if self.nodes[node][4][0] in source_del_list:
					del self.nodes[node][4][0]
				else:
					self.nodes[node][4][0] = srci_map[self.nodes[node][4][0]]
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAfter modification: %s." %(self.nodes[node]))
		"""
		Modify the indices in small_del_indices
		"""
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Modifying the indices in small_del_indices.")
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBefore modification: %s." %(self.small_del_indices))
		for sdi in range(len(self.small_del_indices))[::-1]:
			if self.small_del_indices[sdi] in del_list:
				del self.small_del_indices[sdi]
		for sdi in range(len(self.small_del_indices)):
			self.small_del_indices[sdi] = bpi_map[self.small_del_indices[sdi]]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAfter modification: %s." %(self.small_del_indices))
		"""
		Modify the indices in discordant_edges_pos
		This step has to be done before calling split_seg_bp()
		"""
		for concordant_edge in self.discordant_edges_pos.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Modifying the discordant_edges_pos info for concordant edge %s." %(str(concordant_edge)))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBefore modification: %s." %(self.discordant_edges_pos[concordant_edge]))
			for j in [0, 1]:
				for i in range(len(self.discordant_edges_pos[concordant_edge][j]))[::-1]:
					if self.discordant_edges_pos[concordant_edge][j][i] in del_list:
						del self.discordant_edges_pos[concordant_edge][j][i]
				for i in range(len(self.discordant_edges_pos[concordant_edge][j])):
					self.discordant_edges_pos[concordant_edge][j][i] = bpi_map[self.discordant_edges_pos[concordant_edge][j][i]]
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAfter modification: %s." %(self.discordant_edges_pos[concordant_edge]))
		"""
		Merge sequence edges resulted from breakpoint edge deletion
		Note that this step should be done after breakpoint matching and will not lead to the reset of min_bp_match_cutoff
		"""
		c_del_list = [] # The list of concordant edges to be deleted
		seg_del_list = [] # The list of sequence edges to be deleted
		for ci in range(len(self.concordant_edges)):
			ce = self.concordant_edges[ci]
			node1 = (ce[0], ce[1], ce[2])
			node2 = (ce[3], ce[4], ce[5])
			#print(node1, node2)
			if len(self.nodes[node1][4]) == 0 and len(self.nodes[node2][4]) == 0 and \
				len(self.discordant_edges_pos[(ce[0], ce[1], ce[4])][0]) == 0 and \
				len(self.discordant_edges_pos[(ce[0], ce[1], ce[4])][1]) == 0:
				segi1 = self.nodes[node1][1][0]
				segi2 = self.nodes[node2][1][0]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Merging sequence edges at index %d and %d." %(segi1, segi2))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBefore merging: %s; %s." %(self.seq_edges[segi1], self.seq_edges[segi2]))
				self.seq_edges[segi2][1] = self.seq_edges[segi1][1]
				self.seq_edges[segi2][3] = -1
				self.seq_edges[segi2][-1] = self.seq_edges[segi2][2] - self.seq_edges[segi2][1] + 1
				seg_del_list.append(segi1)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDeleting nodes %s and %s." %(str(node1), str(node2)))
				del self.nodes[node1]
				del self.nodes[node2]
				c_del_list.append(ci)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAfter merging: %s." %(self.seq_edges[segi2]))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will delete the sequence edges at index %s." %(seg_del_list))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will delete the concordant edges at index %s." %(c_del_list))
		for segi in seg_del_list[::-1]:
			del self.seq_edges[segi]
		for ci in c_del_list[::-1]:
			del self.concordant_edges[ci]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "After deletion, there are %d remaining sequence edges." %(len(self.seq_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "After deletion, there are %d remaining concordant edges." %(len(self.concordant_edges)))


	def assign_cov(self, lr_seq_dist = 'poisson'):
		"""
		Extract the short read and long read coverage from bam file, if missing, for each sequence edge 
		"""
		max_seg = -1
		avg_rl = 0.0
		for seg in self.seq_edges:
			if seg[3] == -1:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Finding SR cov for sequence edge %s." %(seg))
				rl_list = [a.infer_read_length() for a in self.sr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if a.infer_read_length()]
				len_rl_list = len(rl_list)
				seg[3] = [len_rl_list, 'full']
				if len_rl_list > max_seg:
					max_seg = len_rl_list
					avg_rl = np.average(rl_list) # Average SR length
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "SR cov assigned for sequence edge %s." %(seg))
			if seg[4] == -1:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Finding LR cov for sequence edge %s." %(seg))
				"""
				For long read, use the total number of nucleotides
				"""
				rl_list = [a.infer_read_length() for a in self.lr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if a.infer_read_length()]
				if lr_seq_dist == 'poisson':
					seg[4] = [len(rl_list), np.average(rl_list)]
					if seg[4][0] == 0:
						seg[4] = [0, 10000] # 
				else:
					seg[4] = [len(rl_list), sum([sum(nc) for nc in self.lr_bamfh.count_coverage(seg[0], seg[1], seg[2] + 1, quality_threshold = 0, read_callback = 'nofilter')])]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "LR cov assigned for sequence edge %s." %(seg))
		if self.sr_length == 0.0:
			if avg_rl == 0.0:
				logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "No new sequence edges found, compute SR read length on the first sequence edge in seq_edges.")
				seg = self.seq_edges[0]
				rl_list = [a.infer_read_length() for a in self.sr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if a.infer_read_length()]
				avg_rl = np.average(rl_list) 
			self.sr_length = avg_rl
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset sr_length to %f." %(self.sr_length))
		"""
		Extract the short read and long read coverage from bam file, if missing, for each concordant edge 
		"""
		for eci in range(len(self.concordant_edges)):
			ec = self.concordant_edges[eci]
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Finding cov for concordant edge %s." %(ec))
			if ec[6] == -1:
				crs_ = set([a.query_name for a in self.sr_bamfh.fetch(ec[0], \
									max(1, ec[1] - self.max_sr_insert), \
									ec[1]) \
						if not a.is_unmapped and a.is_proper_pair and not a.is_reverse and \
						a.next_reference_name == ec[0] and a.next_reference_start >= ec[1] and \
						a.reference_start < ec[1] and a.next_reference_start < a.reference_start + self.max_sr_insert - self.sr_length])
				self.concordant_edges[eci][6] = [len(crs_), 'full']
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSR cov assigned for concordant edge %s." %(self.concordant_edges[eci]))
			rls = set([a.query_name for a in self.lr_bamfh.fetch(contig = ec[0], start = ec[1], stop = ec[1] + 1)])
			rrs = set([a.query_name for a in self.lr_bamfh.fetch(contig = ec[3], start = ec[4], stop = ec[4] + 1)])
			rls1 = set([a.query_name for a in self.lr_bamfh.fetch(contig = ec[0], start = ec[1] - self.min_bp_match_cutoff_ - 1, stop = ec[1] - self.min_bp_match_cutoff_)])
			rrs1 = set([a.query_name for a in self.lr_bamfh.fetch(contig = ec[3], start = ec[4] + self.min_bp_match_cutoff_, stop = ec[4] + self.min_bp_match_cutoff_ + 1)]) 
			rbps = set([])
			for bpi in self.discordant_edges_pos[(ec[0], ec[1], ec[4])][0]:
				if bpi >= len(self.discordant_edges):
					rbps |= set([r_[0] for r_ in self.new_bp_list_[bpi - len(self.discordant_edges)][-1]])
				else:
					rbps |= set([r_[0] for r_ in self.discordant_edges[bpi][-1]])
			for bpi in self.discordant_edges_pos[(ec[0], ec[1], ec[4])][1]:
				if bpi >= len(self.discordant_edges):
					rbps |= set([r_[0] for r_ in self.new_bp_list_[bpi - len(self.discordant_edges)][-1]])
				else:
					rbps |= set([r_[0] for r_ in self.discordant_edges[bpi][-1]])
			self.concordant_edges_reads.append(rls | rrs)
			self.concordant_edges[eci][7] = len((rls & rrs & rls1 & rrs1) - rbps)
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tLR cov assigned for concordant edge %s." %(self.concordant_edges[eci]))


	def output_breakpoint_graph(self, ogfile, aa_downsampled = False):
		"""
		Write a breakpoint graph to file in AA graph format
		"""
		with open(ogfile, 'w') as fp:
			fp.write("SequenceEdge: StartPosition, EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads, Size\n")
			for segi in range(len(self.seq_edges)):
				sseg = self.seq_edges[segi]
				if type(sseg[3]) is list:
					if aa_downsampled and self.normal_cov_sr > 10:
						fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%s\t%d\t%s\n" %(sseg[0], sseg[1], sseg[0], sseg[2], sseg[-1], \
							int(math.round(sseg[3][0] * 10.0 / self.normal_cov_sr)), sseg[4][0], str(sseg[5])))
					else:
						fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%s\t%d\t%s\n" %(sseg[0], sseg[1], sseg[0], sseg[2], sseg[-1], sseg[3][0], sseg[4][0], str(sseg[5])))
				else:
					fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%s\t%d\t%s\n" %(sseg[0], sseg[1], sseg[0], sseg[2], sseg[-1], sseg[3], sseg[4][0], str(sseg[5])))
			fp.write("BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads, ")
			fp.write("HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence\n")
			for es in self.source_edges:
				fp.write("source\t%s:%s%s->%s:%s%s\t%f\t%s\t-1\t%s\t%s\n" %(es[0], es[1], es[2], es[3], 
					es[4], es[5], es[-1], str(es[7]), es[8], es[9]))
			for ec in self.concordant_edges:
				if type(ec[6]) is list:
					if aa_downsampled and self.normal_cov_sr > 10:
						fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\t%s\t%s\n" %(ec[0], ec[1], ec[2], ec[3], 
							ec[4], ec[5], ec[-1], int(math.round(ec[6] * 10.0 / self.normal_cov_sr)), ec[7], ec[8], ec[9]))
					else:
						fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\t%s\t%s\n" %(ec[0], ec[1], ec[2], ec[3], 
							ec[4], ec[5], ec[-1], ec[6][0], ec[7], ec[8], ec[9]))
				else:
					fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\t%s\t%s\n" %(ec[0], ec[1], ec[2], ec[3], 
						ec[4], ec[5], ec[-1], ec[6], ec[7], ec[8], ec[9]))
			for ed in self.discordant_edges:
				fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\t%s\t%s\n" %(ed[0], ed[1], ed[2], ed[3], 
					ed[4], ed[5], ed[-1], ed[6][1], ed[7], ed[8], ed[9]))
			for bp in self.new_bp_list_:
				fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t-1\t%d\tNone\tNone\n" %(bp[0], bp[1], bp[2],
					bp[3], bp[4], bp[5], bp[-1], len(bp[6])))
	

	def output_breakpoint_info(self, obpfile):
		"""
		Write the list of breakpoints to file
		"""
		with open(obpfile, 'w') as fp:
			fp.write("chr1\tpos1\tchr2\tpos2\torientation\tsr_support\tlr_support\tlr_info=[avg1, avg2, std1, std2, mapq1, mapq2]\n")
			for ed in self.discordant_edges:
				fp.write("%s\t%s\t%s\t%s\t%s%s\t%f\t%d\t%d\tN/A\n" %(ed[3], ed[4], ed[0], ed[1], ed[5], ed[2], 
					ed[6][0], ed[6][1], ed[7]))
			for bpi in range(len(self.new_bp_list_)):
				bp = self.new_bp_list_[bpi]
				bp_stats = self.new_bp_stats_[bpi]
				fp.write("%s\t%s\t%s\t%s\t%s%s\t-1\t%d\t%s\n" 
					%(bp[3], bp[4], bp[0], bp[1], bp[5], bp[2], len(bp[-1]), bp_stats))
	

	def closebam(self):
		"""
		Close the short read and long read bam file
		"""
		self.sr_bamfh.close()
		self.lr_bamfh.close()


	def assign_cn(self, lr_seq_dist = 'poisson', aa_downsampled = False):
		"""
		Compute the maximum likelihood CN assignment on each edge
		"""
		lseg = len(self.seq_edges)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lnbp = len(self.new_bp_list_)
		lsrc = len(self.source_edges)

		"""
		Prepare adjacency list
		"""
		#endnodes = [] # Nodes corresponding to interval ends
		for segi in range(lseg):
			sseg = self.seq_edges[segi]
			if len(self.nodes[(sseg[0], sseg[1], '-')][1]) == 0:
				self.nodes[(sseg[0], sseg[1], '-')][1].append(segi)
			else:
				self.nodes[(sseg[0], sseg[1], '-')][1] = [segi]
			if len(self.nodes[(sseg[0], sseg[2], '+')][1]) == 0:
				self.nodes[(sseg[0], sseg[2], '+')][1].append(segi)
			else:
				self.nodes[(sseg[0], sseg[2], '+')][1] = [segi]
			if segi == 0:
				self.endnodes.append((sseg[0], sseg[1], '-'))
			else:
				lastseg = self.seq_edges[segi - 1]
				if lastseg[0] != sseg[0] or lastseg[2] + 1 != sseg[1]:
					self.endnodes.append((lastseg[0], lastseg[2], '+'))
					self.endnodes.append((sseg[0], sseg[1], '-'))
				if segi == lseg - 1:
					self.endnodes.append((sseg[0], sseg[2], '+'))
		if lseg == 1:
			sseg = self.seq_edges[0]
			self.endnodes.append((sseg[0], sseg[2], '+'))
		assert len(self.endnodes) == len(self.amplicon_intervals) * 2
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "The following nodes correspond to interval ends.")
		for node in self.endnodes:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNode %s." %(str(node)))
		for eci in range(lc):
			ec = self.concordant_edges[eci]
			self.nodes[(ec[0], ec[1], ec[2])][2].append(eci)
			self.nodes[(ec[3], ec[4], ec[5])][2].append(eci)
		for edi in range(ld):
			ed = self.discordant_edges[edi]
			self.nodes[(ed[0], ed[1], ed[2])][3].append(edi)
			self.nodes[(ed[3], ed[4], ed[5])][3].append(edi)
		for bpi in range(lnbp):
			bp = self.new_bp_list_[bpi]
			self.nodes[(bp[0], bp[1], bp[2])][3].append(ld + bpi)
			self.nodes[(bp[3], bp[4], bp[5])][3].append(ld + bpi)

		# Delete extra sequence edges; may move to the deletion functions
		del_list = []
		if ld + lnbp > 0:
			for segi in range(lseg):
				sseg = self.seq_edges[segi]
				node1 = (sseg[0], sseg[1], '-')
				node2 = (sseg[0], sseg[2], '+')
				s1 = len(self.nodes[node1][2]) + len(self.nodes[node1][3]) + len(self.nodes[node1][4])
				s2 = len(self.nodes[node2][2]) + len(self.nodes[node2][3]) + len(self.nodes[node2][4])
				
				if s1 + s2 == 0:
					del_list.append(segi)
			#print (del_list)
			#print (self.seq_edges)
			for segi in del_list[::-1]:
				#sseg = self.seq_edges[segi]
				ai = self.seq_edges[segi][:3]
				if ai in self.amplicon_intervals:
					del self.amplicon_intervals[self.amplicon_intervals.index(ai)]
				node1 = (self.seq_edges[segi][0], self.seq_edges[segi][1], '-')
				node2 = (self.seq_edges[segi][0], self.seq_edges[segi][2], '+')
				del self.seq_edges[segi]
				del self.nodes[node1]
				del self.nodes[node2]
				del self.endnodes[self.endnodes.index(node1)]
				del self.endnodes[self.endnodes.index(node2)]
			for segi in range(len(self.seq_edges)):
				node1 = (self.seq_edges[segi][0], self.seq_edges[segi][1], '-')
				node2 = (self.seq_edges[segi][0], self.seq_edges[segi][2], '+')
				self.nodes[node1][1][0] = segi
				self.nodes[node2][1][0] = segi
			#print (self.seq_edges)
		lseg = len(self.seq_edges)

		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Updated adjacent list for %d nodes:" %(len(self.nodes)))
		for node in self.nodes.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Node %s; adjacent list = %s." %(str(node), self.nodes[node]))

		nconstraints = len([node for node in self.nodes.keys() if node not in self.endnodes or len(self.nodes[node][0]) > 0])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d constraints for cvxopt." %(nconstraints))
		nvariables = lseg + lc + ld + lnbp + lsrc
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d variables for cvxopt." %(nvariables))

		wcn = []
		wlncn = []
		wlrseg = []
		if lr_seq_dist == 'poisson':
			wcn = [(self.normal_cov_sr * sseg[-1] / self.sr_length) + \
				(self.normal_cov_lr * sseg[-1] / sseg[4][1]) for sseg in self.seq_edges]
		else:
			wcn = [(self.normal_cov_sr * sseg[-1] / self.sr_length + 0.5 * self.normal_cov_lr * sseg[-1]) for sseg in self.seq_edges]
		wcn += [self.normal_cov_sr * (self.sr_length - 1) / self.sr_length + self.normal_cov_lr \
				for eci in range(lc)]
		wcn += [self.normal_cov_sr * (self.sr_length - 2 * self.min_sr_alignment_length) / \
				self.sr_length + self.normal_cov_lr for edi in range(ld)]
		wcn += [self.normal_cov_lr for bpi in range(lnbp)]
		wcn += [self.normal_cov_sr * (self.sr_length - 2 * self.min_sr_alignment_length) / self.sr_length \
				for es in self.source_edges]
		if lr_seq_dist == 'poisson':
			for sseg in self.seq_edges:
				if type(sseg[3]) is list:
					wlncn.append((sseg[3][0] + sseg[4][0]) * 1.0)
				else:
					if aa_downsampled and self.normal_cov_sr > 10:
						wlncn.append(sseg[3] * self.normal_cov_sr / 10.0 + sseg[4][0])
					else:
						wlncn.append((sseg[3] + sseg[4][0]) * 1.0)
		else:
			for sseg in self.seq_edges:
				if type(sseg[3]) is list:
					wlncn.append(sseg[3][0] - 0.5)
				else:
					if aa_downsampled and self.normal_cov_sr > 10:
						wlncn.append(sseg[3] * self.normal_cov_sr / 10.0 - 0.5)
					else:
						wlncn.append(sseg[3] - 0.5)
		for ec in self.concordant_edges:
			if type(ec[6]) is list:
				wlncn.append((ec[6][0] + ec[7]) * 1.0)
			else:
				if aa_downsampled and self.normal_cov_sr > 10:
					wlncn.append(ec[6] * self.normal_cov_sr / 10.0 + ec[7])
				else:
					wlncn.append((ec[6] + ec[7]) * 1.0)
		for ed in self.discordant_edges:
			if aa_downsampled and self.normal_cov_sr > 10:
				wlncn.append(ed[6][1] * self.normal_cov_sr / 10.0 + ed[7])
			else:
				wlncn.append((ed[6][1] + ed[7]) * 1.0)
		wlncn += [len(bp[6]) * 1.0 for bp in self.new_bp_list_]
		if aa_downsampled and self.normal_cov_sr > 10:
			wlncn += [es[7] * self.normal_cov_sr / 10.0 if es[7] >= 1 else 0.1 for es in self.source_edges]
		else:
			wlncn += [es[7] * 1.0 if es[7] >= 1 else 0.1 for es in self.source_edges]
		if lr_seq_dist != 'poisson':
			wlrseg = [(0.5 * sseg[4][1] * sseg[4][1] / (self.normal_cov_lr * sseg[-1])) for sseg in self.seq_edges]
			wlrseg += [0.0 for ec in self.concordant_edges]
			wlrseg += [0.0 for ed in self.discordant_edges]
			wlrseg += [0.0 for bp in self.new_bp_list_]
			wlrseg += [0.0 for es in self.source_edges]
		wcn = cvxopt.matrix(wcn)
		wlncn = cvxopt.matrix(wlncn)
		wlrseg = cvxopt.matrix(wlrseg)
		
		#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "The gradient of CN vector is: wcn * CN - wlncn * CN:")
		#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\twcn = %s;" %(wcn))
		#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\twlncn = %s." %(wlncn))
		
		ci = 0
		balance_constraints = np.zeros([nconstraints, nvariables])
		for node in self.nodes.keys():
			if node not in self.endnodes or len(self.nodes[node][0]) > 0:
				for segi in self.nodes[node][1]:
					balance_constraints[ci][segi] = 1
				for eci in self.nodes[node][2]:
					balance_constraints[ci][lseg + eci] = -1
				for edi in self.nodes[node][3]:
					balance_constraints[ci][lseg + lc + edi] = -1
				for esi in self.nodes[node][4]:
					balance_constraints[ci][lseg + lc + ld + lnbp + esi] = -1
				ci += 1
		balance_constraints = cvxopt.matrix(balance_constraints)

		# Convex optimization function required by cvxopt
		def F_normal(x = None, z = None):
			if x is None: 
				return 0, cvxopt.matrix(1.0, (nvariables, 1))
			if min(x) <= 0.0: 
				return None
			f = cvxopt.modeling.dot(wlrseg, x**-1) + cvxopt.modeling.dot(wcn, x) - cvxopt.modeling.dot(wlncn, cvxopt.log(x))
			Df = (wcn - cvxopt.mul(wlncn, x**-1) - cvxopt.mul(wlrseg, x**-2)).T
			if z is None: 
				return f, Df
			H = cvxopt.spdiag(z[0] * (cvxopt.mul(wlncn, x**-2) + cvxopt.mul(2.0 * wlrseg, x**-3)))
			return f, Df, H
		def F_poisson(x = None, z = None):
			if x is None: 
				return 0, cvxopt.matrix(1.0, (nvariables, 1))
			if min(x) <= 0.0: 
				return None
			f = cvxopt.modeling.dot(wcn, x) - cvxopt.modeling.dot(wlncn, cvxopt.log(x))
			Df = (wcn - cvxopt.mul(wlncn, x**-1)).T
			if z is None: 
				return f, Df
			H = cvxopt.spdiag(z[0] * cvxopt.mul(wlncn, x**-2))
			return f, Df, H
		
		options = {'maxiters': 1000, 'show_progress': False}
		sol = ''
		if nconstraints > 0:
			if lr_seq_dist == 'poisson':
				sol = cvxopt.solvers.cp(F_poisson, A = balance_constraints, 
							b = cvxopt.matrix([0.0 for i in range(nconstraints)]), 
							kktsolver = 'ldl', options = options)
			else:
				sol = cvxopt.solvers.cp(F_normal, A = balance_constraints, 
							b = cvxopt.matrix([0.0 for i in range(nconstraints)]), 
							kktsolver = 'ldl', options = options)
			if sol['status'] == "optimal" or sol['status'] == "unknown":
				if sol['status'] == "optimal":
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found optimal solution.")
				else:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reached maximum num iterations.")
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tprimal objective = %f" %(sol['primal objective']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tdual objective = %f" %(sol['dual objective']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tgap = %f" %(sol['gap']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\trelative gap = %f" %(sol['relative gap']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tprimal infeasibility = %f" %(sol['primal infeasibility']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tdual infeasibility = %f" %(sol['dual infeasibility']))
				for segi in range(lseg):
					self.seq_edges[segi] += [sol['x'][segi] * 2]
					if sol['x'][segi] * 2 > self.max_CN:
						self.max_CN = sol['x'][segi] * 2
				for eci in range(lc):
					self.concordant_edges[eci] += [sol['x'][lseg + eci] * 2]
					if sol['x'][lseg + eci] * 2 > self.max_CN:
						self.max_CN = sol['x'][lseg + eci] * 2
				for edi in range(ld):
					self.discordant_edges[edi] += [sol['x'][lseg + lc + edi] * 2]
					if sol['x'][lseg + lc + edi] * 2 > self.max_CN:
						self.max_CN = sol['x'][lseg + lc + edi] * 2
				for bpi in range(lnbp):
					self.new_bp_list_[bpi] += [sol['x'][lseg + lc + ld + bpi] * 2]
					if sol['x'][lseg + lc + ld + bpi] * 2 > self.max_CN:
						self.max_CN = sol['x'][lseg + lc + ld + bpi] * 2
				for esi in range(len(self.source_edges)):
					self.source_edges[esi] += [sol['x'][lseg + lc + ld + lnbp + esi] * 2]
					if sol['x'][lseg + lc + ld + lnbp + esi] * 2 > self.max_CN:
						self.max_CN = sol['x'][lseg + lc + ld + lnbp + esi] * 2
		else:
			assert lc == 0 and ld == 0 and lnbp == 0 and lsrc == 0
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Skipped convex optimization.")
			for segi in range(lseg):
				cn_segi = (self.sr_length * sseg[3]) / (self.normal_cov_sr * sseg[-1]) 
				sseg = self.seq_edges[segi]
				if lr_seq_dist == 'poisson':
					cn_segi += (sseg[4][0] * sseg[4][1]) / (self.normal_cov_lr * sseg[-1]) 
					#print ("LR CN: %f" %((sseg[4][0] * sseg[4][1]) / (self.normal_cov_lr * sseg[-1])))
				else:
					cn_segi += self.normal_cov_lr * sseg[-1] / sseg[4][1]
				#cn_segi *= 0.5
				self.seq_edges[segi] += [cn_segi]
				if cn_segi > self.max_CN:
					self.max_CN = cn_segi
		self.max_CN += 1.0


	def valid_path(self, path):
		if len(path) < 3 or len(path) % 2 == 0:
			return False
		if path[0][0] != 's' or path[-1][0] != 's':
			return False
		for i in range(len(path)):
			if i % 2 == 0:
				if len(path[i]) != 2:
					return False
			else:
				if len(path[i]) != 3:
					return False
				e1 = path[i - 1]
				e2 = path[i + 1]
				try:
					if (e1[0] == 's' and e2[0] == 's') or (e1[0] != 's' and e2[0] != 's'):
						return False
					if e1[1] not in self.nodes[path[i]][edge_type_to_index[e1[0]]]:
						return False
					if e2[1] not in self.nodes[path[i]][edge_type_to_index[e2[0]]]:
						return False
				except:
					return False
		return True


	def alignment_to_path(self, rint):
		seq_edge_list = []
		for segi in range(len(self.seq_edges)):
			sseg = self.seq_edges[segi]
			if interval_overlap(rint, sseg):
				seq_edge_list.append(segi)
		if len(seq_edge_list) == 0:
			return []
		# potential path covers two sequence edges
		seq_edge_list = sorted(seq_edge_list, key = lambda item: self.seq_edges[item][1])
		# require the last sequence edge with sufficient coverage
		segi0 = seq_edge_list[0]
		if len(seq_edge_list) > 1 and min(self.seq_edges[segi0][2], rint[2]) - max(self.seq_edges[segi0][1], rint[1]) < 500: # need to parameterize this
			del seq_edge_list[0]
		segi0 = seq_edge_list[-1]
		if len(seq_edge_list) > 1 and min(self.seq_edges[segi0][2], rint[2]) - max(self.seq_edges[segi0][1], rint[1]) < 500: # need to parameterize this
			del seq_edge_list[-1]
		if len(seq_edge_list) <= 2:
			return []
		segi0 = seq_edge_list[0]
		node1 = (self.seq_edges[segi0][0], self.seq_edges[segi0][1], '-')
		segi0 = seq_edge_list[-1]
		node2 = (self.seq_edges[segi0][0], self.seq_edges[segi0][2], '+')
		#print (rint, node1, node2, seq_edge_list)
		path_ = self.traverse_through_sequence_edge(node1, node2)[1:-1] 
		return path_


	def chimeric_alignment_to_path_l(self, rints, ai, bp_node):
		al = rints[ai]
		seq_edge_list = []
		for segi in range(len(self.seq_edges)):
			sseg = self.seq_edges[segi]
			if al[-1] == '+':
				if interval_overlap(al, sseg):
					seq_edge_list.append([segi, '+'])
			else:
				if interval_overlap([al[0], al[2], al[1]], sseg):
					seq_edge_list.append([segi, '-'])
		if len(seq_edge_list) == 0:
			return []
		# potential path covers two sequence edges
		if seq_edge_list[0][1] == '+':
			seq_edge_list = sorted(seq_edge_list, key = lambda item: self.seq_edges[item[0]][1])
			# require the last sequence edge with sufficient coverage
			segi0 = seq_edge_list[0][0]
			if len(seq_edge_list) > 1 and min(self.seq_edges[segi0][2], al[2]) - max(self.seq_edges[segi0][1], al[1]) < 500: # need to parameterize this
				del seq_edge_list[0]
			# check if the rightmost node connects to the breakpoint edge at index edi 
			while len(seq_edge_list) > 0:
				segi_last = seq_edge_list[-1][0]
				rnode = (self.seq_edges[segi_last][0], self.seq_edges[segi_last][2], '+')
				if rnode != bp_node:
					del seq_edge_list[-1]
				else:
					break
		else:
			seq_edge_list = sorted(seq_edge_list, key = lambda item: self.seq_edges[item[0]][1], reverse = True)
			# require the last sequence edge with sufficient coverage
			segi0 = seq_edge_list[0][0]
			if len(seq_edge_list) > 1 and min(self.seq_edges[segi0][2], al[1]) - max(self.seq_edges[segi0][1], al[2]) < 500:
				del seq_edge_list[0]
			# check if the rightmost node connects to the breakpoint edge at index edi 
			while len(seq_edge_list) > 0:
				segi_last = seq_edge_list[-1][0]
				rnode = (self.seq_edges[segi_last][0], self.seq_edges[segi_last][1], '-')
				if rnode != bp_node:
					del seq_edge_list[-1]
				else:
					break
		if len(seq_edge_list) == 0:
			return []
		path_l = []
		for si in range(0, len(seq_edge_list)):
			path_l.append(('s', seq_edge_list[si][0]))
			if seq_edge_list[si][1] == '+':
				path_l.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si][0]][2], '+'))
			else:
				path_l.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si][0]][1], '-'))
			if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == '+':
				if self.seq_edges[seq_edge_list[si][0]][2] + 1 == self.seq_edges[seq_edge_list[si + 1][0]][1]:
					for ci in range(len(self.concordant_edges)):
						if self.concordant_edges[ci][0] == self.seq_edges[seq_edge_list[si][0]][0] and \
							self.seq_edges[seq_edge_list[si][0]][2] == self.concordant_edges[ci][1] and \
							self.seq_edges[seq_edge_list[si + 1][0]][1] == self.concordant_edges[ci][4]:
							path_l.append(('c', ci))
							path_l.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si + 1][0]][1], '-'))
							break
			if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == '-':
				if self.seq_edges[seq_edge_list[si][0]][1] - 1 == self.seq_edges[seq_edge_list[si + 1][0]][2]:
					for ci in range(len(self.concordant_edges)):
						if self.concordant_edges[ci][0] == self.seq_edges[seq_edge_list[si][0]][0] and \
							self.seq_edges[seq_edge_list[si + 1][0]][2] == self.concordant_edges[ci][1] and \
							self.seq_edges[seq_edge_list[si][0]][1] == self.concordant_edges[ci][4]:
							path_l.append(('c', ci))
							path_l.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si + 1][0]][2], '+'))
							break
		return path_l


	def chimeric_alignment_to_path_r(self, rints, ai, bp_node):
		ar = rints[ai]
		seq_edge_list = []
		for segi in range(len(self.seq_edges)):
			sseg = self.seq_edges[segi]
			if ar[-1] == '+':
				if interval_overlap(ar, sseg):
					seq_edge_list.append([segi, '+'])
			else:
				if interval_overlap([ar[0], ar[2], ar[1]], sseg):
					seq_edge_list.append([segi, '-'])
		if len(seq_edge_list) == 0:
			return []
		# potential path covers two sequence edges
		if seq_edge_list[0][1] == '+':
			seq_edge_list = sorted(seq_edge_list, key = lambda item: self.seq_edges[item[0]][1])
			# require the last sequence edge with sufficient coverage
			segi1 = seq_edge_list[-1][0]
			if min(self.seq_edges[segi1][2], ar[2]) - max(self.seq_edges[segi1][1], ar[1]) < 500: # need to parameterize this
				del seq_edge_list[-1]
			# check if the leftmost node connects to the breakpoint edge at index edi 
			while len(seq_edge_list) > 0:
				segi_last = seq_edge_list[0][0]
				lnode = (self.seq_edges[segi_last][0], self.seq_edges[segi_last][1], '-')
				if lnode != bp_node:
					del seq_edge_list[0]
				else:
					break
		else:
			seq_edge_list = sorted(seq_edge_list, key = lambda item: self.seq_edges[item[0]][1], reverse = True)
			#print (seq_edge_list)
			# require the last sequence edge with sufficient coverage
			segi1 = seq_edge_list[-1][0]
			if min(self.seq_edges[segi1][2], ar[1]) - max(self.seq_edges[segi1][1], ar[2]) < 500: # need to parameterize this
				del seq_edge_list[-1]
			#print (seq_edge_list)
			# check if the leftmost node connects to the breakpoint edge at index edi 
			while len(seq_edge_list) > 0:
				segi_last = seq_edge_list[0][0]
				lnode = (self.seq_edges[segi_last][0], self.seq_edges[segi_last][2], '+')
				if lnode != bp_node:
					del seq_edge_list[0]
				else:
					break
		if len(seq_edge_list) == 0:
			return []			
		path_r = []
		for si in range(0, len(seq_edge_list)):
			if seq_edge_list[si][1] == '+':
				path_r.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si][0]][1], '-'))
			else:
				path_r.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si][0]][2], '+'))
			path_r.append(('s', seq_edge_list[si][0]))
			if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == '+':
				if self.seq_edges[seq_edge_list[si][0]][2] + 1 == self.seq_edges[seq_edge_list[si + 1][0]][1]:
					for ci in range(len(self.concordant_edges)):
						if self.concordant_edges[ci][0] == self.seq_edges[seq_edge_list[si][0]][0] and \
							self.seq_edges[seq_edge_list[si][0]][2] == self.concordant_edges[ci][1] and \
							self.seq_edges[seq_edge_list[si + 1][0]][1] == self.concordant_edges[ci][4]:
							path_r.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si][0]][2], '+'))
							path_r.append(('c', ci))
							break
			if si < len(seq_edge_list) - 1 and seq_edge_list[si][1] == '-':
				if self.seq_edges[seq_edge_list[si][0]][1] - 1 == self.seq_edges[seq_edge_list[si + 1][0]][2]:
					for ci in range(len(self.concordant_edges)):
						if self.concordant_edges[ci][0] == self.seq_edges[seq_edge_list[si][0]][0] and \
							self.seq_edges[seq_edge_list[si + 1][0]][2] == self.concordant_edges[ci][1] and \
							self.seq_edges[seq_edge_list[si][0]][1] == self.concordant_edges[ci][4]:
							path_r.append((self.seq_edges[seq_edge_list[si][0]][0], self.seq_edges[seq_edge_list[si][0]][1], '-'))
							path_r.append(('c', ci))
							break
		return path_r


	def chimeric_alignment_to_path_i(self, rints, ai1, ai2, edi):
		path_ = [('d', edi)]
		node1, node2 = (), ()
		if edi < len(self.discordant_edges):
			node1 = (self.discordant_edges[edi][0], self.discordant_edges[edi][1], self.discordant_edges[edi][2])
			node2 = (self.discordant_edges[edi][3], self.discordant_edges[edi][4], self.discordant_edges[edi][5])
		else: 
			bpi = edi - len(self.discordant_edges)
			node1 = (self.new_bp_list_[bpi][0], self.new_bp_list_[bpi][1], self.new_bp_list_[bpi][2])
			node2 = (self.new_bp_list_[bpi][3], self.new_bp_list_[bpi][4], self.new_bp_list_[bpi][5])
		if ai1 > ai2:
			#node1, node2 = node2, node1
			path_ = self.chimeric_alignment_to_path_l(rints, ai2, node2) + path_ + self.chimeric_alignment_to_path_r(rints, ai1, node1)
			#print (node1, node2, path_)
		else:
			path_ = self.chimeric_alignment_to_path_l(rints, ai1, node1) + path_ + self.chimeric_alignment_to_path_r(rints, ai2, node2)
			#print (node1, node2, path_)
		return path_


	def traverse_through_sequence_edge(self, start_node, end_node):
		assert start_node[2] != end_node[2]
		path_ = [start_node]
		seqi = self.nodes[start_node][1][0]
		seq_edge = self.seq_edges[seqi]
		next_end = (seq_edge[0], seq_edge[1], '-')
		if start_node[2] == '-':
			next_end = (seq_edge[0], seq_edge[2], '+')
		path_.append(('s', seqi))
		path_.append(next_end)
		while next_end != end_node:
			ci = self.nodes[next_end][2][0]
			path_.append(('c', ci))
			cedge = self.concordant_edges[ci]
			next_start = (cedge[0], cedge[1], cedge[2])
			if next_start == next_end:
				next_start = (cedge[3], cedge[4], cedge[5])
			path_.append(next_start)
			seqi = self.nodes[next_start][1][0]
			seq_edge = self.seq_edges[seqi]
			next_end = (seq_edge[0], seq_edge[1], '-')
			if next_start[2] == '-':
				next_end = (seq_edge[0], seq_edge[2], '+')
			path_.append(('s', seqi))
			path_.append(next_end)
		return path_			


	def chimeric_alignment_to_path(self, rints, ai_list, bp_list):
		path_ = []
		lastnode = ()
		for i in range(len(bp_list)):
			node1, node2 = (), ()
			if bp_list[i] < len(self.discordant_edges):
				bpi = bp_list[i]
				node1 = (self.discordant_edges[bpi][0], self.discordant_edges[bpi][1], self.discordant_edges[bpi][2])
				node2 = (self.discordant_edges[bpi][3], self.discordant_edges[bpi][4], self.discordant_edges[bpi][5])
			else: 
				bpi = bp_list[i] - len(self.discordant_edges)
				node1 = (self.new_bp_list_[bpi][0], self.new_bp_list_[bpi][1], self.new_bp_list_[bpi][2])
				node2 = (self.new_bp_list_[bpi][3], self.new_bp_list_[bpi][4], self.new_bp_list_[bpi][5])
			if ai_list[i][0] > ai_list[i][1]:
				if i == 0:
					path_ = self.chimeric_alignment_to_path_l(rints, ai_list[i][1], node2) + [('d', bp_list[i])]
					lastnode = node1
				else:
					path_ += self.traverse_through_sequence_edge(lastnode, node2)
					path_.append(('d', bp_list[i]))
					lastnode = node1
					if i == len(bp_list) - 1:
						path_ += self.chimeric_alignment_to_path_r(rints, ai_list[i][0], node1)
			else:
				if i == 0:
					path_ = self.chimeric_alignment_to_path_l(rints, ai_list[i][0], node1) + [('d', bp_list[i])]
					lastnode = node2
				else:
					path_ += self.traverse_through_sequence_edge(lastnode, node1)
					path_.append(('d', bp_list[i]))
					lastnode = node2
					if i == len(bp_list) - 1:
						path_ += self.chimeric_alignment_to_path_r(rints, ai_list[i][1], node2)
		#print (path_)
		return path_


	#def large_del_alignment_to_path(self, rints, sai, eai, bp_list):


	def compute_path_constraints(self):
		lseg = len(self.seq_edges)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lnbp = len(self.new_bp_list_)

		bp_reads = dict()
		concordant_reads = set([])
		for edi in range(ld):
			for r_ in self.discordant_edges[edi][-2]:
				if r_[1] == r_[2]:
					if r_[0] in bp_reads:
						bp_reads[r_[0]][1].append([r_[1], r_[2], edi])
					else:
						bp_reads[r_[0]] = [[], [[r_[1], r_[2], edi]]]
				else:

					if r_[0] in bp_reads:
						bp_reads[r_[0]][0].append([r_[1], r_[2], edi])
					else:
						bp_reads[r_[0]] = [[[r_[1], r_[2], edi]], []]
		for bpi in range(lnbp):
			for r_ in self.new_bp_list_[bpi][-2]:
				if r_[1] == r_[2]:
					if r_[0] in bp_reads:
						bp_reads[r_[0]][1].append([r_[1], r_[2], ld + bpi])
					else:
						bp_reads[r_[0]] = [[], [[r_[1], r_[2], ld + bpi]]]
				else:
					if r_[0] in bp_reads:
						bp_reads[r_[0]][0].append([r_[1], r_[2], ld + bpi])
					else:
						bp_reads[r_[0]] = [[[r_[1], r_[2], ld + bpi]], []]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d reads in total covering at least one breakpoint." %len(bp_reads))

		for rn in bp_reads.keys():
			bp_reads_rn = bp_reads[rn][0]
			bp_reads_rn_sdel = bp_reads[rn][1]
			path = []
			if len(bp_reads_rn) == 1 and len(bp_reads_rn_sdel) == 0:
				rints = self.chimeric_alignments[rn][1]
				ai1 = bp_reads_rn[0][0]
				ai2 = bp_reads_rn[0][1]
				bpi = bp_reads_rn[0][2]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Read %s covers a single breakpoint." %rn)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bp = (%d, %d, %d)" \
						%(rints, self.chimeric_alignments[rn][2], ai1, ai2, bpi))
				path = self.chimeric_alignment_to_path_i(rints, ai1, ai2, bpi)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tResulting subpath = %s" %path)	
			elif len(bp_reads_rn) > 1 and len(bp_reads_rn_sdel) == 0:
				bp_reads_rn = sorted(bp_reads_rn, key = lambda item: min(item[0], item[1]))
				bp_reads_rn_split = [[0]]
				last_ai = max(bp_reads_rn[0][0], bp_reads_rn[0][1])
				for i in range(1, len(bp_reads_rn)):
					if min(bp_reads_rn[i][0], bp_reads_rn[i][1]) == last_ai:
						bp_reads_rn_split[-1].append(i)
					else:
						bp_reads_rn_split.append([i])
					last_ai = max(bp_reads_rn[i][0], bp_reads_rn[i][1])
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Read %s covers multiple breakpoints." %rn)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %bp_reads_rn_split)
				qints = self.chimeric_alignments[rn][0]
				skip = 0
				for qi in range(len(qints) - 1):
					if qints[qi + 1][0] - qints[qi][1] < -200:
						skip = 1
						break
				if skip == 1:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the read due to overlapping local alignments.")
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s." \
							%(self.chimeric_alignments[rn][1], self.chimeric_alignments[rn][2]))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on the read = %s." %qints)
					continue
				for ai_block in bp_reads_rn_split:
					rints = self.chimeric_alignments[rn][1]
					ai_list = [bp_reads_rn[bi][:2] for bi in ai_block]
					bp_list = [bp_reads_rn[bi][2] for bi in ai_block]
					if len(set(bp_list)) < len(bp_list):
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the block due to repeated breakpoints.")
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %ai_block)
						continue
					path = self.chimeric_alignment_to_path(rints, ai_list, bp_list)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bps = %s" \
						%(rints, self.chimeric_alignments[rn][2], bp_reads_rn))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tResulting subpath = %s" %path)
			elif len(bp_reads_rn) == 0 and len(bp_reads_rn_sdel) == 1:
				rints = self.large_indel_alignments[rn][0]
				rq = rints[-1]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Read %s covers a single small del breakpoint." %rn)
				if rints[3] < rints[4]:
					if rints[2] < rints[1]:
						rints = [[rints[0], rints[3], rints[2], '+'], [rints[0], rints[1], rints[4], '+']]
					else:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the read due to inconsistent alignment information.")
						continue
				else:
					if rints[2] > rints[1]:
						rints = [[rints[0], rints[3], rints[2], '-'], [rints[0], rints[1], rints[4], '-']]
					else:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the read due to inconsistent alignment information.")
						continue
				bpi = bp_reads_rn_sdel[0][2]
				if rints[0][3] == '+':
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bp = (1, 0, %d)" \
							%(rints, rq, bpi))
					path = self.chimeric_alignment_to_path_i(rints, 1, 0, bpi)
				else:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bp = (0, 1, %d)" \
							%(rints, rq, bpi))
					path = self.chimeric_alignment_to_path_i(rints, 0, 1, bpi)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tResulting subpath = %s" %path)
			elif len(bp_reads_rn) == 0 and len(bp_reads_rn_sdel) > 1:
				rints = self.large_indel_alignments[rn]
				rq = rints[0][-1]
				rints_ = set([(rint[0], min(rint[3], rint[4]), max(rint[3], rint[4])) for rint in rints])
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Read %s covers multiple small del breakpoints." %rn)
				if len(rints_) > 1 or len(rints) <= 1:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the read due to inconsistent alignment information.")
					continue
				rints_ = [[rint[0], min(rint[3], rint[4]), max(rint[3], rint[4]), '+'] for rint in rints]
				rints = sorted(rints, key = lambda item: min(item[1], item[2]))
				for ri in range(len(rints)):
					rint = rints[ri]
					rints_.append([rint[0], min(rint[3], rint[4]), max(rint[3], rint[4]), '+'])
					rints_[ri][2] = min(rint[1], rint[2])
					rints_[ri + 1][1] = max(rint[1], rint[2])
				bp_reads_rn_sdel_split = [[]]
				bp_reads_rn_sdel = sorted(bp_reads_rn_sdel, key = lambda item: item[0])
				last_ai = 0
				for i in range(len(bp_reads_rn_sdel)):
					if i == 0 or bp_reads_rn_sdel[i][0] == last_ai + 1: 
						bp_reads_rn_sdel_split[-1].append(i)
					else:
						bp_reads_rn_sdel_split.append([i])
					last_ai = bp_reads_rn_sdel[i][0]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %bp_reads_rn_sdel_split)
				for ai_block in bp_reads_rn_sdel_split:
					ai_list = [[bp_reads_rn_sdel[bi][0], bp_reads_rn_sdel[bi][0] + 1] for bi in ai_block]
					bp_list = [bp_reads_rn_sdel[bi][2] for bi in ai_block]
					if len(set(bp_list)) < len(bp_list):
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the block due to repeated breakpoints.")
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %ai_block)
						continue
					path = self.chimeric_alignment_to_path(rints_, ai_list, bp_list)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bps = %s" \
						%(rints_, rq, bp_reads_rn_sdel))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tResulting subpath = %s" %path)
			else:
				rints = self.chimeric_alignments[rn][1]
				rints_ = self.large_indel_alignments[rn]
				rint_split = []
				skip = 0
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Read %s covers breakpoints and small del breakpoints." %rn)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s" \
						%(rints, self.chimeric_alignments[rn][2]))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSmall del alignment intervals on reference = %s" %rints_)
				for rint_ in rints_:
					fount_split_rint = 0
					for ri in range(len(rints)):
						rint = rints[ri]
						if rint_[0] == rint[0] and min(rint_[1], rint_[2]) > min(rint[1], rint[2]) and max(rint_[1], rint_[2]) < max(rint[1], rint[2]):
							fount_split_rint = 1
							rint_split.append(ri)
							break
					if fount_split_rint == 0:
						skip = 1
						break
				if skip == 1:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the read due to inconsistent alignment information.")
					continue
				for rsi in range(len(rint_split)):
					ri = rint_split[rsi]
					rints.insert(ri, rints[ri][:])
					if rints[ri][3] == '+':
						rints[ri][2] = min(rints_[rsi][1], rints_[rsi][2])
						rints[ri + 1][1] = max(rints_[rsi][1], rints_[rsi][2])
					else:
						rints[ri][2] = max(rints_[rsi][1], rints_[rsi][2])
						rints[ri + 1][1] = min(rints_[rsi][1], rints_[rsi][2])
					for i in range(len(bp_reads_rn)):
						if bp_reads_rn[i][0] >= ri and bp_reads_rn[i][1] >= ri:
							bp_reads_rn[i][0] += 1
							bp_reads_rn[i][1] += 1
					for i in range(len(bp_reads_rn_sdel)):
						if bp_reads_rn_sdel[i][0] == rsi:
							if rints[ri][3] == '+':
								bp_reads_rn.append([ri + 1, ri, bp_reads_rn_sdel[i][2]])
							else:
								bp_reads_rn.append([ri, ri + 1, bp_reads_rn_sdel[i][2]])
				#print ("bp_reads_rn = ", bp_reads_rn)
				bp_reads_rn = sorted(bp_reads_rn, key = lambda item: min(item[0], item[1]))
				bp_reads_rn_split = [[0]]
				last_ai = max(bp_reads_rn[0][0], bp_reads_rn[0][1])
				for i in range(1, len(bp_reads_rn)):
					if min(bp_reads_rn[i][0], bp_reads_rn[i][1]) == last_ai:
						bp_reads_rn_split[-1].append(i)
					else:
						bp_reads_rn_split.append([i])
					last_ai = max(bp_reads_rn[i][0], bp_reads_rn[i][1])
				#print ("bp_reads_rn_split = ", bp_reads_rn_split)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %bp_reads_rn_split)
				qints = self.chimeric_alignments[rn][0]
				skip = 0
				for qi in range(len(qints) - 1):
					if qints[qi + 1][0] - qints[qi][1] < -200:
						skip = 1
						break
				if skip == 1:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the read due to overlapping local alignments.")
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s." \
							%(self.chimeric_alignments[rn][1], self.chimeric_alignments[rn][2]))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on the read = %s." %qints)
					continue
				for ai_block in bp_reads_rn_split:
					ai_list = [bp_reads_rn[bi][:2] for bi in ai_block]
					bp_list = [bp_reads_rn[bi][2] for bi in ai_block]
					if len(set(bp_list)) < len(bp_list):
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the block due to repeated breakpoints.")
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %ai_block)
						continue
					path = self.chimeric_alignment_to_path(rints, ai_list, bp_list)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq (unsplit) = %s; bps = %s" \
						%(rints, self.chimeric_alignments[rn][2], bp_reads_rn))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tResulting subpath = %s" %path)
			if len(path) > 5 and self.valid_path(path):
				if path not in self.path_constraints[0]:
					self.path_constraints[0].append(path)
					self.path_constraints[1].append(1)
				else:
					pci = self.path_constraints[0].index(path)
					self.path_constraints[1][pci] += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d distinct subpaths due to reads involving breakpoints." %len(self.path_constraints[0]))
		#extract reads in concordant_edges_reads
		for ci in range(lc):
			for rn in self.concordant_edges_reads[ci]:
				if rn not in self.large_indel_alignments and rn not in self.chimeric_alignments:
					concordant_reads.add(rn)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d concordant reads within amplicon intervals." %len(concordant_reads))
		for aint in self.amplicon_intervals:
			for read in self.lr_bamfh.fetch(aint[0], aint[1], aint[2] + 1):
				rn = read.query_name
				q = read.mapq
				if q >= 20 and rn in concordant_reads:
					path = self.alignment_to_path([read.reference_name, read.reference_start, read.reference_end])
					#print ("case c:", [read.reference_name, read.reference_start, read.reference_end], path)
					#print ()
					if len(path) > 5:
						if path not in self.path_constraints[0]:
							self.path_constraints[0].append(path)
							self.path_constraints[1].append(1)
						else:
							pci = self.path_constraints[0].index(path)
							self.path_constraints[1][pci] += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d distinct subpaths in total." %len(self.path_constraints[0]))
		#print (self.small_del_indices)
		#print (self.small_del_indices_)
		#os.abort()


	def minimize_cycles(self, k, total_weights, node_order, pc_list, pc_indices,
				max_bp_repeat = 1, p_total_weight = 0.9, p_bp_cn = 0.9, num_threads = 16):
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Regular cycle decomposition with at most %d cycles/paths allowed." %k)
		lseg = len(self.seq_edges)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lnbp = len(self.new_bp_list_)
		lsrc = len(self.source_edges)
		nnodes = len(self.nodes)
		nedges = lseg + lc + ld + lnbp + 2 * lsrc + 2 * len(self.endnodes)

		# Gurobi model
		m = gp.Model("cycle_decomposition_" + str(k))
		
		# z[i]: indicating whether cycle or path i exists
		z = m.addVars(k, vtype = GRB.BINARY, name = ["z" + str(i) for i in range(k)])
		
		# w[i]: the weight of cycle or path i, continuous variable
		w = m.addVars(k, lb = 0.0, ub = self.max_CN, vtype = GRB.CONTINUOUS, name = ["w" + str(i) for i in range(k)])
		
		# Relationship between w[i] and z[i]
		for i in range(k):
			m.addConstr(w[i] <= z[i] * self.max_CN)

		# x: the number of times an edge occur in cycle or path i
		x_names = []
		for ei in range(nedges):
			for i in range(k):
				x_names.append("x" + str(ei) + "," + str(i))
		x = m.addVars(k * nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)
	
		# Objective: minimize the total number of cycles
		obj = gp.QuadExpr(0.0)
		for i in range(k):
			obj += z[i]
			for seqi in range(lseg):
				obj -= (x[seqi * k + i] * w[i] * self.seq_edges[seqi][-2] / total_weights)
		m.setObjective(obj, GRB.MINIMIZE)

		# Must include at least 0.9 * total CN weights
		total_weights_expr = gp.QuadExpr(0.0)
		for i in range(k):
			for seqi in range(lseg):
				total_weights_expr += (x[seqi * k + i] * w[i] * self.seq_edges[seqi][-2])
		m.addConstr(total_weights_expr >= p_total_weight * total_weights)

		# Eulerian constraint
		for node in self.nodes.keys():
			if node in self.endnodes and len(self.nodes[node][0]) == 0:
				for i in range(k):
					m.addConstr(x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + i] + \
							x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + k + i] \
							== x[self.nodes[node][1][0] * k + i])
			else:
				for i in range(k):
					ec_expr = gp.LinExpr(0.0)
					for seqi in self.nodes[node][1]:
						ec_expr += x[seqi * k + i]
					for ci in self.nodes[node][2]:
						ec_expr -= x[(lseg + ci) * k + i]
					for di in self.nodes[node][3]:
						dedge = []
						if di < len(self.discordant_edges):
							dedge = self.discordant_edges[di]
						else:
							dedge = self.new_bp_list_[di - len(self.discordant_edges)]
						if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]:
							ec_expr -= 2 * x[(lseg + lc + di) * k + i]
						else:
							ec_expr -= x[(lseg + lc + di) * k + i]
					for srci in self.nodes[node][4]:
						ec_expr -= x[(lseg + lc + ld + lnbp + 2 * srci) * k + i] # connected to s
						ec_expr -= x[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i] # connected to t
					m.addConstr(ec_expr == 0.0)
		for i in range(k):
			path_expr = gp.LinExpr(0.0)
			for node in self.endnodes:
				if len(self.nodes[node][0]) == 0:
					path_expr += x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + i] # (s, v)
					path_expr -= x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + k + i] # (v, t)
			for srci in range(lsrc): 
				path_expr += x[(lseg + lc + ld + lnbp + 2 * srci) * k + i] # (s, v)
				path_expr -= x[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i] # (v, t)
			m.addConstr(path_expr == 0.0)

		# CN constraint
		for seqi in range(lseg):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[seqi * k + i]
			m.addQConstr(cn_expr <= self.seq_edges[seqi][-1])
		for ci in range(lc):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + ci) * k + i]
			m.addQConstr(cn_expr <= self.concordant_edges[ci][-1])
		for di in range(ld):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + lc + di) * k + i]
			m.addQConstr(cn_expr <= self.discordant_edges[di][-1])
			m.addQConstr(cn_expr >= p_bp_cn * self.discordant_edges[di][-1])
		for bpi in range(lnbp):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + lc + ld + bpi) * k + i]
			m.addQConstr(cn_expr <= self.new_bp_list_[bpi][-1])
			m.addQConstr(cn_expr >= p_bp_cn * self.new_bp_list_[bpi][-1])
		for srci in range(lsrc):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + lc + ld + lnbp + 2 * srci) * k + i]
				cn_expr += w[i] * x[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i]
			m.addQConstr(cn_expr <= self.source_edges[srci][-1])
			
		# Occurrence of breakpoints in each cycle/path
		for i in range(k):
			for edi in range(ld):
				if edi not in self.small_del_indices:
					m.addConstr(x[(lseg + lc + edi) * k + i] <= max_bp_repeat)
			for bpi in range(lnbp):
				if bpi not in self.small_del_indices_:
					m.addConstr(x[(lseg + lc + ld + bpi) * k + i] <= max_bp_repeat)
			
		# c: decomposition i is a cycle, and start at particular node
		c_names = []
		for ni in range(nnodes):
			for i in range(k):
				c_names.append("c" + str(ni) + "," + str(i))
		c = m.addVars(k * nnodes, vtype = GRB.BINARY, name = c_names)
		
		# Relationship between c and x
		for i in range(k):
			cycle_expr = gp.LinExpr(0.0)
			for ni in range(nnodes):
				cycle_expr += c[ni * k + i]
			for node in self.endnodes:
				if len(self.nodes[node][0]) == 0:
					cycle_expr += x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + i] # (s, v)
			for srci in range(lsrc): 
				cycle_expr += x[(lseg + lc + ld + lnbp + 2 * srci) * k + i] # (s, v)
			m.addConstr(cycle_expr <= 1.0)
			
		# d: BFS/spanning tree order of the nodes in decomposition i
		d_names = []
		for ni in range(nnodes):
			for i in range(k):
				d_names.append("d" + str(ni) + "," + str(i))
		for i in range(k):
			d_names.append("ds," + str(i))
		for i in range(k):
			d_names.append("dt," + str(i))
		d = m.addVars(k * (nnodes + 2), lb = 0.0, ub = nnodes + 2, vtype = GRB.INTEGER, name = d_names) # including s, t, at the end
			
		# y: spanning tree indicator (directed)
		y1_names = []
		y2_names = []
		for ei in range(nedges):
			for i in range(k):
				y1_names.append("y1-" + str(ei) + "," + str(i))
				y2_names.append("y2-" + str(ei) + "," + str(i))
		y1 = m.addVars(k * nedges, vtype = GRB.BINARY, name = y1_names) #small -> large
		y2 = m.addVars(k * nedges, vtype = GRB.BINARY, name = y2_names) #large -> small
			
		# Relationship between c and d
		for i in range(k):
			for ni in range(nnodes):
				m.addConstr(d[ni * k + i] >= c[ni * k + i])
		for i in range(k):
			d_expr = gp.LinExpr(0.0)
			for ni in range(nnodes):
				d_expr += c[ni * k + i]
			d_expr += d[k * nnodes + i] # d_s,i
			m.addConstr(d_expr <= 1.0)
			
		# Relationship between y and z:
		for i in range(k):
			for j in range(nedges):
				m.addConstr(y1[j * k + i] <= z[i])
				m.addConstr(y2[j * k + i] <= z[i])

		# Relationship between x, y and d
		for i in range(k * nedges):
			m.addConstr(y1[i] + y2[i] <= x[i])
		for i in range(k):
			for di in range(ld):
				dedge = self.discordant_edges[di]
				if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
					m.addConstr(y1[(lseg + lc + di) * k + i] == 0)
					m.addConstr(y2[(lseg + lc + di) * k + i] == 0)
			for bpi in range(lnbp):
				dedge = self.new_bp_list_[bpi]
				if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
					m.addConstr(y1[(lseg + lc + ld + bpi) * k + i] == 0)
					m.addConstr(y2[(lseg + lc + ld + bpi) * k + i] == 0)
		for i in range(k):
			t_expr_x = gp.LinExpr(0.0)
			t_expr_y = gp.LinExpr(0.0)
			t_expr_yd = gp.QuadExpr(0.0)
			for node in self.endnodes:
				if len(self.nodes[node][0]) == 0:
					t_expr_x += x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + k + i]
					t_expr_y += y1[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + k + i]
					t_expr_yd += y1[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + k + i] * \
							(d[k * nnodes + k + i] - d[k * node_order[node] + i]) # node -> t
					expr_x = gp.LinExpr(0.0)
					expr_y = gp.LinExpr(0.0)
					expr_xc = gp.QuadExpr(0.0)
					expr_yd = gp.QuadExpr(0.0)
					for seqi in self.nodes[node][1]:
						sseg = self.seq_edges[seqi]
						node_ = (sseg[0], sseg[1], '-')
						if node_ == node:
							node_ = (sseg[0], sseg[2], '+')
						expr_x += x[seqi * k + i]
						expr_xc += x[seqi * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[seqi * k + i]
							expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[seqi * k + i]
							expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						
					expr_x += x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + i] # from s
					expr_xc += x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + i] * c[k * node_order[node] + i]
					expr_x += x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + k + i] # to t
					expr_xc += x[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + k + i] * c[k * node_order[node] + i]
					expr_y += y1[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + i] # from s
					expr_yd += y1[(lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)) * k + i] * \
								(d[k * node_order[node] + i] - d[k * nnodes + i])
					m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
					m.addConstr(expr_y <= 1.0)
					m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
					m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)
			
			for srci in range(lsrc):
				srce = self.source_edges[srci]
				t_expr_x += x[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i]
				t_expr_y += y1[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i]
				t_expr_yd += y1[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i] * \
						(d[k * nnodes + k + i] - d[k * node_order[(srce[3], srce[4], srce[5])] + i])
			m.addConstr(t_expr_x * (nnodes + 2) >= d[k * nnodes + k + i])
			m.addConstr(t_expr_y <= 1.0)
			m.addConstr(t_expr_y * nedges * k >= t_expr_x)
			m.addConstr(t_expr_yd >= t_expr_x)
			
			for node in self.nodes.keys():
				if node not in self.endnodes or len(self.nodes[node][0]) > 0:
					expr_x = gp.LinExpr(0.0)
					expr_y = gp.LinExpr(0.0)
					expr_xc = gp.QuadExpr(0.0)
					expr_yd = gp.QuadExpr(0.0)
					for seqi in self.nodes[node][1]:
						sseg = self.seq_edges[seqi]
						node_ = (sseg[0], sseg[1], '-')
						if node_ == node:
							node_ = (sseg[0], sseg[2], '+')
						expr_x += x[seqi * k + i]
						expr_xc += x[seqi * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[seqi * k + i]
							expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[seqi * k + i]
							expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for ci in self.nodes[node][2]:
						cedge = self.concordant_edges[ci]
						node_ = (cedge[0], cedge[1], cedge[2])
						if node_ == node:
							node_ = (cedge[3], cedge[4], cedge[5])
						expr_x += x[(lseg + ci) * k + i]
						expr_xc += x[(lseg + ci) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[(lseg + ci) * k + i]
							expr_yd += y1[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[(lseg + ci) * k + i]
							expr_yd += y2[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for di in self.nodes[node][3]:
						if di < len(self.discordant_edges):
							dedge = self.discordant_edges[di]
						else:
							dedge = self.new_bp_list_[di - len(self.discordant_edges)]
						node_ = (dedge[0], dedge[1], dedge[2])
						if node_ == node:
							node_ = (dedge[3], dedge[4], dedge[5])
						expr_x += x[(lseg + lc + di) * k + i]
						expr_xc += x[(lseg + lc + di) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[(lseg + lc + di) * k + i]
							expr_yd += y1[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[(lseg + lc + di) * k + i]
							expr_yd += y2[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for srci in self.nodes[node][4]:
						expr_x += x[(lseg + lc + ld + lnbp + 2 * srci) * k + i]
						expr_x += x[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i]
						expr_xc += x[(lseg + lc + ld + lnbp + 2 * srci) * k + i] * c[k * node_order[node] + i]
						expr_xc += x[(lseg + lc + ld + lnbp + 2 * srci) * k + k + i] * c[k * node_order[node] + i]
						expr_y += y1[(lseg + lc + ld + lnbp + 2 * srci) * k + i]
						expr_yd += y1[(lseg + lc + ld + lnbp + 2 * srci) * k + i] * \
								(d[k * node_order[node] + i] - d[k * nnodes + i])
					m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
					m.addConstr(expr_y <= 1.0)
					m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
					m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)

		# Subpath constraints
		r_names = []
		for pi in range(len(pc_list)):
			for i in range(k):
				r_names.append("r" + str(pi) + "," + str(i))
		r = m.addVars(k * len(pc_list), vtype = GRB.BINARY, name = r_names)
		for pi in range(len(pc_list)):
			path_constraint_ = pc_list[pi]
			sum_r = gp.LinExpr(0.0)
			for ri in range(pi * k, (pi + 1) * k):
				sum_r += r[ri]
			m.addConstr(sum_r >= 1.0)	
			for edge in path_constraint_.keys():
				for i in range(k):
					if edge[0] == 's':
						m.addConstr(x[edge[1] * k + i] >= r[pi * k + i] * path_constraint_[edge])
					elif edge[0] == 'c':
						m.addConstr(x[(lseg + edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])
					else:
						m.addConstr(x[(lseg + lc + edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])
			
		m.setParam(GRB.Param.Threads, num_threads)
		m.setParam(GRB.Param.NonConvex, 2)
		m.setParam(GRB.Param.TimeLimit, max(3600, (ld + lnbp) * 300)) # each breakpoint edge is assigned 5 minutes 
		m.write("model.lp")
		m.optimize()
		print('MS:', m.Status)
		if m.Status == GRB.INFEASIBLE:
			print('Optimization was stopped with status %d' % m.Status)
			return m.Status
		else:
			sol_z = m.getAttr('X', z)
			sol_w = m.getAttr('X', w)
			sol_d = m.getAttr('X', d)
			sol_r = m.getAttr('X', r)
				
			total_weights_included = 0.0
			for i in range(k):
				if sol_z[i] >= 0.9:
					print ("Cycle/Path %d exists; weight = %f" %(i, sol_w[i]))
					sol_x = m.getAttr('X', x)
					sol_c = m.getAttr('X', c)
					cycle_flag = -1
					for ci in range(len(sol_c)):
						if ci % k == i and sol_c[ci] >= 0.9:
							cycle_flag = ci // k
							break
					if cycle_flag == -1:
						cycle = dict()
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if xi % k == i and sol_x[xi] >= 0.9:
								xi_ = xi // k
								x_xi = int(round(sol_x[xi]))
								if xi_ < lseg:
									cycle[('e', xi_)] = x_xi
								elif xi_ < lseg + lc:
									cycle[('c', xi_ - lseg)] = x_xi
								elif xi_ < lseg + lc + ld + lnbp:
									cycle[('d', xi_ - lseg - lc)] = x_xi
								elif xi_ < lseg + lc + ld + lnbp + 2 * lsrc:
									assert x_xi == 1
									if (xi_ - lseg - lc - ld - lnbp) % 2 == 0:
										cycle[('s', (xi_ - lseg - lc - ld - lnbp) // 2)] = 1 # source edge connected to s
									else:
										cycle[('t', (xi_ - lseg - lc - ld - lnbp - 1) // 2)] = 1 # source edge connected to t
								else:
									assert x_xi == 1
									if (xi_ - lseg - lc - ld - lnbp - 2 * lsrc) % 2 == 0:
										cycle[('ns', (xi_ - lseg - lc - ld - lnbp - 2 * lsrc) // 2)] = 1 # source edge connected to s
									else:
										cycle[('nt', (xi_ - lseg - lc - ld - lnbp - 2 * lsrc - 1) // 2)] = 1 # source edge connected to t
						for pi in range(len(pc_list)):
							if sol_r[pi * k + i] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
						self.cycles[1].append(cycle)
						self.cycle_weights[1].append(sol_w[i])
						self.path_constraints_satisfied[1].append(path_constraints_s)
					else:
						cycle = dict()
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if xi % k == i and sol_x[xi] >= 0.9:
								xi_ = xi // k
								x_xi = int(round(sol_x[xi]))
								if xi_ < lseg:
									cycle[('e', xi_)] = x_xi
								elif xi_ < lseg + lc:
									cycle[('c', xi_ - lseg)] = x_xi
								elif xi_ < lseg + lc + ld + lnbp:
									cycle[('d', xi_ - lseg - lc)] = x_xi
								else:
									print ("Cyclic path cannot connect to source nodes.")
									os.abort()
						for pi in range(len(pc_list)):
							if sol_r[pi * k + i] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
						self.cycles[0].append(cycle)
						self.cycle_weights[0].append(sol_w[i])
						self.path_constraints_satisfied[0].append(path_constraints_s)
					for seqi in range(lseg):
						total_weights_included += (sol_x[seqi * k + i] * sol_w[i] * self.seq_edges[seqi][-2])
			print ("Total weights = ", total_weights_included, total_weights)
			return m.Status
			

	def maximize_weights_greedy(self, total_weights, node_order, pc_list, pc_indices, alpha = 0.01,
				max_bp_repeat = 1, p_total_weight = 0.9, resolution = 0.1, num_threads = 16):
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Integer program too large, perform greedy cycle decomposition.")
		lseg = len(self.seq_edges)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lnbp = len(self.new_bp_list_)
		lsrc = len(self.source_edges)
		nnodes = len(self.nodes)
		nedges = lseg + lc + ld + lnbp + 2 * lsrc + 2 * len(self.endnodes)

		remaining_weights = total_weights
		unsatisfied_pc = [pc_indices[i] for i in range(len(pc_indices))]
		remaining_CN = dict()
		for segi in range(lseg):
			remaining_CN[('s', segi)] = self.seq_edges[segi][-1]
		for ci in range(lc):
			remaining_CN[('c', ci)] = self.concordant_edges[ci][-1]
		for di in range(ld):
			remaining_CN[('d', di)] = self.discordant_edges[di][-1]
		for bpi in range(lnbp):
			remaining_CN[('d', ld + bpi)] = self.new_bp_list_[bpi][-1]
		for srci in range(lsrc):
			remaining_CN[('src', srci)] = self.source_edges[srci][-1]
		print (remaining_CN)

		next_w = resolution * 1.1
		cycle_id = 0
		while next_w >= resolution and remaining_weights > (1.0 - p_total_weight) * total_weights:
			num_unsatisfied_pc = 0
			for i in range(len(pc_indices)):
				if unsatisfied_pc[i] >= 0:
					num_unsatisfied_pc += 1
			pp = 1.0
			if alpha > 0:
				pp = alpha * remaining_weights / num_unsatisfied_pc # multi - objective optimization parameter
			print ("Cycle id = ", cycle_id)
			print ("Remaining weights = ", remaining_weights, total_weights)
			print ("Num unsatisfied path constraints = ", num_unsatisfied_pc, len(pc_indices))
			print ("Path constraints factor = ", pp)

			# Gurobi model
			m = gp.Model("cycle_decomposition_greedy_" + str(cycle_id))

			# z[i]: indicating whether cycle or path i exists
			z = m.addVars(1, vtype = GRB.BINARY, name = ["z0"])

			# w[i]: the weight of cycle or path i, continuous variable
			w = m.addVars(1, lb = 0.0, ub = self.max_CN, vtype = GRB.CONTINUOUS, name = ["w0"])

			# Relationship between w[i] and z[i]
			m.addConstr(w[0] <= z[0] * self.max_CN)
			#m.addConstr(w[0] >= z[0] * resolution)

			# x: the number of times an edge occur in cycle or path i
			x_names = []
			for ei in range(nedges):
				x_names.append("x" + str(ei))
			x = m.addVars(nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)

			# Subpath constraints
			r_names = []
			for pi in range(len(pc_list)):
				r_names.append("r" + str(pi))
			r = m.addVars(len(pc_list), vtype = GRB.BINARY, name = r_names)
			
			# Objective: maximize the total weight + num subpath constraints satisfied
			obj = gp.QuadExpr(0.0)
			for seqi in range(lseg):
				obj += (x[seqi] * w[0] * self.seq_edges[seqi][-2])
			for pi in range(len(pc_list)):
				if unsatisfied_pc[pi] >= 0: 
					obj += (r[pi] * max(pp, 1.0))
			m.setObjective(obj, GRB.MAXIMIZE)

			# Eulerian constraint
			for node in self.nodes.keys():
				if node in self.endnodes and len(self.nodes[node][0]) == 0:
					m.addConstr(x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)] + \
							x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node) + 1] \
							== x[self.nodes[node][1][0]])
				else:
					ec_expr = gp.LinExpr(0.0)
					for seqi in self.nodes[node][1]:
						ec_expr += x[seqi]
					for ci in self.nodes[node][2]:
						ec_expr -= x[lseg + ci]
					for di in self.nodes[node][3]:
						dedge = []
						if di < len(self.discordant_edges):
							dedge = self.discordant_edges[di]
						else:
							dedge = self.new_bp_list_[di - len(self.discordant_edges)]
						if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]:
							ec_expr -= 2 * x[lseg + lc + di]
						else:
							ec_expr -= x[lseg + lc + di]
					for srci in self.nodes[node][4]:
						ec_expr -= x[lseg + lc + ld + lnbp + 2 * srci] # connected to s
						ec_expr -= x[lseg + lc + ld + lnbp + 2 * srci + 1] # connected to t
					m.addConstr(ec_expr == 0.0)
			path_expr = gp.LinExpr(0.0)
			for node in self.endnodes:
				if len(self.nodes[node][0]) == 0:
					path_expr += x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)] # (s, v)
					path_expr -= x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node) + 1] # (v, t)
			for srci in range(lsrc): 
				path_expr += x[lseg + lc + ld + lnbp + 2 * srci] # (s, v)
				path_expr -= x[lseg + lc + ld + lnbp + 2 * srci + 1] # (v, t)
			m.addConstr(path_expr == 0.0)

			# CN constraint
			for seqi in range(lseg):
				m.addQConstr(w[0] * x[seqi] <= remaining_CN[('s', seqi)])
			for ci in range(lc):
				m.addQConstr(w[0] * x[lseg + ci] <= remaining_CN[('c', ci)])
			for di in range(ld):
				m.addQConstr(w[0] * x[lseg + lc + di] <= remaining_CN[('d', di)])
				if self.discordant_edges[di][-1] < resolution:
					m.addConstr(x[lseg + lc + di] == 0.0)
					print ("Set coverage of bp edge at index %d to 0." %(di))
			for bpi in range(lnbp):
				m.addQConstr(w[0] * x[lseg + lc + ld + bpi] <= remaining_CN[('d', ld + bpi)])
				if self.new_bp_list_[bpi][-1] < resolution:
					m.addConstr(x[lseg + lc + ld + bpi] == 0.0)
					print ("Set coverage of bp edge at index %d to 0." %(ld + bpi))
			for srci in range(lsrc):
				cn_expr = gp.QuadExpr(0.0)
				cn_expr += w[0] * x[lseg + lc + ld + lnbp + 2 * srci]
				cn_expr += w[0] * x[lseg + lc + ld + lnbp + 2 * srci + 1]
				m.addQConstr(cn_expr <= remaining_CN[('src', srci)])
			
			# Occurrence of breakpoints in each cycle/path
			for edi in range(ld):
				if edi not in self.small_del_indices:
					m.addConstr(x[lseg + lc + edi] <= max_bp_repeat)
			for bpi in range(lnbp):
				if bpi not in self.small_del_indices_:
					m.addConstr(x[lseg + lc + ld + bpi] <= max_bp_repeat)
			
			# c: decomposition i is a cycle, and start at particular node
			c_names = []
			for ni in range(nnodes):
				c_names.append("c" + str(ni))
			c = m.addVars(nnodes, vtype = GRB.BINARY, name = c_names)

			# Relationship between c and x
			cycle_expr = gp.LinExpr(0.0)
			for ni in range(nnodes):
				cycle_expr += c[ni]
			for node in self.endnodes:
				if len(self.nodes[node][0]) == 0:
					cycle_expr += x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)] # (s, v)
			for srci in range(lsrc): 
				cycle_expr += x[lseg + lc + ld + lnbp + 2 * srci] # (s, v)
			m.addConstr(cycle_expr <= 1.0)
			
			# d: BFS/spanning tree order of the nodes in decomposition i
			d_names = []
			for ni in range(nnodes):
				d_names.append("d" + str(ni))
			d_names.append("ds")
			d_names.append("dt")
			d = m.addVars(nnodes + 2, lb = 0.0, ub = nnodes + 2, vtype = GRB.INTEGER, name = d_names) # including s, t, at the end
			
			# y: spanning tree indicator (directed)
			y1_names = []
			y2_names = []
			for ei in range(nedges):
				y1_names.append("y1-" + str(ei))
				y2_names.append("y2-" + str(ei))
			y1 = m.addVars(nedges, vtype = GRB.BINARY, name = y1_names) #small -> large
			y2 = m.addVars(nedges, vtype = GRB.BINARY, name = y2_names) #large -> small
			
			# Relationship between c and d
			for ni in range(nnodes):
				m.addConstr(d[ni] >= c[ni])
			d_expr = gp.LinExpr(0.0)
			for ni in range(nnodes):
				d_expr += c[ni]
			d_expr += d[nnodes]
			m.addConstr(d_expr <= 1.0)
			
			# Relationship between y and z:
			for ei in range(nedges):
				m.addConstr(y1[ei] <= z[0])
				m.addConstr(y2[ei] <= z[0])

			# Relationship between x, y and d
			for ei in range(nedges):
				m.addConstr(y1[ei] + y2[ei] <= x[ei])
			for di in range(ld):
				dedge = self.discordant_edges[di]
				if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
					m.addConstr(y1[lseg + lc + di] == 0)
					m.addConstr(y2[lseg + lc + di] == 0)
			for bpi in range(lnbp):
				dedge = self.new_bp_list_[bpi]
				if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
					m.addConstr(y1[lseg + lc + ld + bpi] == 0)
					m.addConstr(y2[lseg + lc + ld + bpi] == 0)
			t_expr_x = gp.LinExpr(0.0)
			t_expr_y = gp.LinExpr(0.0)
			t_expr_yd = gp.QuadExpr(0.0)
			for node in self.endnodes:
				if len(self.nodes[node][0]) == 0:
					t_expr_x += x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node) + 1]
					t_expr_y += y1[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node) + 1]
					t_expr_yd += y1[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node) + 1] * \
							(d[nnodes + 1] - d[node_order[node]]) # node -> t
					expr_x = gp.LinExpr(0.0)
					expr_y = gp.LinExpr(0.0)
					expr_xc = gp.QuadExpr(0.0)
					expr_yd = gp.QuadExpr(0.0)
					for seqi in self.nodes[node][1]:
						sseg = self.seq_edges[seqi]
						node_ = (sseg[0], sseg[1], '-')
						if node_ == node:
							node_ = (sseg[0], sseg[2], '+')
						expr_x += x[seqi]
						expr_xc += x[seqi] * c[node_order[node]]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[seqi]
							expr_yd += y1[seqi] * (d[node_order[node]] - d[node_order[node_]])
						else:
							expr_y += y2[seqi]
							expr_yd += y2[seqi] * (d[node_order[node]] - d[node_order[node_]])
						
					expr_x += x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)] # from s
					expr_xc += x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)] * c[node_order[node]]
					expr_x += x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node) + 1] # to t
					expr_xc += x[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node) + 1] * c[node_order[node]]
					expr_y += y1[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)] # from s
					expr_yd += y1[lseg + lc + ld + lnbp + 2 * lsrc + 2 * self.endnodes.index(node)] * \
								(d[node_order[node]] - d[nnodes])
					m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
					m.addConstr(expr_y <= 1.0)
					m.addConstr(expr_y * nedges + expr_xc >= expr_x)
					m.addConstr(expr_yd * nedges + expr_xc >= expr_x)
			for srci in range(lsrc):
				srce = self.source_edges[srci]
				t_expr_x += x[lseg + lc + ld + lnbp + 2 * srci + 1]
				t_expr_y += y1[lseg + lc + ld + lnbp + 2 * srci + 1]
				t_expr_yd += y1[lseg + lc + ld + lnbp + 2 * srci + 1] * \
						(d[nnodes + 1] - d[node_order[(srce[3], srce[4], srce[5])]])
			m.addConstr(t_expr_x * (nnodes + 2) >= d[nnodes + 1])
			m.addConstr(t_expr_y <= 1.0)
			m.addConstr(t_expr_y * nedges >= t_expr_x)
			m.addConstr(t_expr_yd >= t_expr_x)
			
			for node in self.nodes.keys():
				if node not in self.endnodes or len(self.nodes[node][0]) > 0:
					expr_x = gp.LinExpr(0.0)
					expr_y = gp.LinExpr(0.0)
					expr_xc = gp.QuadExpr(0.0)
					expr_yd = gp.QuadExpr(0.0)
					for seqi in self.nodes[node][1]:
						sseg = self.seq_edges[seqi]
						node_ = (sseg[0], sseg[1], '-')
						if node_ == node:
							node_ = (sseg[0], sseg[2], '+')
						expr_x += x[seqi]
						expr_xc += x[seqi] * c[node_order[node]]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[seqi]
							expr_yd += y1[seqi] * (d[node_order[node]] - d[node_order[node_]])
						else:
							expr_y += y2[seqi]
							expr_yd += y2[seqi] * (d[node_order[node]] - d[node_order[node_]])
					for ci in self.nodes[node][2]:
						cedge = self.concordant_edges[ci]
						node_ = (cedge[0], cedge[1], cedge[2])
						if node_ == node:
							node_ = (cedge[3], cedge[4], cedge[5])
						expr_x += x[lseg + ci]
						expr_xc += x[lseg + ci] * c[node_order[node]]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[lseg + ci]
							expr_yd += y1[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
						else:
							expr_y += y2[lseg + ci]
							expr_yd += y2[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
					for di in self.nodes[node][3]:
						if di < len(self.discordant_edges):
							dedge = self.discordant_edges[di]
						else:
							dedge = self.new_bp_list_[di - len(self.discordant_edges)]
						node_ = (dedge[0], dedge[1], dedge[2])
						if node_ == node:
							node_ = (dedge[3], dedge[4], dedge[5])
						expr_x += x[lseg + lc + di]
						expr_xc += x[lseg + lc + di] * c[node_order[node]]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[lseg + lc + di]
							expr_yd += y1[lseg + lc + di] * (d[node_order[node]] - d[node_order[node_]])
						else:
							expr_y += y2[lseg + lc + di]
							expr_yd += y2[lseg + lc + di] * (d[node_order[node]] - d[node_order[node_]])
					for srci in self.nodes[node][4]:
						expr_x += x[lseg + lc + ld + lnbp + 2 * srci]
						expr_x += x[lseg + lc + ld + lnbp + 2 * srci + 1]
						expr_xc += x[lseg + lc + ld + lnbp + 2 * srci] * c[node_order[node]]
						expr_xc += x[lseg + lc + ld + lnbp + 2 * srci + 1] * c[node_order[node]]
						expr_y += y1[lseg + lc + ld + lnbp + 2 * srci]
						expr_yd += y1[lseg + lc + ld + lnbp + 2 * srci] * (d[node_order[node]] - d[nnodes])
					m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
					m.addConstr(expr_y <= 1.0)
					m.addConstr(expr_y * nedges + expr_xc >= expr_x)
					m.addConstr(expr_yd * nedges + expr_xc >= expr_x)

			# Subpath constraints
			for pi in range(len(pc_list)):
				path_constraint_ = pc_list[pi]
				for edge in path_constraint_.keys():
					if edge[0] == 's':
						m.addConstr(x[edge[1]] >= r[pi] * path_constraint_[edge])
					elif edge[0] == 'c':
						m.addConstr(x[lseg + edge[1]] >= r[pi] * path_constraint_[edge])
					else:
						m.addConstr(x[lseg + lc + edge[1]] >= r[pi] * path_constraint_[edge])
			
			m.setParam(GRB.Param.Threads, num_threads)
			m.setParam(GRB.Param.NonConvex, 2)
			m.setParam(GRB.Param.TimeLimit, max(3600, (ld + lnbp) * 300)) # each breakpoint edge is assigned 5 minutes
			m.write("model" + str(cycle_id) + "_alpha=" + str(alpha) + ".lp") 
			m.optimize()
			print('MS:', m.Status)
			if m.Status == GRB.INFEASIBLE:
				print('Optimization was stopped with status %d' % m.Status)
				break
			else:
				cycle_id += 1
				sol_z = m.getAttr('X', z)
				sol_w = m.getAttr('X', w)
				sol_d = m.getAttr('X', d)
				sol_r = m.getAttr('X', r)
				total_weights_included = 0.0
				if sol_z[0] >= 0.9:
					print ("Next cycle/Path exists; weight = %f" %(sol_w[0]))
					next_w = sol_w[0]
					sol_x = m.getAttr('X', x)
					sol_c = m.getAttr('X', c)
					cycle_flag = -1
					for ci in range(len(sol_c)):
						if sol_c[ci] >= 0.9:
							cycle_flag = ci
							break
					if cycle_flag == -1:
						cycle = dict()
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if sol_x[xi] >= 0.9:
								x_xi = int(round(sol_x[xi]))
								if xi < lseg:
									cycle[('e', xi)] = x_xi
									print (cycle_id, 'path', 'seq', x_xi, self.seq_edges[xi][:3])
									remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('s', xi)] < resolution:
										remaining_CN[('s', xi)] = 0.0
								elif xi < lseg + lc:
									cycle[('c', xi - lseg)] = x_xi
									print (cycle_id, 'path', 'concordant', x_xi, self.concordant_edges[xi - lseg][:6])
									remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('c', xi - lseg)] < resolution:
										remaining_CN[('c', xi - lseg)] = 0.0
								elif xi < lseg + lc + ld + lnbp:
									cycle[('d', xi - lseg - lc)] = x_xi
									if xi < lseg + lc + ld:
										print (cycle_id, 'path', 'discordant', x_xi, self.discordant_edges[xi - lseg - lc][:6])
									else:
										print (cycle_id, 'path', 'discordant', x_xi, self.new_bp_list_[xi - lseg - lc - ld][:6])
									remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('d', xi - lseg - lc)] < resolution:
										remaining_CN[('d', xi - lseg - lc)] = 0.0
								elif xi < lseg + lc + ld + lnbp + 2 * lsrc:
									assert x_xi == 1
									if (xi - lseg - lc - ld - lnbp) % 2 == 0:
										cycle[('s', (xi - lseg - lc - ld - lnbp) // 2)] = 1 # source edge connected to s
										print (cycle_id, 'path', 'source', x_xi, self.source_edges[(xi - lseg - lc - ld - lnbp) // 2][:6])
										remaining_CN[('src', (xi - lseg - lc - ld - lnbp) // 2)] -= sol_w[0]
										if remaining_CN[('src', (xi - lseg - lc - ld - lnbp) // 2)] < resolution:
											remaining_CN[('src', (xi - lseg - lc - ld - lnbp) // 2)] = 0.0
									else:
										cycle[('t', (xi - lseg - lc - ld - lnbp - 1) // 2)] = 1 # source edge connected to t
										print (cycle_id, 'path', 'source', x_xi, self.source_edges[(xi - lseg - lc - ld - lnbp - 1) // 2][:6])
										remaining_CN[('src', (xi - lseg - lc - ld - lnbp - 1) // 2)] -= sol_w[0]
										if remaining_CN[('src', (xi - lseg - lc - ld - lnbp - 1) // 2)] < resolution:
											remaining_CN[('src', (xi - lseg - lc - ld - lnbp - 1) // 2)] = 0.0
								else:
									assert x_xi == 1
									if (xi - lseg - lc - ld - lnbp - 2 * lsrc) % 2 == 0:
										cycle[('ns', (xi - lseg - lc - ld - lnbp - 2 * lsrc) // 2)] = 1 # source edge connected to s
										print (cycle_id, 'path', 'source', x_xi, self.endnodes[(xi - lseg - lc - ld - lnbp - 2 * lsrc) // 2])
									else:
										cycle[('nt', (xi - lseg - lc - ld - lnbp - 2 * lsrc - 1) // 2)] = 1 # source edge connected to t
										print (cycle_id, 'path', 'source', x_xi, self.endnodes[(xi - lseg - lc - ld - lnbp - 2 * lsrc - 1) // 2])
						for pi in range(len(pc_list)):
							if sol_r[pi] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
								unsatisfied_pc[pi] = -1
						self.cycles[1].append(cycle)
						self.cycle_weights[1].append(sol_w[0])
						self.path_constraints_satisfied[1].append(path_constraints_s)
					else:
						cycle = dict()
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if sol_x[xi] >= 0.9:
								x_xi = int(round(sol_x[xi]))
								if xi < lseg:
									cycle[('e', xi)] = x_xi
									print (cycle_id, 'cycle', 'seq', x_xi, self.seq_edges[xi][:3])
									remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('s', xi)] < resolution:
										remaining_CN[('s', xi)] = 0.0
								elif xi < lseg + lc:
									cycle[('c', xi - lseg)] = x_xi
									print (cycle_id, 'cycle', 'concordant', x_xi, self.concordant_edges[xi - lseg][:6])
									remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('c', xi - lseg)] < resolution:
										remaining_CN[('c', xi - lseg)] = 0.0
								elif xi < lseg + lc + ld + lnbp:
									cycle[('d', xi - lseg - lc)] = x_xi
									if xi < lseg + lc + ld:
										print (cycle_id, 'cycle', 'discordant', x_xi, self.discordant_edges[xi - lseg - lc][:6])
									else:
										print (cycle_id, 'cycle', 'discordant', x_xi, self.new_bp_list_[xi - lseg - lc - ld][:6])
									remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('d', xi - lseg - lc)] < resolution:
										remaining_CN[('d', xi - lseg - lc)] = 0.0
								else:
									print ("Cyclic path cannot connect to source nodes.")
									os.abort()
						for pi in range(len(pc_list)):
							if sol_r[pi] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
								unsatisfied_pc[pi] = -1
						self.cycles[0].append(cycle)
						self.cycle_weights[0].append(sol_w[0])
						self.path_constraints_satisfied[0].append(path_constraints_s)
					for seqi in range(lseg):
						total_weights_included += (sol_x[seqi] * sol_w[0] * self.seq_edges[seqi][-2])
					print ("Total weights = ", total_weights_included, total_weights)
					if total_weights_included < 0.01 * total_weights:
						break
					remaining_weights -= total_weights_included
				else:
					break


	def cycle_decomposition(self, alpha = 0.01, max_bp_repeat = 1, p_total_weight = 0.9, resolution = 0.1):
		lseg = len(self.seq_edges)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lnbp = len(self.new_bp_list_)
		lsrc = len(self.source_edges)
		total_weights = 0.0
		for segi in range(lseg):
			sseg = self.seq_edges[segi]
			total_weights += sseg[-2] * sseg[-1]
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Begin cycle decomposition.")
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total CN weights = %f." %total_weights)

		"""
		# endnodes = [] # Nodes corresponding to interval ends
		for segi in range(lseg):
			sseg = self.seq_edges[segi]
			total_weights += sseg[-2] * sseg[-1]
			if segi == 0:
				endnodes.append((sseg[0], sseg[1], '-'))
			else:
				lastseg = self.seq_edges[segi - 1]
				if lastseg[0] != sseg[0] or lastseg[2] + 1 != sseg[1]:
					endnodes.append((lastseg[0], lastseg[2], '+'))
					endnodes.append((sseg[0], sseg[1], '-'))
				if segi == lseg - 1:
					endnodes.append((sseg[0], sseg[2], '+'))
		if lseg == 1:
			sseg = self.seq_edges[0]
			endnodes.append((sseg[0], sseg[2], '+'))
		"""
		#assert len(endnodes) == len(self.amplicon_intervals) * 2
		#print (len(self.path_constraints[0]))
		len_path_constraints_ = len(self.path_constraints[0])
		for pathi in range(len_path_constraints_):
			path = self.path_constraints[0][pathi]
			path_constraint = dict()
			for ei in range(len(path)):
				if ei % 2 == 0:
					try:
						path_constraint[path[ei]] += 1
					except:
						path_constraint[path[ei]] = 1
			self.valid_path_constraints[0].append(path_constraint)
			self.valid_path_constraints[1].append(pathi)
			self.valid_path_constraints[2].append(self.path_constraints[1][pathi])
		for pathi in range(len_path_constraints_)[::-1]:
			path_constraint = self.valid_path_constraints[0][pathi]
			subpath_flag = -1
			for pathi_ in range(len(self.valid_path_constraints[0])):
				path_constraint_ = self.valid_path_constraints[0][pathi_]
				s1 = 1
				for edge in path_constraint.keys():
					if edge not in path_constraint_:
						s1 = 0
						break
					elif path_constraint_[edge] < path_constraint[edge]:
						s1 = 0
						break
				if s1 == 1 and pathi_ != pathi:
					subpath_flag = pathi_
					break
			if subpath_flag >= 0:
				del self.valid_path_constraints[0][pathi]
				del self.valid_path_constraints[1][pathi]
				self.valid_path_constraints[2][subpath_flag] = max(self.valid_path_constraints[2][subpath_flag], self.valid_path_constraints[2][pathi])
				del self.valid_path_constraints[2][pathi]

		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total num subpath constraints = %d." %len(self.valid_path_constraints[0]))
		for pathi in self.valid_path_constraints[1]:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSubpath constraints %d = %s" %(pathi, self.path_constraints[0][pathi]))

		k = max(10, (ld + lnbp) // 2) # Initial num cycles/paths
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total num initial cycles / paths = %d." %k)
		nnodes = len(self.nodes) # Does not include s and t
		node_order = dict()
		ni_ = 0
		for node in self.nodes.keys():
			node_order[node] = ni_
			ni_ += 1
		nedges = lseg + lc + ld + lnbp + 2 * lsrc + 2 * len(self.endnodes)
		if nedges < k:
			k = nedges
			logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset num cycles/paths to %d." %k)
		while k <= nedges:
			if nedges > 100 or (3 * k + 3 * k * nedges + 2 * k * nnodes + k * len(self.valid_path_constraints[0])) >= 10000:
				self.maximize_weights_greedy(total_weights, node_order, self.valid_path_constraints[0], \
						self.valid_path_constraints[1], alpha, max_bp_repeat, p_total_weight, resolution, 16)
				break 
			else:
				if self.minimize_cycles(k, total_weights, node_order, self.valid_path_constraints[0], self.valid_path_constraints[1], \
						max_bp_repeat, p_total_weight, 0.9, 16) == GRB.INFEASIBLE:
					k *= 2
				else:
					break


	def eulerian_cycle_t(self, edges_next_cycle, path_constraints_next_cycle):
		lseg = len(self.seq_edges)
		ld = len(self.discordant_edges)

		eulerian_cycle = [] # A cycle is edge - node list starting and ending with the same edge
					# Since Eulerian, there could be subcycles in the middle of a cycle
		eulerian_cycle_ = [] # Cycle in AA cycle format
		valid = 0
		num_trials = 0
		while valid == 0 and num_trials < 2000:
			valid = 1
			num_trials += 1
			eulerian_cycle = []
			eulerian_cycle_ = []
			edges_cur = edges_next_cycle.copy()
			path_constraints_cur = path_constraints_next_cycle.copy()
			last_seq_edge = lseg # Start with the edge with smallest index and on the positive strand
			for edge in edges_cur.keys():
				if edge[0] == 'e':
					last_seq_edge = min(last_seq_edge, edge[1])
			last_edge_dir = '+'
			eulerian_cycle.append(('s', last_seq_edge))
			eulerian_cycle_.append(str(last_seq_edge + 1) + '+')
			while len(edges_cur) > 0:
				seq_edge = self.seq_edges[last_seq_edge]
				node = (seq_edge[0], seq_edge[2], '+')
				if last_edge_dir == '-':
					node = (seq_edge[0], seq_edge[1], '-')
				eulerian_cycle.append(node)
				next_bp_edges = [] # Since cycle, only consider discordant edges and concordant edges
				for ci in self.nodes[node][2]:
					next_bp_edges.append(('c', ci))
				for di in self.nodes[node][3]:
					next_bp_edges.append(('d', di))
				del_list = [i for i in range(len(next_bp_edges)) if next_bp_edges[i] not in edges_cur]
				for i in del_list[::-1]:
					del next_bp_edges[i]
				if len(next_bp_edges) == 0:
					valid = 0
					break
				if len(next_bp_edges) == 1: # No branching on the path
					eulerian_cycle.append(next_bp_edges[0])
					edges_cur[next_bp_edges[0]] = int(edges_cur[next_bp_edges[0]]) - 1
					if edges_cur[next_bp_edges[0]] == 0:
						del edges_cur[next_bp_edges[0]]
						
					bp_edge = []
					if next_bp_edges[0][0] == 'c':
						bp_edge = self.concordant_edges[next_bp_edges[0][1]][:6]
					elif next_bp_edges[0][1] < ld:
						bp_edge = self.discordant_edges[next_bp_edges[0][1]][:6]
					else:
						bp_edge = self.new_bp_list_[next_bp_edges[0][1] - ld][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_cycle.append(node_)
					last_seq_edge = self.nodes[node_][1][0]
					eulerian_cycle.append(('s', last_seq_edge))
					#print ('cycle', eulerian_path, edges_cur)
					if node_[2] == '-':
						last_edge_dir = '+'
						eulerian_cycle_.append(str(last_seq_edge + 1) + '+')
					else:
						last_edge_dir = '-'
						eulerian_cycle_.append(str(last_seq_edge + 1) + '-')
					edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
					if edges_cur[('e', last_seq_edge)] == 0:
						del edges_cur[('e', last_seq_edge)]
				else:
					r = random.randint(0, len(next_bp_edges) - 1)
					eulerian_cycle.append(next_bp_edges[r])
					edges_cur[next_bp_edges[r]] = int(edges_cur[next_bp_edges[r]]) - 1
					if edges_cur[next_bp_edges[r]] == 0:
						del edges_cur[next_bp_edges[r]]
					bp_edge = []
					if next_bp_edges[r][0] == 'c':
						bp_edge = self.concordant_edges[next_bp_edges[r][1]][:6]
					elif next_bp_edges[r][1] < ld:
						bp_edge = self.discordant_edges[next_bp_edges[r][1]][:6]
					else:
						bp_edge = self.new_bp_list_[next_bp_edges[r][1] - ld][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_cycle.append(node_)
					last_seq_edge = self.nodes[node_][1][0]
					eulerian_cycle.append(('s', last_seq_edge))
					if node_[2] == '-':
						last_edge_dir = '+'
						eulerian_cycle_.append(str(last_seq_edge + 1) + '+')
					else:
						last_edge_dir = '-'
						eulerian_cycle_.append(str(last_seq_edge + 1) + '-')
					edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
					if edges_cur[('e', last_seq_edge)] == 0:
						del edges_cur[('e', last_seq_edge)]
			# check if the remaining path constraints are satisfied
			for pathi in path_constraints_cur:
				path_ = self.path_constraints[0][pathi]
				path0 = path_[0]
				s = 0
				for ei in range(len(eulerian_cycle) - 1):
					obj = eulerian_cycle[ei]
					if obj == path0:
						s_ = 1
						for i in range(len(path_)):
							if eulerian_cycle[:-1][(ei + i) % (len(eulerian_cycle) - 1)] != path_[i]:
								s_ = 0 
								break
						if s_ == 1:
							s = 1
							break
						else:
							s_ = 1
							for i in range(len(path_)):
								if eulerian_cycle[:-1][ei - i] != path_[i]:
									s_ = 0 
									break
							if s_ == 1:
								s = 1
								break
				if s == 0 and valid == 1 and num_trials < 1000:
					valid = 0
					print('not satisfied', num_trials, eulerian_cycle, path_)
		return eulerian_cycle_
				

	def eulerian_path_t(self, edges_next_path, path_constraints_next_path):
		lseg = len(self.seq_edges)
		ld = len(self.discordant_edges)

		eulerian_path = [] # A path is edge - node list starting and ending with edges
					# Since Eulerian, there could be subcycles in the middle of a path
		eulerian_path_ = [] # Path in AA cycle format
		valid = 0
		num_trials = 0
		while valid == 0 and num_trials < 2000:
			valid = 1
			num_trials += 1
			eulerian_path = [] 
			eulerian_path_ = []
			edges_cur = edges_next_path.copy()
			path_constraints_cur = path_constraints_next_path.copy()
			src_edge = ()
			last_seq_edge = lseg
			last_edge_dir = '+'
			for edge in edges_cur.keys(): #Start with the edge with smallest index
				if (edge[0] == 's' or edge[0] == 't'):
					src_edge = edge
					node = (self.source_edges[edge[1]][3], self.source_edges[edge[1]][4], self.source_edges[edge[1]][5])
					if len(eulerian_path) == 0:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path.append(('$', -1))
						eulerian_path.append(node)
						last_seq_edge = self.nodes[node][1][0]
					elif self.nodes[node][1][0] < last_seq_edge:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path[-1] = node
						last_seq_edge = self.nodes[node][1][0]
				elif (edge[0] == 'ns' or edge[0] == 'nt'):
					src_edge = edge
					node = self.endnodes[edge[1]]
					if len(eulerian_path) == 0:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path.append(('$', -1))
						eulerian_path.append(node)
						last_seq_edge = self.nodes[node][1][0]
					elif self.nodes[node][1][0] < last_seq_edge:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path[-1] = node
						last_seq_edge = self.nodes[node][1][0]					
			del edges_cur[src_edge]
			eulerian_path.append(('s', last_seq_edge))
			if last_edge_dir == '+':
				eulerian_path_.append(str(last_seq_edge + 1) + '+')
			else:
				eulerian_path_.append(str(last_seq_edge + 1) + '-')
			edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
			if edges_cur[('e', last_seq_edge)] == 0:
				del edges_cur[('e', last_seq_edge)]
			while len(edges_cur) > 0:
				seq_edge = self.seq_edges[last_seq_edge]
				node = (seq_edge[0], seq_edge[2], '+')
				if last_edge_dir == '-':
					node = (seq_edge[0], seq_edge[1], '-')
				eulerian_path.append(node)
				if len(edges_cur) == 1 and (list(edges_cur.keys())[0][0] == 's' or list(edges_cur.keys())[0][0] == 'ns' or \
					list(edges_cur.keys())[0][0] == 't' or list(edges_cur.keys())[0][0] == 'nt'):
					eulerian_path.append(('$', -1))
					break
				next_bp_edges = [] # Since cycle, only consider discordant edges and concordant edges
				for ci in self.nodes[node][2]:
					next_bp_edges.append(('c', ci))
				for di in self.nodes[node][3]:
					next_bp_edges.append(('d', di))
				del_list = [i for i in range(len(next_bp_edges)) if next_bp_edges[i] not in edges_cur]
				for i in del_list[::-1]:
					del next_bp_edges[i]
				if len(next_bp_edges) == 0:
					valid = 0
					break
				if len(next_bp_edges) == 1: # No branching on the path
					eulerian_path.append(next_bp_edges[0])
					edges_cur[next_bp_edges[0]] = int(edges_cur[next_bp_edges[0]]) - 1
					if edges_cur[next_bp_edges[0]] == 0:
						del edges_cur[next_bp_edges[0]]
					bp_edge = []
					if next_bp_edges[0][0] == 'c':
						bp_edge = self.concordant_edges[next_bp_edges[0][1]][:6]
					elif next_bp_edges[0][1] < ld:
						bp_edge = self.discordant_edges[next_bp_edges[0][1]][:6]
					else:
						bp_edge = self.new_bp_list_[next_bp_edges[0][1] - ld][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_path.append(node_)
					last_seq_edge = self.nodes[node_][1][0]
					eulerian_path.append(('s', last_seq_edge))
					if node_[2] == '-':
						last_edge_dir = '+'
						eulerian_path_.append(str(last_seq_edge + 1) + '+')
					else:
						last_edge_dir = '-'
						eulerian_path_.append(str(last_seq_edge + 1) + '-')
					edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
					if edges_cur[('e', last_seq_edge)] == 0:
						del edges_cur[('e', last_seq_edge)]
				else:
					r = random.randint(0, len(next_bp_edges) - 1)
					eulerian_path.append(next_bp_edges[r])
					edges_cur[next_bp_edges[r]] = int(edges_cur[next_bp_edges[r]]) - 1
					if edges_cur[next_bp_edges[r]] == 0:
						del edges_cur[next_bp_edges[r]]
					bp_edge = []
					if next_bp_edges[r][0] == 'c':
						bp_edge = self.concordant_edges[next_bp_edges[r][1]][:6]
					elif next_bp_edges[r][1] < ld:
						bp_edge = self.discordant_edges[next_bp_edges[r][1]][:6]
					else:
						bp_edge = self.new_bp_list_[next_bp_edges[r][1] - ld][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_path.append(node_)
					last_seq_edge = self.nodes[node_][1][0]
					eulerian_path.append(('s', last_seq_edge))
					if node_[2] == '-':
						last_edge_dir = '+'
						eulerian_path_.append(str(last_seq_edge + 1) + '+')
					else:
						last_edge_dir = '-'
						eulerian_path_.append(str(last_seq_edge + 1) + '-')
					edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
					if edges_cur[('e', last_seq_edge)] == 0:
						del edges_cur[('e', last_seq_edge)]
			# check if the remaining path constraints are satisfied
			for pathi in path_constraints_cur:
				path_ = self.path_constraints[0][pathi]
				s = 0
				for ei in range(2, len(eulerian_path) - 1 - len(path_)):
					if eulerian_path[ei: ei + len(path_)] == path_[:] or eulerian_path[ei: ei + len(path_)] == path_[::-1]:
						s = 1
						break
				if s == 0 and valid == 1 and num_trials < 1000:
					valid = 0
					print('not satisfied', num_trials, eulerian_path, path_)
		return eulerian_path_


	def output_cycles(self, cycle_file):
		lseg = len(self.seq_edges)

		fp = open(cycle_file, 'w')
		interval_num = 1
		for segi in range(lseg):
			sseg = self.seq_edges[segi]
			if segi == 0:
				fp.write("Interval\t%d\t%s\t%d\t" %(interval_num, sseg[0], sseg[1]))
			else:
				lastseg = self.seq_edges[segi - 1]
				if lastseg[0] != sseg[0] or lastseg[2] + 1 != sseg[1]:
					fp.write("%d\n" %(lastseg[2]))
					interval_num += 1
					fp.write("Interval\t%d\t%s\t%d\t" %(interval_num, sseg[0], sseg[1]))
				if segi == lseg - 1:
					fp.write("%d\n" %(sseg[2]))
		if lseg == 1:
			sseg = self.seq_edges[0]
			fp.write("%d\n" %(sseg[2]))
		fp.write("List of cycle segments\n")
		for segi in range(lseg):
			sseg = self.seq_edges[segi]
			fp.write("Segment\t%d\t%s\t%d\t%d\n" %(segi + 1, sseg[0], sseg[1], sseg[2]))
		fp.write("List of subpath constraints\n")
		path_constraint_indices_ = []
		for paths in (self.path_constraints_satisfied[0] + self.path_constraints_satisfied[1]):
			for pathi in paths:
				if pathi not in path_constraint_indices_:
					path_constraint_indices_.append(pathi)
		for constraint_i in range(len(self.valid_path_constraints[1])):
			fp.write("Path constraint\t%d\t" %(constraint_i + 1))
			pathi = self.valid_path_constraints[1][constraint_i]
			path_ = self.path_constraints[0][pathi]
			if path_[0][1] > path_[-1][1]:
				path_ = path_[::-1]
			for i in range(len(path_)):
				if i % 4 == 0:
					if i < len(path_) - 1:
						if path_[i + 1][2] == '+':
							fp.write("%d+," %(path_[i][1] + 1))
						else:
							fp.write("%d-," %(path_[i][1] + 1))
					else:
						if path_[i - 1][2] == '+':
							fp.write("%d-\t" %(path_[i][1] + 1))
						else:
							fp.write("%d+\t" %(path_[i][1] + 1))
			fp.write("Support=%d\t" %(self.valid_path_constraints[2][constraint_i]))
			if pathi in path_constraint_indices_:
				fp.write("Satisfied\n")
			else:
				fp.write("Unsatisfied\n")
		# sort cycles according to weights
		cycle_indices = sorted([(0, i) for i in range(len(self.cycle_weights[0]))] + [(1, i) for i in range(len(self.cycle_weights[1]))], 
				key = lambda item: self.cycle_weights[item[0]][item[1]], reverse = True)
		for cycle_i in cycle_indices: 
			cycle_edge_list = []
			if cycle_i[0] == 0: # cycles
				cycle_seg_list = self.eulerian_cycle_t(self.cycles[cycle_i[0]][cycle_i[1]], self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]])
				assert cycle_seg_list[0] == cycle_seg_list[-1]
				fp.write("Cycle=%d;" %(cycle_indices.index(cycle_i) + 1))
				fp.write("Copy_count=%s;" %str(self.cycle_weights[cycle_i[0]][cycle_i[1]]))
				fp.write("Segments=")
				for segi in range(len(cycle_seg_list) - 2):
					fp.write("%s," %cycle_seg_list[segi])
				fp.write("%s;" %cycle_seg_list[-2])
				fp.write("Path_constraints_satisfied=")
				for pathi in range(len(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]]) - 1):
					fp.write("%d," %(self.valid_path_constraints[1].index(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]][pathi]) + 1))
				if len(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]]) > 0:
					fp.write("%d\n" %(self.valid_path_constraints[1].index(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]][-1]) + 1))
				else:
					fp.write("\n")
			else: # paths
				cycle_seg_list = self.eulerian_path_t(self.cycles[cycle_i[0]][cycle_i[1]], self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]])
				fp.write("Cycle=%d;" %(cycle_indices.index(cycle_i) + 1))
				fp.write("Copy_count=%s;" %str(self.cycle_weights[cycle_i[0]][cycle_i[1]]))
				fp.write("Segments=0+,")
				for segi in range(len(cycle_seg_list) - 1):
					fp.write("%s," %cycle_seg_list[segi])
				fp.write("%s,0-;" %cycle_seg_list[-1])
				fp.write("Path_constraints_satisfied=")
				for pathi in range(len(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]]) - 1):
					fp.write("%d," %(self.valid_path_constraints[1].index(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]][pathi]) + 1))
				if len(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]]) > 0:
					fp.write("%d\n" %(self.valid_path_constraints[1].index(self.path_constraints_satisfied[cycle_i[0]][cycle_i[1]][-1]) + 1))
				else:
					fp.write("\n")
		fp.close()



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("--sr_bam", help = "Sorted indexed short read bam file.", required = True)
	parser.add_argument("--lr_bam", help = "Sorted indexed long read bam file.", required = True)
	parser.add_argument("--aa_graph", help = "AA-formatted graph file.", required = True)
	parser.add_argument("--aa_cycle", help = "AA-formatted cycle file.", required = True)
	parser.add_argument("--output_bp", help = "If specified, only output the list of breakpoints.", action = 'store_true')
	parser.add_argument("--output_cycle_fn", help = ".")
	parser.add_argument("--sr_cnseg", help = "Short read *.cns file.")
	parser.add_argument("--lr_cnseg", help = "Long read *.cns file.")
	parser.add_argument("--lr_seq_normal", help = "Use normal distribution on the total number of nucleotides on sequence edges for long read.", action = 'store_true')
	parser.add_argument("--aa_downsampled", help = "AA breakpoint graph was constructed with downsampled bam.", action = 'store_true')
	parser.add_argument("--ilp_alpha", help = ".", type = float)
	parser.add_argument("--log_fn", help = ".")
	args = parser.parse_args()

	log_fn = "refine_breakpoint_graph.log"
	if args.log_fn:
		log_fn = args.log_fn
	logging.basicConfig(filename = log_fn, filemode = 'w', level = logging.DEBUG, 
			format = '[%(name)s:%(levelname)s]\t%(message)s')
	logging.info("Python version " + sys.version + "\n")
	commandstring = 'Commandline: '
	for arg in sys.argv:
		if ' ' in arg:
			commandstring += '"{}" '.format(arg)
		else:
			commandstring += "{} ".format(arg)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + commandstring)

	b2bn = bam_to_breakpoint_hybrid(args.sr_bam, args.lr_bam, args.aa_cycle)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Opened SR and LR bam files.")
	b2bn.read_cycle(args.aa_cycle)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed parsing AA cycle file.")
	b2bn.read_graph(args.aa_graph)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed parsing AA graph file.")
	
	b2bn.fetch()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed fetching reads containing breakpoints.")
	if args.sr_cnseg is not None and args.sr_cnseg is not None:
		b2bn.read_cns(args.sr_cnseg, args.lr_cnseg)
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed parsing CN segment files.")
	b2bn.find_smalldel_breakpoints()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed finding small del breakpoints.")
	b2bn.find_breakpoints()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed finding all breakpoints.")
	if args.aa_downsampled:
		b2bn.del_sr_bp(aa_downsampled = True)
	else:
		b2bn.del_sr_bp(aa_downsampled = False)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Deleted SR only breakpoints with low cov or CN.")
	b2bn.split_seg_bp() # Split the sequence edges with strong long read breakpoint(s) in the middle
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed breakpoint graph structure construction.")
	if args.output_bp:
		b2bn.output_breakpoint_info(args.aa_graph.split('/')[-1][:-9] + 'breakpoints.tsv')
	else:
		if args.sr_cnseg is None or args.sr_cnseg is None:
			print("Please specify the copy number segment files.")
			os.abort()
		if args.lr_seq_normal:
			b2bn.assign_cov(lr_seq_dist = 'normal')
		else:
			b2bn.assign_cov(lr_seq_dist = 'poisson')
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Assigned read coverage for all sequence and concordant edges.")
		if args.lr_seq_normal:
			if args.aa_downsampled:
				b2bn.assign_cn(lr_seq_dist = 'normal', aa_downsampled = True)
			else:
				b2bn.assign_cn(lr_seq_dist = 'normal', aa_downsampled = False)
		else:
			if args.aa_downsampled:
				b2bn.assign_cn(lr_seq_dist = 'poisson', aa_downsampled = True)
			else:
				b2bn.assign_cn(lr_seq_dist = 'poisson', aa_downsampled = False)
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Assigned CN for all edges.")
		b2bn.compute_path_constraints()
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Computed all subpath constraints.")
		if args.ilp_alpha:
			b2bn.cycle_decomposition(alpha = args.ilp_alpha)
		else:
			b2bn.cycle_decomposition()
		b2bn.output_breakpoint_graph(args.aa_graph.split('/')[-1][:-4] + '_.txt')
		if args.output_cycle_fn:
			b2bn.output_cycles(args.output_cycle_fn)
		else:
			b2bn.output_cycles(args.aa_cycle.split('/')[-1][:-4] + '_.txt')
	b2bn.closebam()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total runtime.")

