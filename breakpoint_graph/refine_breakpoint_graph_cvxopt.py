import time
start_time = time.time()
import logging
import math
import pysam
import argparse
import sys
import os
import numpy as np
import cvxopt
import cvxopt.modeling

import chimeric_reads


neg_plus_minus = {"+": "-", "-": "+"}


chr_idx = {'chr1': 0, 'chr2': 1, 'chr3': 2, 'chr4': 3,
		'chr5': 4, 'chr6': 5, 'chr7': 6, 'chr8': 7,
		'chr9': 8, 'chr10': 9, 'chr11': 10, 'chr12': 11,
		'chr13': 12, 'chr14': 13, 'chr15': 14, 'chr16': 15,
		'chr17': 16, 'chr18': 17, 'chr19': 18, 'chr20': 19,
		'chr21': 20, 'chr22': 21, 'chrX': 22, 'chrY': 23, 'chrM': 24}


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


def interval2bp(R1, R2, rn = '', rgap = 0):
	"""
	Convert split/chimeric alignment to breakpoint
	"""
	if (chr_idx[R2[0]] < chr_idx[R1[0]]) or (chr_idx[R2[0]] == chr_idx[R1[0]] and R2[1] < R1[2]):
		return [R1[0], R1[2], R1[3], R2[0], R2[1], neg_plus_minus[R2[3]], rn, rgap, 0]
	return [R2[0], R2[1], neg_plus_minus[R2[3]], R1[0], R1[2], R1[3], rn, rgap, 1]


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
	for bp_ in bp_cluster:
		if bp_[1] <= bp_stats[0] + 3 * bp_stats[2] and bp_[1] >= bp_stats[0] - 3 * bp_stats[2] and \
			bp_[4] <= bp_stats[1] + 3 * bp_stats[3] and bp_[4] >= bp_stats[1] - 3 * bp_stats[3]:
			if (bp_[2] == '+' and bp_[1] > bp[1]) or (bp_[2] == '-' and bp_[1] < bp[1]):
				bp[1] = bp_[1]
			if (bp_[5] == '+' and bp_[4] > bp[4]) or (bp_[5] == '-' and bp_[4] < bp[4]):
				bp[4] = bp_[4]
	#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp = %s" %(bp))
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
	if len(bpr) == 0:
		return bp, bpr, [0, 0, 0, 0, 0, 0]
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
	return bp, bpr, bp_stats_


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
	"""
	nodes: adjacent list - keys with format (chr, pos); 
	vals = [[flags], [sequence edges], [concordant edges], [discordant edges], [source edges]]  
	"""
	nodes = dict()
	special_nodes = dict() 
	"""
	discordant_edges_pos: For each concordant edge, which side has a discordant edge
	[[], []] = No discordant edges, i.e., CN boundary breakpoint
	[[], [edge indices]] = Only incoming breakpoint edge 
	[[edge indices], []] = Only outgoing breakpoint edge
	[[edge indices], [edge indices]] = Both incoming and outgoing breakpoint edges
	"""
	discordant_edges_pos = dict() 
	small_del_indices = [] # Indices of +- breakpoints on the same chr and < small_del_cutoff in discordant_edges
	source_edges = [] # Special AA discordant edges connected to a source node 
	new_bp_list_ = [] # Long read only breakpoints (discordant edges)
	new_bp_stats_ = [] # Statistics of long read only breakpoints (discordant edges)
	
	sr_length = 0.0 # Average short read length
	min_sr_alignment_length = 30.0 # Minimum alignment length of short reads - used for CN assignment
	max_sr_insert = 2000
	normal_cov_sr = 0 # Normal short read coverage - used for CN assignment
	normal_cov_lr = 0 # Normal long read coverage - used for CN assignment
	

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
			cns_intervals_median.append(cns_intervals[log2_cn_order[ip]])
			cns_intervals_median.append(cns_intervals[log2_cn_order[im]])
			log2_cn_median.append(log2_cn[log2_cn_order[ip]])
			log2_cn_median.append(log2_cn[log2_cn_order[im]])
			total_int_len += (cns_intervals[log2_cn_order[ip]][2] - cns_intervals[log2_cn_order[ip]][1] + 1)
			total_int_len += (cns_intervals[log2_cn_order[im]][2] - cns_intervals[log2_cn_order[im]][1] + 1)
			i = 1
			while total_int_len < 10000000:
				cns_intervals_median.append(cns_intervals[log2_cn_order[ip + i]])
				cns_intervals_median.append(cns_intervals[log2_cn_order[im - i]])
				log2_cn_median.append(log2_cn[log2_cn_order[ip]])
				log2_cn_median.append(log2_cn[log2_cn_order[im]])
				total_int_len += (cns_intervals[log2_cn_order[ip + i]][2] - cns_intervals[log2_cn_order[ip + i]][1] + 1)
				total_int_len += (cns_intervals[log2_cn_order[im - i]][2] - cns_intervals[log2_cn_order[im - i]][1] + 1)
				i += 1
			if cns == sr_cns:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Use %d SR copy number segments." %(len(cns_intervals_median)))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total length of SR copy number segments: %d." %(total_int_len))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Average SR copy number: %f." %(np.average(log2_cn_median)))
				nnc = 0
				for i in range(len(cns_intervals_median)):
					nnc += sum([sum(nc) for nc in self.sr_bamfh.count_coverage(cns_intervals_median[i][0], \
												cns_intervals_median[i][1], \
												cns_intervals_median[i][2] + 1)])
				self.normal_cov_sr = nnc * 1.0 / total_int_len
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
		while (chr, pos_) in self.nodes.keys():
			if pos_ != pos and (chr, pos_ - 1, pos_) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_ - 1, pos_)][1]) > 0:
				break
			if cr >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_)][1][0]][-1]
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
		while (chr, pos_ - 1) in self.nodes.keys():
			if pos_ != pos and (chr, pos_ - 1, pos_) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_ - 1, pos_)][1]) > 0:
				break
			if cl >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_ - 1)][1][0]][-1]
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
		while (chr, pos_ + 1) in self.nodes.keys():
			if pos_ != pos and (chr, pos_, pos_ + 1) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_, pos_ + 1)][0]) > 0:
				break
			if cr >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_ + 1)][1][0]][-1]
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
		while (chr, pos_) in self.nodes.keys():
			if pos_ != pos and (chr, pos_, pos_ + 1) in self.discordant_edges_pos and \
				len(self.discordant_edges_pos[(chr, pos_, pos_ + 1)][0]) > 0:
				break
			if cl >= self.min_bp_match_cutoff_:
				break
			seglen = self.seq_edges[self.nodes[(chr, pos_)][1][0]][-1]
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
					self.nodes[(chr, l)] = [[], [len(self.seq_edges)], [], [], []]
					self.nodes[(chr, r)] = [[], [len(self.seq_edges)], [], [], []]
					if l == r:
						self.special_nodes[(chr, l)] = 2
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
							self.nodes[(t0_[0], int(t0_[1][:-1]))][0].append(len(self.discordant_edges) - 1)
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge connected to interval ends.")
						self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1)][0].append(len(self.discordant_edges) - 1)	
					else:
						if (t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1])) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1]))] = [[], []]
							self.nodes[(t0_[0], int(t0_[1][:-1]))][0].append(len(self.discordant_edges) - 1)
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge connected to interval ends.")
						self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1]))][1].append(len(self.discordant_edges) - 1)
					if t1_[1][-1] == '+':
						if (t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1)] = [[], []]
							self.nodes[(t1_[0], int(t1_[1][:-1]))][0].append(len(self.discordant_edges) - 1)
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscordant edge connected to interval ends.")
						self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1)][0].append(len(self.discordant_edges) - 1)
					else:
						if (t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1])) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1]))] = [[], []]
							self.nodes[(t1_[0], int(t1_[1][:-1]))][0].append(len(self.discordant_edges) - 1)
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
						self.nodes[(t1_[0], int(t1_[1][:-1]))][4].append(len(self.source_edges))
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
				blocks = read.get_blocks()
				for bi in range(len(blocks) - 1):
					if abs(blocks[bi + 1][0] - blocks[bi][1]) > self.min_del_len:
						try:
							self.large_indel_alignments[rn].append([ai[0], blocks[bi + 1][0], blocks[bi][1]])
						except:
							self.large_indel_alignments[rn] = [[ai[0], blocks[bi + 1][0], blocks[bi][1]]]
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
			rcs = read.get_cigar_stats()[0]
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Fetched %d chimeric reads." %(len(self.chimeric_alignments)))
		for r in self.chimeric_alignments.keys():
			rl = self.read_length[r]
			r_int = []
			rr_int = []
			q_ = []
			for sa in self.chimeric_alignments[r]:
				t = sa.split(',')
				if 'S' not in t[3] or 'M' not in t[3]:
					"""
					Require a chimeric alignment record having at least some (soft)clips and matches 
					"""
					logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found chimeric alignment without match or soft clips.")
					logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tRead name: %s; Read length: %d." %(r, rl))
					logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAll CIGAR strings: %s." %(self.chimeric_alignments[r]))
					continue
				op = ''.join(c for c in t[3] if not c.isdigit())
				rs, re, rrl = chimeric_reads.cigar2pos_ops[op](t[3], t[2], rl)
				r_int.append([rs, re])
				if t[2] == '+':
					rr_int.append([t[0], int(t[1]) - 1, int(t[1]) + rrl - 2, '+']) # converted to 0 based coordinates
				else:
					rr_int.append([t[0], int(t[1]) + rrl - 2, int(t[1]) - 1, '-']) # converted to 0 based coordinates
				q_.append(int(t[4]))
			r_int_ind = sorted(range(len(r_int)), key = lambda i: (r_int[i][0], r_int[i][1]))
			r_int = [r_int[i] for i in r_int_ind]
			rr_int = [rr_int[i] for i in r_int_ind]
			q_ = [q_[i] for i in r_int_ind]
			self.chimeric_alignments[r] = [r_int, rr_int, q_]
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
							rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]]], d, 
							int(r_int[i + 1][0]) - int(r_int[i][1]), self.min_bp_match_cutoff[di]):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add(r)
						overlap_i = 1
						bassigned[i] = 1
					elif bp_match([rr_int[i + 1][0], rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]], 
							rr_int[i][0], rr_int[i][2], rr_int[i][3]], d, 
							int(r_int[i + 1][0]) - int(r_int[i][1]), self.min_bp_match_cutoff[di]):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add(r)
						overlap_i = 1
						bassigned[i] = 1
				"""
				Add unmatched breakpoint to new_bp_list
				"""
				if overlap_i == 0 and interval_overlap_l(rr_int[i], self.amplicon_intervals) >= 0 and \
					interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) >= 0:
					if rr_int[i + 1][0] != rr_int[i][0] or rr_int[i + 1][3] != rr_int[i][3]:
						if q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list.append(interval2bp(rr_int[i], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '+':
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list.append(interval2bp(rr_int[i], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '-':
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list.append(interval2bp(rr_int[i], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
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
								rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]]], d, 
								int(r_int[i + 1][0]) - int(r_int[i - 1][1]), self.min_bp_match_cutoff[di]):
							self.discordant_edges[di][7] += 1
							self.discordant_edges[di][-1].add(r)
							overlap_i = 1
						elif bp_match([rr_int[i + 1][0], rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]], 
								rr_int[i - 1][0], rr_int[i - 1][2], rr_int[i - 1][3]], d, 
								int(r_int[i + 1][0]) - int(r_int[i - 1][1]), self.min_bp_match_cutoff[di]):
							self.discordant_edges[di][7] += 1
							self.discordant_edges[di][-1].add(r)
							overlap_i = 1
				"""
				Add unmatched breakpoint to new_bp_list
				"""
				if overlap_i == 0 and bassigned[i - 1] == 0 and bassigned[i] == 0 and q_[i] < 10 and q_[i - 1] >= 20 and q_[i + 1] >= 20 and \
					interval_overlap_l(rr_int[i - 1], self.amplicon_intervals) >= 0 and interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) >= 0:
					if rr_int[i + 1][0] != rr_int[i - 1][0] or rr_int[i + 1][3] != rr_int[i - 1][3]:
						new_bp_list.append(interval2bp(rr_int[i - 1], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '+':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i - 1][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list.append(interval2bp(rr_int[i - 1], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '-':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i - 1][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list.append(interval2bp(rr_int[i - 1], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found %d reads with new breakpoints." %(len(new_bp_list)))
	
		new_bp_clusters = cluster_bp_list(new_bp_list, self.min_cluster_cutoff, self.max_bp_distance_cutoff)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "These reads formed %d clusters." %(len(new_bp_clusters)))
		for c in new_bp_clusters:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "New cluster of size %d." %(len(c)))
			if len(c) >= self.min_cluster_cutoff:
				bp, bpr, bp_stats_ = bpc2bp(c, self.min_bp_match_cutoff_)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp = %s" %(bp))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNum long read support = %d" %(len(set(bpr))))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp_stats = %s" %(bp_stats_))
				if len(set(bpr)) >= self.min_cluster_cutoff:
					if interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals) >= 0 and \
						interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals) >= 0:
						self.new_bp_list_.append(bp[:-3] + [set(bpr)])
						self.new_bp_stats_.append(bp_stats_)
				else:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the cluster.")
			else:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the cluster.")


	def find_smalldel_breakpoints(self):
		"""
		For each alignment record with large indels, first try to match the resulting breakpoint to the list of AA breakpoints
		Cluster the unmatched breakpoints
		"""
		new_bp_list = []
		for r in self.large_indel_alignments.keys():
			for rr_gap in self.large_indel_alignments[r]:
				overlap_ = 0
				rr_gap_ = rr_gap
				if rr_gap[2] > rr_gap[1]:
					rr_gap_[2] = rr_gap[1]
					rr_gap_[1] = rr_gap[2]
				for di in range(len(self.discordant_edges)):
					d = self.discordant_edges[di]
					if bp_match([rr_gap_[0], rr_gap_[1], '-', rr_gap_[0], rr_gap_[2], '+'], d, 0, \
							self.min_bp_match_cutoff[di]):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add(r)
						overlap_ = 1
				if overlap_ == 0:
					new_bp_list.append([rr_gap_[0], rr_gap_[1], '-', rr_gap_[0], rr_gap_[2], '+', r, 0, 0, -1, -1])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found %d reads with new small del breakpoints." %(len(new_bp_list)))

		new_bp_clusters = cluster_bp_list(new_bp_list, self.min_cluster_cutoff, self.max_bp_distance_cutoff)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "These reads formed %d clusters." %(len(new_bp_clusters)))
		for c in new_bp_clusters:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "New cluster of size %d." %(len(c)))
			if len(c) >= self.min_cluster_cutoff:
				bp, bpr, bp_stats_ = bpc2bp(c, self.min_bp_match_cutoff_)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp = %s" %(bp))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNum long read support = %d" %(len(set(bpr))))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tbp_stats = %s" %(bp_stats_))
				if len(set(bpr)) >= self.min_cluster_cutoff:
					if interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals) >= 0 and \
						interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals) >= 0:
						self.new_bp_list_.append(bp[:-3] + [set(bpr)])
						self.new_bp_stats_.append(bp_stats_)
				else:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the cluster.")
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
			sseg = self.seq_edges[segi]	
			new_segs[segi] = []
			lssi = 0
			del_list = []
			for ssi in range(len(split_seg[segi])):
				spp = split_seg[segi][ssi]
				if ssi > 0 and int(spp[0]) - int(split_seg[segi][lssi][1]) <= self.min_bp_match_cutoff_ and \
					spp[3] == split_seg[segi][lssi][3]:
					if self.new_bp_list_[spp[2]][spp[3] + 1] == '+':
						self.new_bp_list_[spp[2]][spp[3]] = split_seg[segi][lssi][0]
					else:
						self.new_bp_list_[spp[2]][spp[3]] = split_seg[segi][lssi][1]
					del_list.append(ssi)
				else:
					lssi = ssi
			for dssi in del_list[::-1]:
				del split_seg[segi][dssi]
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
					self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi][2])] = [[], [], [], [], []]
					self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi + 1][1])] = [[], [], [], [], []]
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNew node %s; adjacent list = %s." 
							%(str((new_segs[segi][nsi][0], new_segs[segi][nsi][2])), self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi][2])]))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNew node %s; adjacent list = %s." 
							%(str((new_segs[segi][nsi][0], new_segs[segi][nsi + 1][1])), self.nodes[(new_segs[segi][nsi][0], new_segs[segi][nsi + 1][1])]))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num sequence edges in breakpoint graph = %d." %(len(self.seq_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num concordant edges in breakpoint graph = %d." %(len(self.concordant_edges)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num discordant edges in breakpoint graph = %d." %(len(self.discordant_edges) + len(self.new_bp_list_)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num source edges in breakpoint graph = %d." %(len(self.source_edges)))
				

	def del_sr_bp(self):
		"""
		Delete short read breakpoints without long read coverage;
		Merge the corresponding sequence edges where no breakpoints left at the end after breakpoint edge deletion.
		"""
		del_list = [] # The list of short read ONLY breakpoints to be deleted
		source_del_list = [] # The list of source edges to be deleted
		bpi_map = dict()
		bpi_ = 0
		for bpi in range(len(self.discordant_edges)):
			if self.discordant_edges[bpi][7] == 0 and (self.discordant_edges[bpi][6][0] < 0.5 or self.discordant_edges[bpi][6][1] < 0.5 * self.normal_cov_sr):
				del_list.append(bpi)
			else:
				bpi_map[bpi] = bpi_
				bpi_ += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will delete the discordant edges at index %s." %(del_list))
		srci_map = dict()
		srci_ = 0
		for srci in range(len(self.source_edges)):
			if self.source_edges[srci][7] > 0.5 and (self.source_edges[srci][6] < 0.5 or self.source_edges[srci][7] < 0.5 * self.normal_cov_sr):
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
			node1 = (ce[0], ce[1])
			node2 = (ce[3], ce[4])
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
				if node1 in self.special_nodes:
					self.special_nodes[node1] -= 1
					if self.special_nodes[node1] == 0:
						del self.nodes[node1]
						del self.special_nodes[node1]
				else:
					del self.nodes[node1]
				if node2 in self.special_nodes:
					self.special_nodes[node2] -= 1
					if self.special_nodes[node2] == 0:
						del self.nodes[node2]
						del self.special_nodes[node2]
				else:
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
				rl_list = [read.infer_read_length() for read in self.sr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if read.infer_read_length()]
				len_rl_list = len(rl_list)
				seg[3] = len_rl_list
				if len_rl_list > max_seg:
					max_seg = len_rl_list
					avg_rl = np.average(rl_list) # Average SR length
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "SR cov assigned for sequence edge %s." %(seg))
			if seg[4] == -1:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Finding LR cov for sequence edge %s." %(seg))
				"""
				For long read, use the total number of nucleotides
				"""
				rl_list = [read.infer_read_length() for read in self.lr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if read.infer_read_length()]
				if lr_seq_dist == 'poisson':
					seg[4] = [len(rl_list), np.average(rl_list)]
				else:
					seg[4] = [len(rl_list), sum([sum(nc) for nc in self.lr_bamfh.count_coverage(seg[0], seg[1], seg[2], quality_threshold = 0, read_callback = 'nofilter')])]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "LR cov assigned for sequence edge %s." %(seg))
		if self.sr_length == 0.0:
			if avg_rl == 0.0:
				logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "No new sequence edges found, compute SR read length on the first sequence edge in seq_edges.")
				seg = self.seq_edges[0]
				rl_list = [read.infer_read_length() for read in self.sr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if read.infer_read_length()]
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
				rls_ = set([read.query_name for read in self.sr_bamfh.fetch(contig = ec[0], start = ec[1], stop = ec[1] + 1)
						if not read.is_unmapped and read.is_proper_pair and read.next_reference_name == ec[0] and
						abs(read.next_reference_start - read.reference_start) <= self.max_sr_insert])
				rrs_ = set([read.query_name for read in self.sr_bamfh.fetch(contig = ec[3], start = ec[4], stop = ec[4] + 1)
						if not read.is_unmapped and read.is_proper_pair and read.next_reference_name == ec[0] and
						abs(read.next_reference_start - read.reference_start) <= self.max_sr_insert])
				self.concordant_edges[eci][6] = len(rls_ & rrs_)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSR cov assigned for concordant edge %s." %(self.concordant_edges[eci]))
			rls = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[0], start = ec[1], stop = ec[1] + 1)])
			rrs = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[3], start = ec[4], stop = ec[4] + 1)])
			rls1 = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[0], start = ec[1] - self.min_bp_match_cutoff_ - 1, stop = ec[1] - self.min_bp_match_cutoff_)])
			rrs1 = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[3], start = ec[4] + self.min_bp_match_cutoff_, stop = ec[4] + self.min_bp_match_cutoff_ + 1)]) 
			rbps = set([])
			for bpi in self.discordant_edges_pos[(ec[0], ec[1], ec[4])][0]:
				if bpi >= len(self.discordant_edges):
					rbps |= self.new_bp_list_[bpi - len(self.discordant_edges)][-1]
				else:
					rbps |= self.discordant_edges[bpi][-1]
			for bpi in self.discordant_edges_pos[(ec[0], ec[1], ec[4])][1]:
				if bpi >= len(self.discordant_edges):
					rbps |= self.new_bp_list_[bpi - len(self.discordant_edges)][-1]
				else:
					rbps |= self.discordant_edges[bpi][-1]
			self.concordant_edges[eci][7] = len((rls & rrs & rls1 & rrs1) - rbps)
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tLR cov assigned for concordant edge %s." %(self.concordant_edges[eci]))


	def output_breakpoint_graph(self, ogfile):
		"""
		Write a breakpoint graph to file in AA graph format
		"""
		with open(ogfile, 'w') as fp:
			fp.write("SequenceEdge: StartPosition, EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads, Size\n")
			for segi in range(len(self.seq_edges)):
				sseg = self.seq_edges[segi]
				fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%s\t%d\t%s\n" %(sseg[0], sseg[1], sseg[0], sseg[2], sseg[-1], sseg[3], sseg[4][0], str(sseg[5])))
			fp.write("BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads, ")
			fp.write("HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence\n")
			for es in self.source_edges:
				fp.write("source\t%s:%s%s->%s:%s%s\t%f\t%s\t-1\t%s\t%s\n" %(es[0], es[1], es[2], es[3], 
					es[4], es[5], es[-1], str(es[7]), es[8], es[9]))
			for ec in self.concordant_edges:				
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


	def assign_cn(self, lr_seq_dist = 'poisson'):
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
		endnodes = [] # Nodes corresponding to interval ends
		for segi in range(lseg):
			sseg = self.seq_edges[segi]
			if len(self.nodes[(sseg[0], sseg[1])][1]) == 0:
				self.nodes[(sseg[0], sseg[1])][1].append(segi)
			else:
				self.nodes[(sseg[0], sseg[1])][1] = [segi]
			if len(self.nodes[(sseg[0], sseg[2])][1]) == 0:
				self.nodes[(sseg[0], sseg[2])][1].append(segi)
			else:
				self.nodes[(sseg[0], sseg[2])][1] = [segi]
			if segi == 0:
				endnodes.append((sseg[0], sseg[1]))
			else:
				lastseg = self.seq_edges[segi - 1]
				if lastseg[0] != sseg[0] or lastseg[2] + 1 != sseg[1]:
					endnodes.append((lastseg[0], lastseg[2]))
					endnodes.append((sseg[0], sseg[1]))
				if segi == lseg - 1:
					endnodes.append((sseg[0], sseg[2]))
		assert len(endnodes) == len(self.amplicon_intervals) * 2
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "The following nodes correspond to interval ends.")
		for node in endnodes:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNode %s." %(str(node)))
		for eci in range(lc):
			ec = self.concordant_edges[eci]
			self.nodes[(ec[0], ec[1])][2].append(eci)
			self.nodes[(ec[3], ec[4])][2].append(eci)
		for edi in range(ld):
			ed = self.discordant_edges[edi]
			self.nodes[(ed[0], ed[1])][3].append(edi)
			self.nodes[(ed[3], ed[4])][3].append(edi)
		for bpi in range(lnbp):
			bp = self.new_bp_list_[bpi]
			self.nodes[(bp[0], bp[1])][3].append(ld + bpi)
			self.nodes[(bp[3], bp[4])][3].append(ld + bpi)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Updated adjacent list for %d nodes:" %(len(self.nodes)))
		for node in self.nodes.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Node %s; adjacent list = %s." %(str(node), self.nodes[node]))
		
		nconstraints = len([node for node in self.nodes.keys() if node not in endnodes])
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
			wlncn = [sseg[3] + sseg[4][0] for sseg in self.seq_edges]
		else:
			wlncn = [sseg[3] - 0.5 for sseg in self.seq_edges]
		wlncn += [ec[6] + ec[7] for ec in self.concordant_edges]
		wlncn += [ed[6][1] + ed[7] for ed in self.discordant_edges]
		wlncn += [len(bp[6]) for bp in self.new_bp_list_]
		wlncn += [es[7] if es[7] >= 1 else 0.1 for es in self.source_edges]
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
			if node not in endnodes:
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
				self.seq_edges[segi] += [sol['x'][segi]]
			for eci in range(lc):
				self.concordant_edges[eci] += [sol['x'][lseg + eci]]
			for edi in range(ld):
				self.discordant_edges[edi] += [sol['x'][lseg + lc + edi]]
			for bpi in range(lnbp):
				self.new_bp_list_[bpi] += [sol['x'][lseg + lc + ld + bpi]]
			for esi in range(len(self.source_edges)):
				self.source_edges[esi] += [sol['x'][lseg + lc + ld + lnbp + esi]]
			


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("--sr_bam", help = "Sorted indexed short read bam file.", required = True)
	parser.add_argument("--lr_bam", help = "Sorted indexed long read bam file.", required = True)
	parser.add_argument("--aa_graph", help = "AA-formatted graph file.", required = True)
	parser.add_argument("--aa_cycle", help = "AA-formatted cycle file.", required = True)
	parser.add_argument("--output_bp", help = "If specified, only output the list of breakpoints.", action = 'store_true')
	parser.add_argument("--sr_cnseg", help = "Short read *.cns file.")
	parser.add_argument("--lr_cnseg", help = "Long read *.cns file.")
	parser.add_argument("--lr_seq_normal", help = "Use normal distribution on the total number of nucleotides on sequence edges for long read.", action = 'store_true')
	args = parser.parse_args()

	logging.basicConfig(filename = 'refine_breakpoint_graph.log', filemode = 'w', level = logging.DEBUG, 
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
	b2bn.del_sr_bp()
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
			b2bn.assign_cn(lr_seq_dist = 'normal')
		else:
			b2bn.assign_cn(lr_seq_dist = 'poisson')
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Assigned CN for all edges.")
		b2bn.output_breakpoint_graph(args.aa_graph.split('/')[-1][:-4] + '_.txt')
	b2bn.closebam()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total runtime.")


