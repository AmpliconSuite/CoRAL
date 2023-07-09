import time
start_time = time.time()
import logging
import math
import pysam
import argparse
import sys
import os
import numpy as np
from scipy import stats
from collections import Counter
import cvxopt
import cvxopt.modeling
import intervaltree

import cigar_parsing
from long_read_aa.breakpoints import global_names


def interval_overlap(int1, int2):
	"""
	Check if two intervals in the form of [chr, s, e] overlap
	"""
	return (int1[0] == int2[0] and int(int1[1]) <= int(int2[2]) and int(int2[1]) <= int(int1[2]))


def interval_adjacent(int1, int2):
	"""
	Check if two intervals in the form of [chr, s, e] are adjacent
	"""
	if int1[0] != int2[0]:
		return False
	if int1[1] <= int2[1]:
		return (int2[1] == int1[2] + 1)
	else:
		return (int1[1] == int2[2] + 1)


def interval_overlap_l(int1, intl):
	"""
	Check if an interval in the form of [chr, s, e] overlaps with a list of intervals
	"""
	for int2i in range(len(intl)):
		if interval_overlap(int1, intl[int2i]):
			return int2i
	return -1


def interval_exclusive(int1, intl):
	overlap_ints = set([])
	intl_ = [[intj for intj in int1]]
	for int2i in range(len(intl)):
		for inti_ in range(len(intl_))[::-1]:
			int_ = intl_[inti_]
			if interval_overlap(int_, intl[int2i]):
				overlap_ints.add(int2i)
				del intl_[inti_]
				if int_[1] < intl[int2i][1]:
					intl_.append([int_[0], int_[1], intl[int2i][1] - 1, -1])
				if int_[2] > intl[int2i][2]:
					intl_.append([int_[0], intl[int2i][2] + 1, int_[2], -1])
	return overlap_ints, intl_


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


def interval2bp(R1, R2, rn = '', rgap = 0):
	"""
	Convert split/chimeric alignment to breakpoint
	"""
	if (global_names.chr_idx[R2[0]] < global_names.chr_idx[R1[0]]) or (global_names.chr_idx[R2[0]] == global_names.chr_idx[R1[0]] and R2[1] < R1[2]):
		return [R1[0], R1[2], R1[3], R2[0], R2[1], global_names.neg_plus_minus[R2[3]], rn, rgap, 0]
	return [R2[0], R2[1], global_names.neg_plus_minus[R2[3]], R1[0], R1[2], R1[3], rn, rgap, 1]


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


class bam_to_breakpoint_nanopore():

	lr_bamfh = "" # Long read bam file

	max_interval_cutoff = 2000000
	min_bp_match_cutoff_ = 100
	interval_delta = 100000 # Not in use for now
	min_cluster_cutoff = 3 # Hard cutoff for considering a long read breakpoint cluster
	max_breakpoint_distance_cutoff = 2000 # Used for breakpoint clustering - if the distance of two breakpoint positions are greater than this cutoff, then start a new cluster 
	small_del_cutoff = 10000 # +- breakpoints (small deletions) with the two ends less than this cutoff are treated specially
	min_del_len = 600 # The minimum length of all +- (deletion) breakpoints returned by AA  

	read_length = dict() # Map read name -> read length
	chimeric_alignments = dict() # Map read name -> chimeric alignments (two or more records for one read)
	chimeric_alignments_seg = dict()
	large_indel_alignments = dict() # Map read name -> alignments with one record per read but large indels showing in CIGAR string
	
	amplicon_intervals = [] # AA amplicon intervals
	amplicon_interval_connections = dict()
	
	cns_intervals = [] # CN segments
	cns_intervals_by_chr = dict()
	log2_cn = [] # log2 CN for each CN segment
	cns_tree = dict() # Interval tree structure for each chromosome
	normal_cov = 0.0 # Normal long read coverage - used for CN assignment

	new_bp_list = [] # List of breakpoints (discordant edges)
	new_bp_stats = [] # Statistics of breakpoints (discordant edges)
	new_bp_ccids = []

	seq_edges = [] # Sequence edges
	seq_edge_ccids = []
	concordant_edges = [] # Concordant edges
	concordant_edge_ccids = []
	"""
	nodes: adjacent list - keys with format (chr, pos); 
	vals = [[flags], [sequence edges], [concordant edges], [discordant edges], [source edges]]  
	"""
	nodes = dict()
	endnodes = [] 
	small_del_indices = [] # Indices of +- breakpoints on the same chr and < small_del_cutoff in discordant_edges
	source_edges = [] # Special AA discordant edges connected to a source node
	source_edge_ccids = []
	max_CN = dict()



	def __init__(self, lr_bamfile, seedfile):
		self.lr_bamfh = pysam.AlignmentFile(lr_bamfile, 'rb')
		with open(seedfile, 'r') as fp:
			for line in fp:
				s = line.strip().split()
				self.amplicon_intervals.append([s[0], int(s[1]), int(s[2]), -1])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Parsed %d seed amplicon intervals." %(len(self.amplicon_intervals)))
		for ai in self.amplicon_intervals:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Seed interval %s" %(ai[:3]))


	def read_cns(self, cns):
		"""
		Read in (cnvkit) *.cns file and estimate the normal long read coverage
		"""
		self.cns_intervals = []
		self.log2_cn = []
		idx = 0		
		with open(cns, 'r') as fp:
			for line in fp:
				s = line.strip().split()
				if s[0] != "chromosome":
					self.cns_intervals.append([s[0], int(s[1]), int(s[2]) - 1]) 
					if s[0] not in self.cns_tree:
						self.cns_tree[s[0]] = intervaltree.IntervalTree()
						self.cns_intervals_by_chr[s[0]] = []
						idx = 0
					self.cns_intervals_by_chr[s[0]].append([s[0], int(s[1]), int(s[2]) - 1])
					self.cns_tree[s[0]][int(s[1]): int(s[2])] = idx
					idx += 1
					self.log2_cn.append(float(s[4]))
		#print (self.cns_tree)
		
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total num LR copy number segments: %d." %(len(self.log2_cn)))
		log2_cn_order = np.argsort(self.log2_cn)
		cns_intervals_median = []
		log2_cn_median = []
		im = int(len(log2_cn_order) / 2.4)
		ip = im + 1
		total_int_len = 0
		cns_intervals_median.append(self.cns_intervals[log2_cn_order[ip]])
		cns_intervals_median.append(self.cns_intervals[log2_cn_order[im]])
		log2_cn_median.append(self.log2_cn[log2_cn_order[ip]])
		log2_cn_median.append(self.log2_cn[log2_cn_order[im]])
		total_int_len += (self.cns_intervals[log2_cn_order[ip]][2] - self.cns_intervals[log2_cn_order[ip]][1] + 1)
		total_int_len += (self.cns_intervals[log2_cn_order[im]][2] - self.cns_intervals[log2_cn_order[im]][1] + 1)
		i = 1
		while total_int_len < 10000000:
			cns_intervals_median.append(self.cns_intervals[log2_cn_order[ip + i]])
			cns_intervals_median.append(self.cns_intervals[log2_cn_order[im - i]])
			log2_cn_median.append(self.log2_cn[log2_cn_order[ip]])
			log2_cn_median.append(self.log2_cn[log2_cn_order[im]])
			total_int_len += (self.cns_intervals[log2_cn_order[ip + i]][2] - self.cns_intervals[log2_cn_order[ip + i]][1] + 1)
			total_int_len += (self.cns_intervals[log2_cn_order[im - i]][2] - self.cns_intervals[log2_cn_order[im - i]][1] + 1)
			i += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Use %d LR copy number segments." %(len(cns_intervals_median)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total length of LR copy number segments: %d." %(total_int_len))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Average LR copy number: %f." %(np.average(log2_cn_median)))
		nnc = 0
		for i in range(len(cns_intervals_median)):
			nnc += sum([sum(nc) for nc in self.lr_bamfh.count_coverage(cns_intervals_median[i][0], cns_intervals_median[i][1], \
						cns_intervals_median[i][2] + 1, quality_threshold = 0, read_callback = 'nofilter')])
		self.normal_cov = nnc * 1.0 / total_int_len
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "LR normal cov = %f." %(self.normal_cov))
		self.min_cluster_cutoff = max(self.min_cluster_cutoff, 0.5 * self.normal_cov)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset min_cluster_cutoff to %f." %(self.min_cluster_cutoff))


	def fetch(self):
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
			self.chimeric_alignments[r] = cigar_parsing.alignment_from_satags(self.chimeric_alignments[r], rl)
		for r in reads_wo_primary_alignment:
			del self.chimeric_alignments[r]
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Computed alignment intervals on all chimeric reads.")


	def pos2cni(self, chr, pos):
		return self.cns_tree[chr][pos]


	def hash_alignment_to_seg(self):
		"""
		"""
		# TBD: improve the running time by a better algorithmic strategy
		for read in self.chimeric_alignments.keys():
			for ri in range(len(self.chimeric_alignments[read][1])):
				rint = self.chimeric_alignments[read][1][ri]
				if rint[0] in self.cns_tree:
					lcni_ = self.pos2cni(rint[0], min(rint[1], rint[2]))
					rcni_ = self.pos2cni(rint[0], max(rint[1], rint[2]))
					assert len(lcni_) <= 1 and len(rcni_) <= 1
					lcni, rcni = -1, -1
					if len(lcni_) == 1:
						lcni = list(lcni_)[0].data
					if len(rcni_) == 1:
						rcni = list(rcni_)[0].data
					cniset = set([lcni, rcni])
					if len(cniset) > 1 and -1 in cniset:
						cniset.remove(-1)
					self.chimeric_alignments[read][1][ri].append(cniset)
					if rint[0] not in self.chimeric_alignments_seg:
						self.chimeric_alignments_seg[rint[0]] = dict()
					for cni in cniset:
						if cni != -1:
							try:
								self.chimeric_alignments_seg[rint[0]][cni].append(read)
							except:
								self.chimeric_alignments_seg[rint[0]][cni] = [read]
				else:
					self.chimeric_alignments[read][1][ri].append(set([-1]))
	

	def find_amplicon_intervals(self):
		# Reset seed intervals
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Updating seed amplicon intervals according to CN segments.")
		for ai in range(len(self.amplicon_intervals)):
			chr = self.amplicon_intervals[ai][0]
			lcni = list(self.pos2cni(chr, self.amplicon_intervals[ai][1]))[0].data
			rcni = list(self.pos2cni(chr, self.amplicon_intervals[ai][2]))[0].data
			self.amplicon_intervals[ai][1] = self.cns_intervals_by_chr[chr][lcni][1]
			self.amplicon_intervals[ai][2] = self.cns_intervals_by_chr[chr][rcni][2]
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tUpdated amplicon interval %s" %self.amplicon_intervals[ai])
 
		ccid = 0
		for ai in range(len(self.amplicon_intervals)):
			if self.amplicon_intervals[ai][3] == -1:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Begin processing amplicon interval %d" %ai)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAmplicon interval %s" %self.amplicon_intervals[ai])
				self.find_interval_i(ai, ccid)
				ccid += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Identified %d amplicon intervals in total." %len(self.amplicon_intervals))
		sorted_ai_indices = sorted(range(len(self.amplicon_intervals)), key = lambda i: (global_names.chr_idx[self.amplicon_intervals[i][0]], self.amplicon_intervals[i][1]))
		amplicon_intervals_sorted = [self.amplicon_intervals[i] for i in sorted_ai_indices]
		for ai in range(len(amplicon_intervals_sorted)):
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAmplicon interval %s" %amplicon_intervals_sorted[ai])
		# Merge amplicon intervals
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Begin merging adjacent intervals.")
		lastai = 0
		intervals_to_merge = []
		for ai in range(len(amplicon_intervals_sorted) - 1):
			if not interval_adjacent(amplicon_intervals_sorted[ai + 1], amplicon_intervals_sorted[ai]):
				if ai > lastai: 
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Merging intervals from %d to %d." %(lastai, ai))
					intervals_to_merge.append([lastai, ai])
				lastai = ai + 1
		if len(amplicon_intervals_sorted) > 0 and lastai < len(amplicon_intervals_sorted) - 1:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Merging intervals from %d to %d." %(lastai, len(amplicon_intervals_sorted) - 1))
			intervals_to_merge.append([lastai, len(amplicon_intervals_sorted) - 1])
		for int_range in intervals_to_merge[::-1]:
			# Reset interval 
			amplicon_intervals_sorted[int_range[0]][2] = amplicon_intervals_sorted[int_range[1]][2]
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset amplicon interval %d to %s." %(int_range[0], amplicon_intervals_sorted[int_range[0]]))
			# Modify ccid 
			for ai in range(int_range[0] + 1, int_range[1] + 1):
				if amplicon_intervals_sorted[ai][3] != amplicon_intervals_sorted[int_range[0]][3]:
					ccid = amplicon_intervals_sorted[ai][3]
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset amplicon intervals with ccid %d." %ccid)
					for ai_ in range(len(amplicon_intervals_sorted)):
						if amplicon_intervals_sorted[ai_][3] == ccid:
							amplicon_intervals_sorted[ai_][3] = amplicon_intervals_sorted[int_range[0]][3]
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset ccid of amplicon interval %d to %d." %(ai_, amplicon_intervals_sorted[int_range[0]][3]))
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Updated amplicon interval: %s" %amplicon_intervals_sorted[ai_])
			# Modify interval connections
			connection_map = dict()
			for connection in self.amplicon_interval_connections.keys():
				connection_map[connection] = connection
			#print (connection_map)
			for ai in range(int_range[0] + 1, int_range[1] + 1):
				ai_unsorted_ = sorted_ai_indices[int_range[0]]
				ai_unsorted = sorted_ai_indices[ai]
				
				for connection in connection_map.keys():
					#print ('a', ai_unsorted, ai_unsorted_, connection, connection_map[connection])
					if ai_unsorted == connection_map[connection][0]:
						connection_map[connection] = (ai_unsorted_, connection_map[connection][1])
					if ai_unsorted == connection_map[connection][1]:
						connection_map[connection] = (connection_map[connection][0], ai_unsorted_)
					if connection_map[connection][1] < connection_map[connection][0]:
						connection_map[connection] = (connection_map[connection][1], connection_map[connection][0])
					#print ('a-', connection_map[connection])
			#print (connection_map)
			for connection in connection_map.keys():
				if connection != connection_map[connection]:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset connection between amplicon intervals %s to %s." %(connection, connection_map[connection]))
					if connection_map[connection] not in self.amplicon_interval_connections:
						self.amplicon_interval_connections[connection_map[connection]] = self.amplicon_interval_connections[connection]
					else:
						self.amplicon_interval_connections[connection_map[connection]] |= self.amplicon_interval_connections[connection]
					del self.amplicon_interval_connections[connection]
					if connection_map[connection][0] == connection_map[connection][1]:
						del self.amplicon_interval_connections[connection_map[connection]]
						
			# Delete intervals
			for ai in range(int_range[0] + 1, int_range[1] + 1)[::-1]:
				ai_unsorted = sorted_ai_indices[ai]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Delete amplicon intervals %d - %s." %(ai_unsorted, self.amplicon_intervals[ai_unsorted]))
				del amplicon_intervals_sorted[ai]
				del sorted_ai_indices[ai]
		
		self.amplicon_intervals = [amplicon_intervals_sorted[ai] for ai in range(len(amplicon_intervals_sorted))]
		ind_map = {sorted_ai_indices[i]: i for i in range(len(sorted_ai_indices))}
		#print (ind_map, self.amplicon_interval_connections)
		connection_map = {connection: (min(ind_map[connection[0]], ind_map[connection[1]]), max(ind_map[connection[0]], ind_map[connection[1]])) for connection in self.amplicon_interval_connections.keys()}
		self.amplicon_interval_connections = {connection_map[connection]: self.amplicon_interval_connections[connection] for connection in self.amplicon_interval_connections.keys()}
		#print (self.amplicon_interval_connections)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d amplicon intervals after merging." %len(self.amplicon_intervals))
		for ai in range(len(self.amplicon_intervals)):
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAmplicon interval %s after merging." %self.amplicon_intervals[ai])
				 
		#print (len(self.amplicon_intervals))
		#print (self.amplicon_intervals)
		#print (self.amplicon_interval_connections.keys())

	
	def addbp(self, bp_, bpr_, bp_stats_, ccid):
		for bpi in range(len(self.new_bp_list)):
			bp = self.new_bp_list[bpi]
			if bp[0] == bp_[0] and bp[3] == bp_[3] and bp[2] == bp_[2] and bp[5] == bp_[5] and abs(bp[1] - bp_[1]) < 200 and abs(bp[4] - bp_[4]) < 200:
				self.new_bp_list[bpi][-1] |= set(bpr_)
				return bpi
		bpi = len(self.new_bp_list)
		self.new_bp_list.append(bp_ + [bpr_])
		self.new_bp_ccids.append(ccid)
		self.new_bp_stats.append(bp_stats_)
		return bpi


	def find_interval_i(self, ai, ccid):
		"""
		Given an amplification interval I indexed by ai, search for amplification intervals connected with I iteratively (with BFS)
			by a breakpoint edge
		Assign I a connected component id ccid if not already assigned
		"""
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tStart BFS on amplicon interval %d." %ai)
		L = [ai] # BFS queue
		while len(L) > 0:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tBFS queue: %s" %L)
			ai_ = L.pop(0)
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tNext amplicon interval %d: %s." %(ai_, self.amplicon_intervals[ai_]))
			chr = self.amplicon_intervals[ai_][0]
			s = self.amplicon_intervals[ai_][1]
			e = self.amplicon_intervals[ai_][2]
			if self.amplicon_intervals[ai_][3] == -1:
				self.amplicon_intervals[ai_][3] = ccid
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tReset connected component ID to %d" %ccid)

			# Identify all amplification intervals connected to interval indexed by ai_ with a breakpoint edge
			si = list(self.pos2cni(chr, s))[0].data
			ei = list(self.pos2cni(chr, e))[0].data
			
			d1_segs = dict() # All CN segments which share a chimeric alignment to the given interval 
			for i in range(si, ei + 1):
				if i in self.chimeric_alignments_seg[chr]:
					for r in self.chimeric_alignments_seg[chr][i]:
						rint = self.chimeric_alignments[r][1]
						for int_ in rint:
							for i_ in int_[-1]:
								if (int_[0] != chr) or (i_ < si or i_ > ei):
									try:
										if i_!= -1:
											d1_segs[int_[0]][i_].add(r)
									except:
										if int_[0] not in d1_segs:
											d1_segs[int_[0]] = dict()
										if i_ != -1:
											d1_segs[int_[0]][i_] = set([r])
			# Initial filtering of potential breakpoints of insufficient support
			chr_del_list = []
			for chr_ in d1_segs.keys():
				segs = d1_segs[chr_].keys()
				del_list = []
				for segi in segs:
					if len(d1_segs[chr_][segi]) < self.min_cluster_cutoff:
						del_list.append(segi)
				for segi in del_list:
					del d1_segs[chr_][segi]
				if len(d1_segs[chr_]) == 0:
					chr_del_list.append(chr_)
			for chr_ in chr_del_list:
				del d1_segs[chr_]

			new_intervals_refined = []
			new_intervals_connections = []
			for chr_ in d1_segs.keys():
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFound new intervals on chr %s" %chr_)
				new_intervals = [] # Initial list of new amplicon intervals
				sorted_bins_chr_ = sorted(d1_segs[chr_].keys())
				nir = set([])
				lasti = 0
				for i in range(0, len(sorted_bins_chr_) - 1):
					nil = self.cns_intervals_by_chr[chr_][sorted_bins_chr_[i + 1]][1]
					lir = self.cns_intervals_by_chr[chr_][sorted_bins_chr_[i]][2]
					if sorted_bins_chr_[i + 1] - sorted_bins_chr_[i] > 2 or nil - lir > self.max_interval_cutoff:
						nir |= d1_segs[chr_][sorted_bins_chr_[i]]
						new_intervals.append([chr_, sorted_bins_chr_[lasti], sorted_bins_chr_[i], nir])
						lasti = i + 1
						nir = set([])
					else:
						nir |= d1_segs[chr_][sorted_bins_chr_[i]]
				nir |= d1_segs[chr_][sorted_bins_chr_[-1]]
				new_intervals.append([chr_, sorted_bins_chr_[lasti], sorted_bins_chr_[-1], nir])
				
				# Refine initial intervals
				for nint_ in new_intervals:
					ns = self.cns_intervals_by_chr[nint_[0]][nint_[1]][1]
					ne = self.cns_intervals_by_chr[nint_[0]][nint_[2]][2]
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tRefining new interval %s." %[chr_, ns, ne])
					new_bp_list = []
					for r in nint_[-1]:
						r_int = self.chimeric_alignments[r][0]
						rr_int = self.chimeric_alignments[r][1]
						q_ = self.chimeric_alignments[r][2]
						bassigned = [0 for i in range(len(rr_int) - 1)]

						# Breakpoint from local alignment i and i + 1
						for ri in range(len(rr_int) - 1):
							if interval_overlap(rr_int[ri], [nint_[0], ns, ne]) and \
								interval_overlap(rr_int[ri + 1], self.amplicon_intervals[ai_]) and q_[ri] >= 20 and q_[ri + 1] >= 20:
								new_bp_list.append(interval2bp(rr_int[ri], rr_int[ri + 1], r, int(r_int[ri + 1][0]) - int(r_int[ri][1])) + [q_[ri], q_[ri + 1]])
								bassigned[ri] = 1
								#bp_flag_r = 1
							elif interval_overlap(rr_int[ri + 1], [nint_[0], ns, ne]) and \
								interval_overlap(rr_int[ri], self.amplicon_intervals[ai_]) and q_[ri] >= 20 and q_[ri + 1] >= 20:
								new_bp_list.append(interval2bp(rr_int[ri], rr_int[ri + 1], r, int(r_int[ri + 1][0]) - int(r_int[ri][1])) + [q_[ri], q_[ri + 1]])
								bassigned[ri] = 1

						# Breakpoint from local alignment i - 1 and i + 1
						for ri in range(1, len(rr_int) - 1):
							if bassigned[ri - 1] == 0 and bassigned[ri] == 0 and q_[ri] < 10 and q_[ri - 1] >= 20 and q_[ri + 1] >= 20 and \
								interval_overlap(rr_int[ri - 1], [nint_[0], ns, ne]) and \
								interval_overlap(rr_int[ri + 1], self.amplicon_intervals[ai_]):
								new_bp_list.append(interval2bp(rr_int[ri - 1], rr_int[ri + 1], r, int(r_int[ri + 1][0]) - int(r_int[ri - 1][1])) + [q_[ri - 1], q_[ri + 1]])
							elif bassigned[ri - 1] == 0 and bassigned[ri] == 0 and q_[ri] < 10 and q_[ri - 1] >= 20 and q_[ri + 1] >= 20 and \
								interval_overlap(rr_int[ri + 1], [nint_[0], ns, ne]) and \
								interval_overlap(rr_int[ri - 1], self.amplicon_intervals[ai_]):
								new_bp_list.append(interval2bp(rr_int[ri - 1], rr_int[ri + 1], r, int(r_int[ri + 1][0]) - int(r_int[ri - 1][1])) + [q_[ri - 1], q_[ri + 1]])

					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFound %d reads connecting the two intervals." %len(new_bp_list))
					new_bp_clusters = cluster_bp_list(new_bp_list, self.min_cluster_cutoff, self.max_breakpoint_distance_cutoff)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tThese reads formed %d clusters." %(len(new_bp_clusters)))
					new_bp_refined = []
					for c in new_bp_clusters:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\tNew cluster of size %d." %(len(c)))
						if len(c) >= self.min_cluster_cutoff:
							num_subcluster = 0
							bp_cluster_r = c
							while len(bp_cluster_r) >= self.min_cluster_cutoff:
								bp, bpr, bp_stats_, bp_cluster_r = bpc2bp(bp_cluster_r, self.min_bp_match_cutoff_)
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSubcluster %d" %(num_subcluster))
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\t\tbp = %s" %(bp))
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\t\tNum long read support = %d" %(len(set(bpr))))
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\t\tbp_stats = %s" %(bp_stats_))
								if (num_subcluster == 0 and len(set(bpr)) >= self.min_cluster_cutoff) or \
									(len(set(bpr)) >= max(self.normal_cov, 3.0)):
									bpi = self.addbp(bp, set(bpr), bp_stats_, ccid)
									if bpi not in new_bp_refined:
										new_bp_refined.append(bpi)
									logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\tKeeped the cluster.")
								else:
									logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\tDiscarded the cluster.")
						else:
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\tDiscarded the cluster.")
					
					nint_segs = []
					nint_segs_ = []
					if len(new_bp_refined) > 0:
						for bpi in new_bp_refined:
							bp = self.new_bp_list[bpi][:6]
							if interval_overlap([bp[0], bp[1], bp[1]], self.amplicon_intervals[ai_]) and interval_overlap([bp[3], bp[4], bp[4]], [nint_[0], ns, ne]):
								nint_segs.append([list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi])
							elif interval_overlap([bp[3], bp[4], bp[4]], self.amplicon_intervals[ai_]) and interval_overlap([bp[0], bp[1], bp[1]], [nint_[0], ns, ne]):
								nint_segs.append([list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi])
							else:
								logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tExact breakpoint outside amplicon interval.")
								logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tBreakpoint %s." %bp)
								logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tCurrent interval %s." %self.amplicon_intervals[ai_])
								logging.warning("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tNew interval %s." %[nint_[0], ns, ne])
								o1 = interval_overlap([bp[0], bp[1], bp[1]], [nint_[0], ns, ne])
								o2 = interval_overlap([bp[3], bp[4], bp[4]], [nint_[0], ns, ne])
								if o1 and o2:
									nint_segs.append([list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi])
									nint_segs.append([list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi])
								elif o1:
									nint_segs.append([list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi])
									nint_segs_.append([bp[3], list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi])
								elif o2:
									nint_segs_.append([bp[0], list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi])
									nint_segs.append([list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi])
								else:
									nint_segs_.append([bp[0], list(self.pos2cni(bp[0], bp[1]))[0].data, bp[1], bpi])
									nint_segs_.append([bp[3], list(self.pos2cni(bp[3], bp[4]))[0].data, bp[4], bpi])
						nint_segs = sorted(nint_segs, key = lambda item: (item[0], item[1]))
						nint_segs_ = sorted(nint_segs_, key = lambda item: (global_names.chr_idx[item[0]], item[1], item[2]))
						lasti = 0
						for i in range(0, len(nint_segs) - 1):
							nil = self.cns_intervals_by_chr[chr_][nint_segs[i + 1][0]][1]
							lir = self.cns_intervals_by_chr[chr_][nint_segs[i][0]][2]
							if nint_segs[i + 1][0] - nint_segs[i][0] > 2 or nil - lir > self.max_interval_cutoff or \
								nint_segs[i + 1][1] - nint_segs[i][1] > self.max_interval_cutoff:
								new_intervals_refined.append([chr_, max(self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1], nint_segs[lasti][1] - int(self.max_interval_cutoff / 2)),
											min(lir, nint_segs[i][1] + int(self.max_interval_cutoff / 2)), -1])
								new_intervals_connections.append([])
								for i_ in range(lasti, i + 1):
									new_intervals_connections[-1].append(nint_segs[i_][2])
								lasti = i + 1
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFixed new interval: %s." %new_intervals_refined[-1])
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tList of breakpoints connected to the new interval:")
								for bpi in new_intervals_connections[-1]:
									logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\t%s" %self.new_bp_list[bpi][:6])
						if len(nint_segs) > 0:
							new_intervals_refined.append([chr_, max(self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1], nint_segs[lasti][1] - int(self.max_interval_cutoff / 2)), 
									min(self.cns_intervals_by_chr[chr_][nint_segs[-1][0]][2], nint_segs[-1][1] + int(self.max_interval_cutoff / 2)), -1])
							new_intervals_connections.append([])
							for i_ in range(lasti, len(nint_segs)):
								new_intervals_connections[-1].append(nint_segs[i_][2])
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFixed new interval: %s." %new_intervals_refined[-1])
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tList of breakpoints connected to the new interval:")
							for bpi in new_intervals_connections[-1]:
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\t%s" %self.new_bp_list[bpi][:6])
						lasti = 0
						for i in range(0, len(nint_segs_) - 1):
							nil = self.cns_intervals_by_chr[nint_segs_[i + 1][0]][nint_segs_[i + 1][1]][1]
							lir = self.cns_intervals_by_chr[nint_segs_[i][0]][nint_segs_[i][1]][2]
							if (nint_segs_[i + 1][0] != nint_segs_[i][0]) or (nint_segs_[i + 1][0] == nint_segs_[i][0] and \
								(nint_segs_[i + 1][1] - nint_segs_[i][1] > 2 or nil - lir > self.max_interval_cutoff)) or \
								(nint_segs_[i + 1][0] == nint_segs_[i][0] and nint_segs_[i + 1][2] - nint_segs_[i][2] > self.max_interval_cutoff):
								new_intervals_refined.append([nint_segs_[lasti][0], max(self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1], nint_segs_[lasti][2] - int(self.max_interval_cutoff / 2)),
											min(lir, nint_segs_[i][2] + int(self.max_interval_cutoff / 2)), -1])
								new_intervals_connections.append([])
								lasti = i + 1
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFixed new interval: %s." %new_intervals_refined[-1])
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tSkip breakpoints connected to the new interval.")
						if len(nint_segs_) > 0:
							new_intervals_refined.append([nint_segs_[lasti][0], max(self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1], nint_segs_[lasti][2] - int(self.max_interval_cutoff / 2)), 
									min(self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[-1][1]][2], nint_segs_[-1][2] + int(self.max_interval_cutoff / 2)), -1])
							new_intervals_connections.append([])
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFixed new interval: %s." %new_intervals_refined[-1])
							logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tSkip breakpoints connected to the new interval:")

			# BFS
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tProcessing new intervals.")
			#print self.amplicon_intervals[ai_], new_intervals_refined
			for ni in range(len(new_intervals_refined)):
				ei, intl = interval_exclusive(new_intervals_refined[ni], self.amplicon_intervals)
				#print (new_intervals_refined[ni], ei, intl)
				if len(intl) == 0:
					ei_str = ""
					for ei_ in ei:
						ei_str += "%s " %self.amplicon_intervals[ei_]
					ei_str = ei_str.rstrip()	
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tNew interval %s overlaps with existing interval %s." 
									%(new_intervals_refined[ni], ei_str))
					for bpi in new_intervals_connections[ni]:
						bp = self.new_bp_list[bpi][:6]
						for ei_ in ei:
							connection = (min(ai_, ei_), max(ai_, ei_))
							if ei_!= ai_ and interval_overlap([bp[0], bp[1], bp[1]], self.amplicon_intervals[ei_]) or \
								interval_overlap([bp[3], bp[4], bp[4]], self.amplicon_intervals[ei_]):
								try:
									self.amplicon_interval_connections[connection].add(bpi)
								except:
									self.amplicon_interval_connections[connection] = set([bpi])
					for ei_ in ei:
						if ei_ != ai_ and self.amplicon_intervals[ei_][3] < 0:
							L.append(ei_)
				else:
					for int_ in intl:
						nai = len(self.amplicon_intervals)
						self.amplicon_intervals.append(int_)
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tAdded new interval %s to the amplicon interval list." 
								%int_)
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tNew interval index: %d." %nai) 
						self.amplicon_interval_connections[(ai_, nai)] = set([])
						if len(ei) == 0:
							#lbp_pos = 0
							#rbp_pos = 0
							for bpi in new_intervals_connections[ni]:
								"""
								bp = self.new_bp_list[bpi][:6]
								if interval_overlap([bp[0], bp[1], bp[1]], self.amplicon_intervals[nai]):
									if lbp_pos == 0:
										lbp_pos = bp[1]
									elif bp[1] < lbp_pos:
										lbp_pos = bp[1]
									if rbp_pos == 0:
										rbp_pos = bp[1]
									elif bp[1] > rbp_pos:
										rbp_pos = bp[1]
								if interval_overlap([bp[3], bp[4], bp[4]], self.amplicon_intervals[nai]):
									if lbp_pos == 0:
										lbp_pos = bp[4]
									elif bp[4] < lbp_pos:
										lbp_pos = bp[4]
									if rbp_pos == 0:
										rbp_pos = bp[4]
									elif bp[4] > rbp_pos:
										rbp_pos = bp[4]
								"""
								self.amplicon_interval_connections[(ai_, nai)].add(bpi)
							#if rbp_pos > 0 and self.amplicon_intervals[nai][2] > rbp_pos + self.max_interval_cutoff / 2:
							#	self.amplicon_intervals[nai][2] = rbp_pos + self.max_interval_cutoff / 2
							#if lbp_pos > 0 and self.amplicon_intervals[nai][1] < lbp_pos - self.max_interval_cutoff / 2:
							#	self.amplicon_intervals[nai][1] = lbp_pos - self.max_interval_cutoff / 2
						else:
							#lbp_pos = 0
							#rbp_pos = 0
							for bpi in new_intervals_connections[ni]:
								bp = self.new_bp_list[bpi][:6]
								for ei_ in ei:
									connection = (min(ai_, ei_), max(ai_, ei_))
									if interval_overlap([bp[0], bp[1], bp[1]], self.amplicon_intervals[ei_]) or \
										interval_overlap([bp[3], bp[4], bp[4]], self.amplicon_intervals[ei_]):
										try:
											self.amplicon_interval_connections[connection].add(bpi)
										except:
											self.amplicon_interval_connections[connection] = set([bpi])
									else:
										"""
										if interval_overlap([bp[0], bp[1], bp[1]], self.amplicon_intervals[nai]):
											if lbp_pos == 0:
												lbp_pos = bp[1]
											elif bp[1] < lbp_pos:
												lbp_pos = bp[1]
											if rbp_pos == 0:
												rbp_pos = bp[1]
											elif bp[1] > rbp_pos:
												rbp_pos = bp[1]
										if interval_overlap([bp[3], bp[4], bp[4]], self.amplicon_intervals[nai]):
											if lbp_pos == 0:
												lbp_pos = bp[4]
											elif bp[4] < lbp_pos:
												lbp_pos = bp[4]
											if rbp_pos == 0:
												rbp_pos = bp[4]
											elif bp[4] > rbp_pos:
												rbp_pos = bp[4]
										"""
										self.amplicon_interval_connections[(ai_, nai)].add(bpi)
							#if rbp_pos > 0 and self.amplicon_intervals[nai][2] > rbp_pos + self.max_interval_cutoff / 2:
							#	self.amplicon_intervals[nai][2] = rbp_pos + self.max_interval_cutoff / 2
							#if lbp_pos > 0 and self.amplicon_intervals[nai][1] < lbp_pos - self.max_interval_cutoff / 2:
							#	self.amplicon_intervals[nai][1] = lbp_pos - self.max_interval_cutoff / 2
						#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tAdded new interval %s to the amplicon interval list." 
						#		%int_)
						L.append(nai)
				

	def find_breakpoints(self):
		"""
		For each chimeric alignment, first try to match the resulting breakpoint to the list of AA breakpoints
		Cluster the unmatched breakpoints
		"""
		new_bp_list_ = []
		for r in self.chimeric_alignments.keys():
			r_int = self.chimeric_alignments[r][0]
			rr_int = self.chimeric_alignments[r][1]
			q_ = self.chimeric_alignments[r][2]
			bassigned = [0 for i in range(len(rr_int) - 1)]
			"""
			Breakpoint from local alignment i and i + 1
			"""
			for i in range(len(rr_int) - 1):
				"""
				Add unmatched breakpoint to new_bp_list
				"""
				io1 = interval_overlap_l(rr_int[i], self.amplicon_intervals)
				io2 = interval_overlap_l(rr_int[i + 1], self.amplicon_intervals)
				if io1 >= 0 and io2 >= 0 and io1 == io2:
					if rr_int[i + 1][3] != rr_int[i][3]:
						if q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list_.append(interval2bp(rr_int[i], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '+':
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list_.append(interval2bp(rr_int[i], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '-':
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list_.append(interval2bp(rr_int[i], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
			"""
			Breakpoint from local alignment i - 1 and i + 1
			"""		
			for i in range(1, len(rr_int) - 1):
				"""
				Add unmatched breakpoint to new_bp_list
				"""
				io1 = interval_overlap_l(rr_int[i - 1], self.amplicon_intervals)
				io2 = interval_overlap_l(rr_int[i + 1], self.amplicon_intervals)
				if bassigned[i - 1] == 0 and bassigned[i] == 0 and q_[i] < 10 and q_[i - 1] >= 20 and q_[i + 1] >= 20 and \
					io1 >= 0 and io2 >= 0 and io1 == io2:
					if rr_int[i + 1][3] != rr_int[i - 1][3]:
						new_bp_list_.append(interval2bp(rr_int[i - 1], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '+':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i - 1][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list_.append(interval2bp(rr_int[i - 1], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '-':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i - 1][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list_.append(interval2bp(rr_int[i - 1], rr_int[i + 1], r, int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found %d reads with new breakpoints." %(len(new_bp_list_)))
	
		new_bp_clusters = cluster_bp_list(new_bp_list_, self.min_cluster_cutoff, self.max_breakpoint_distance_cutoff)
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
					if (num_subcluster == 0 and len(set(bpr)) >= self.min_cluster_cutoff) or (len(set(bpr)) >= max(self.normal_cov, 3.0)):
						io1 = interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals)
						io2 = interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals)
						if io1 >= 0 and io2 >= 0:
							assert self.amplicon_intervals[io1][3] == self.amplicon_intervals[io2][3] 
							bpi = self.addbp(bp, set(bpr), bp_stats_, self.amplicon_intervals[io1][3])
							try:
								self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))].add(bpi)
							except:
								self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))] = set([bpi])
					else:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the subcluster %d." %(num_subcluster))
					num_subcluster += 1
			else:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the cluster.")



	def find_smalldel_breakpoints(self):
		"""
		For each alignment record with large indels, first try to match the resulting breakpoint to the list of AA breakpoints
		Cluster the unmatched breakpoints
		"""
		new_bp_list_ = []
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
		
		for r in self.large_indel_alignments.keys():
			for rr_gap in self.large_indel_alignments[r]:
				rr_gap_ = rr_gap
				if rr_gap[2] > rr_gap[1]:
					rr_gap_[2] = rr_gap[1]
					rr_gap_[1] = rr_gap[2]
				new_bp_list_.append([rr_gap_[0], rr_gap_[1], '-', rr_gap_[0], rr_gap_[2], '+', r, 0, 0, -1, -1])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Found %d reads with new small del breakpoints." %(len(new_bp_list_)))

		new_bp_clusters = cluster_bp_list(new_bp_list_, self.min_cluster_cutoff, self.max_breakpoint_distance_cutoff)
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
					if (num_subcluster == 0 and len(set(bpr)) >= self.min_cluster_cutoff) or (len(set(bpr)) >= max(self.normal_cov, 3.0)):
						io1 = interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals)
						io2 = interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals)
						if io1 >= 0 and io2 >= 0:
							assert self.amplicon_intervals[io1][3] == self.amplicon_intervals[io2][3]
							bpi = self.addbp(bp, set(bpr), bp_stats_, self.amplicon_intervals[io1][3])
							try:
								self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))].add(bpi)
							except:
								self.amplicon_interval_connections[(min(io1, io2), max(io1, io2))] = set([bpi])
					else:
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the subcluster %d." %(num_subcluster))
					num_subcluster += 1
			else:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the cluster.")


	
	def find_cn_breakpoints(self, b = 300, n = 50):
		"""
		"""
		cns_boundaries = []
		for ai in range(len(self.amplicon_intervals)):
			chr = self.amplicon_intervals[ai][0]
			s = self.amplicon_intervals[ai][1]
			e = self.amplicon_intervals[ai][2]
			si = list(self.pos2cni(chr, s))[0].data
			ei = list(self.pos2cni(chr, e))[0].data
			for i in range(si, ei):
				cns_boundaries.append((ai, chr, self.cns_intervals_by_chr[chr][i][1], \
							self.cns_intervals_by_chr[chr][i][2], self.cns_intervals_by_chr[chr][i + 1][2]))
		cb_del_list = []
		for cbi in range(len(cns_boundaries)):
			cb = cns_boundaries[cbi]
			for bpi in range(len(self.new_bp_list)):
				bp = self.new_bp_list[bpi]			
				if (bp[0] == cb[1] and bp[1] - 6001 < cb[3] < bp[1] + 6000) or \
					(bp[3] == cb[1] and bp[4] - 6001 < cb[3] < bp[4] + 6000):
					cb_del_list.append(cbi)
					break
		for cbi in range(len(cns_boundaries)):
			if cbi not in cb_del_list:
				cb = cns_boundaries[cbi]
				nl, nr = n, n
				nl = min(nl, (cb[3] - cb[2] + 1) // b)
				nr = min(nr, (cb[4] - cb[3]) // b)
				print (nl, nr, cb) 
				cov_profile = [sum([sum(nc) for nc in self.lr_bamfh.count_coverage(cb[1], cb[3] - (i + 1) * b + 1, cb[3] - i * b + 1)]) * 1.0 / b for i in range(nl)][::-1] + \
						[sum([sum(nc) for nc in self.lr_bamfh.count_coverage(cb[1], cb[3] + i * b + 1, cb[3] + (i + 1) * b + 1)]) * 1.0 / b for i in range(nr)]
				print(cov_profile)
				cb_refined = [-1, 0.0]
				for i in range(max(1, nl - 6000 // b), nl + min(nr - 1, 6000 // b)):
					dmu = np.mean(cov_profile[:i]) - np.mean(cov_profile[i:])
					if abs(dmu) > abs(cb_refined[1]):
						cb_refined[0] = i
						cb_refined[1] = dmu
				print (cb_refined)
				pval = 1.0
				cov_profile_split = [cov_profile[:cb_refined[0]], cov_profile[cb_refined[0]:]]
				if len(cov_profile_split[0]) > 1 and len(cov_profile_split[1]) > 1:
					pval = stats.ttest_ind(cov_profile_split[0], cov_profile_split[1], equal_var = False)[1]
				elif len(cov_profile_split[0]) == 1:
					zscore = abs(cov_profile_split[0][0] - np.mean(cov_profile)) / np.std(cov_profile)
					pval = stats.norm.sf(zscore)
				elif len(cov_profile_split[1]) == 1:
					zscore = abs(cov_profile_split[1][0] - np.mean(cov_profile)) / np.std(cov_profile)
					pval = stats.norm.sf(zscore)
				print ("P val = ", pval)
				if pval <= 0.01 and abs(cb_refined[1]) >= self.normal_cov:
					if cb_refined[0] < nl:
						self.source_edges.append(['source', -1, '-', cb[1], cb[3] - (nl - cb_refined[0]) * b, '+', abs(cb_refined[1])])
					else:
						self.source_edges.append(['source', -1, '-', cb[1], cb[3] + (cb_refined[0] - nl) * b, '+', abs(cb_refined[1])])
					if cb_refined[1] < 0:
						self.source_edges[-1][4] += 1
						self.source_edges[-1][5] = '-'
					self.source_edge_ccids.append(self.amplicon_intervals[cb[0]][3])
		#print (self.source_edges)
		#print (self.source_edge_ccids)
		#for bp in self.new_bp_list:
		#	print (bp[:-1])
			


	def build_graph(self):
		"""
		"""
		split_int = dict()
		for bpi in range(len(self.new_bp_list)):
			bp = self.new_bp_list[bpi]
			for ai in range(len(self.amplicon_intervals)):	
				seg = self.amplicon_intervals[ai]
				if bp[0] == seg[0] and seg[1] < bp[1] < seg[2]:
					if bp[2] == '+':
						try:
							split_int[ai].append((bp[1], bp[1] + 1, bpi, 1, '+'))
						except:
							split_int[ai] = [(bp[1], bp[1] + 1, bpi, 1, '+')]
					if bp[2] == '-':
						try:
							split_int[ai].append((bp[1] - 1, bp[1], bpi, 1, '-'))
						except:
							split_int[ai] = [(bp[1] - 1, bp[1], bpi, 1, '-')]
				if bp[3] == seg[0] and seg[1] < bp[4] < seg[2]:
					if bp[5] == '+':
						try:
							split_int[ai].append((bp[4], bp[4] + 1, bpi, 4, '+'))
						except:
							split_int[ai] = [(bp[4], bp[4] + 1, bpi, 4, '+')]
					if bp[5] == '-':
						try:
							split_int[ai].append((bp[4] - 1, bp[4], bpi, 4, '-'))
						except:
							split_int[ai] = [(bp[4] - 1, bp[4], bpi, 4, '-')]
		for srci in range(len(self.source_edges)):
			srce = self.source_edges[srci]
			for ai in range(len(self.amplicon_intervals)):	
				seg = self.amplicon_intervals[ai]
				if srce[3] == seg[0] and seg[1] < srce[4] < seg[2]:
					if srce[5] == '+':
						try:
							split_int[ai].append((srce[4], srce[4] + 1, len(self.new_bp_list) + srci, 4, '+'))
						except:
							split_int[ai] = [(srce[4], srce[4] + 1, len(self.new_bp_list) + srci, 4, '+')]
					if srce[5] == '-':
						try:
							split_int[ai].append((srce[4] - 1, srce[4], len(self.new_bp_list) + srci, 4, '-'))
						except:
							split_int[ai] = [(srce[4] - 1, srce[4], len(self.new_bp_list) + srci, 4, '-')]

		#for ai in split_int.keys():
		#	split_int[ai].sort(key = lambda item: item[0])
		"""
		for ai in range(len(self.amplicon_intervals)):
			chr = self.amplicon_intervals[ai][0]
			lcni = list(self.pos2cni(chr, self.amplicon_intervals[ai][1]))[0].data
			rcni = list(self.pos2cni(chr, self.amplicon_intervals[ai][2]))[0].data
			assert self.amplicon_intervals[ai][1] == self.cns_intervals_by_chr[chr][lcni][1]
			assert self.amplicon_intervals[ai][2] == self.cns_intervals_by_chr[chr][rcni][2]
			#print (lcni, rcni)
			if lcni != rcni:
				for cni in range (lcni, rcni):
					if ai not in split_int:
						print ('not in', self.cns_intervals_by_chr[chr][cni], self.cns_intervals_by_chr[chr][cni + 1])
					else:
						cn_bp_covered = 0
						for bp_ in split_int[ai]:
							if 0 <= self.cns_intervals_by_chr[chr][cni][2] - bp_[1] <= 5000 or 0<= bp_[0] - self.cns_intervals_by_chr[chr][cni + 1][1] <= 5000:
								cn_bp_covered = 1
								break
						if cn_bp_covered == 0:
							print ('not covered', self.cns_intervals_by_chr[chr][cni], self.cns_intervals_by_chr[chr][cni + 1])
		"""
		#source_edges = [] # Special AA discordant edges connected to a source node
		#source_edge_ccids = []
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will split the following %d amplicon intervals into sequence edges." %(len(split_int)))
		for ai in split_int.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Will split the amplicon interval at index %d." %(ai))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSplit interval at %s." %(split_int[ai]))
		for ai in split_int.keys():
			split_int[ai].sort(key = lambda item: item[0])
			#sseg = self.seq_edges[ai]
			# Cut amplicon intervals
			print (self.amplicon_intervals[ai])
			print ('split_int:', split_int[ai])
			if self.amplicon_intervals[ai][2] - self.amplicon_intervals[ai][1] > self.max_interval_cutoff:
				if split_int[ai][0][0] - self.amplicon_intervals[ai][1] > 100000:
					self.amplicon_intervals[ai][1] = split_int[ai][0][0] - 100000
				if self.amplicon_intervals[ai][2] - split_int[ai][-1][1] > 100000:
					self.amplicon_intervals[ai][2] = split_int[ai][-1][1] + 100000
				print ('corrected', self.amplicon_intervals[ai])
			
			sseg = self.amplicon_intervals[ai]
			for ssi in range(len(split_int[ai])):
				if ssi == 0:
					self.seq_edges.append([sseg[0], sseg[1], split_int[ai][ssi][0], -1, 
								int(split_int[ai][ssi][0]) - int(sseg[1]) + 1])
					self.seq_edge_ccids.append(sseg[3])
					self.concordant_edges.append([sseg[0], split_int[ai][ssi][0], '+', sseg[0], split_int[ai][ssi][1], '-', -1])
					self.concordant_edge_ccids.append(sseg[3])
				elif split_int[ai][ssi][0] > split_int[ai][ssi - 1][0]:
					self.seq_edges.append([sseg[0], split_int[ai][ssi - 1][1], split_int[ai][ssi][0], 
								-1, int(split_int[ai][ssi][0]) - int(split_int[ai][ssi - 1][1]) + 1])
					self.seq_edge_ccids.append(sseg[3])
					self.concordant_edges.append([sseg[0], split_int[ai][ssi][0], '+', sseg[0], split_int[ai][ssi][1], '-', -1])
					self.concordant_edge_ccids.append(sseg[3])
			self.seq_edges.append([sseg[0], split_int[ai][-1][1], sseg[2], -1, 
					int(sseg[2]) - int(split_int[ai][-1][1]) + 1])
			self.seq_edge_ccids.append(sseg[3])
		for ai in range(len(self.amplicon_intervals)):
			if ai not in split_int:
				sseg = self.amplicon_intervals[ai]
				self.seq_edges.append([sseg[0], sseg[1], sseg[2], -1, sseg[2] - sseg[1] + 1])
				self.seq_edge_ccids.append(sseg[3])
				#print ('---', ai, self.amplicon_intervals[ai])
		#print (self.seq_edges, self.seq_edge_ccids)
		seq_edge_idx_sorted = sorted(range(len(self.seq_edges)), key=lambda i: (global_names.chr_idx[self.seq_edges[i][0]], self.seq_edges[i][1]))
		self.seq_edges = [self.seq_edges[i] for i in seq_edge_idx_sorted]
		self.seq_edge_ccids = [self.seq_edge_ccids[i] for i in seq_edge_idx_sorted]
		concordant_edge_idx_sorted = sorted(range(len(self.concordant_edges)), key=lambda i: (global_names.chr_idx[self.concordant_edges[i][0]], self.concordant_edges[i][1]))
		self.concordant_edges = [self.concordant_edges[i] for i in concordant_edge_idx_sorted]
		self.concordant_edge_ccids = [self.concordant_edge_ccids[i] for i in concordant_edge_idx_sorted]
		
		for ai in range(len(self.amplicon_intervals)):
			self.endnodes.append((self.amplicon_intervals[ai][0], self.amplicon_intervals[ai][1], '-'))
			self.endnodes.append((self.amplicon_intervals[ai][0], self.amplicon_intervals[ai][2], '+'))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "The following nodes correspond to interval ends.")
		for node in self.endnodes:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tNode %s." %(str(node)))
		for seqi in range(len(self.seq_edge_ccids)):
			seq_edge = self.seq_edges[seqi]
			seq_edge_ccid = self.seq_edge_ccids[seqi]
			if seq_edge_ccid not in self.nodes:
				self.nodes[seq_edge_ccid] = dict()
			try:
				self.nodes[seq_edge_ccid][(seq_edge[0], seq_edge[1], '-')][0].append(seqi)
			except:
				self.nodes[seq_edge_ccid][(seq_edge[0], seq_edge[1], '-')] = [[seqi], [], [], []]
			try:
				self.nodes[seq_edge_ccid][(seq_edge[0], seq_edge[2], '+')][0].append(seqi)
			except:
				self.nodes[seq_edge_ccid][(seq_edge[0], seq_edge[2], '+')] = [[seqi], [], [], []]
		for ci in range(len(self.concordant_edges)):
			concordant_edge = self.concordant_edges[ci]
			concordant_edge_ccid = self.concordant_edge_ccids[ci]
			self.nodes[concordant_edge_ccid][(concordant_edge[0], concordant_edge[1], concordant_edge[2])][1].append(ci)
			self.nodes[concordant_edge_ccid][(concordant_edge[3], concordant_edge[4], concordant_edge[5])][1].append(ci)

		for bpi in range(len(self.new_bp_list)):
			bp = self.new_bp_list[bpi]
			bp_ccid = self.new_bp_ccids[bpi]
			io1 = interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals)
			io2 = interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals)
			assert self.amplicon_intervals[io1][3] == self.amplicon_intervals[io2][3]
			if self.amplicon_intervals[io1][3] != bp_ccid:
				bp_ccid = self.amplicon_intervals[io1][3]
			self.nodes[bp_ccid][(bp[0], bp[1], bp[2])][2].append(bpi)
			self.nodes[bp_ccid][(bp[3], bp[4], bp[5])][2].append(bpi)

		for srci in range(len(self.source_edges)):
			srce = self.source_edges[srci]
			src_ccid = self.source_edge_ccids[srci]
			self.nodes[src_ccid][(srce[3], srce[4], srce[5])][3].append(srci)

		for ccid in self.nodes.keys():
			nseq, nc, nd, nsrc = 0, 0, 0, 0 
			for seqi in range(len(self.seq_edge_ccids)):
				if self.seq_edge_ccids[seqi] == ccid:
					nseq += 1
			for ci in range(len(self.concordant_edges)):
				if self.concordant_edge_ccids[ci] == ccid:
					nc += 1
			for bpi in range(len(self.new_bp_list)):
				if self.new_bp_ccids[bpi] == ccid:
					nd += 1
			for srci in range(len(self.source_edges)):
				if self.source_edge_ccids[srci] == ccid:
					nsrc += 1

			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num sequence edges in amplicon %d = %d." %(ccid, nseq))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num concordant edges in amplicon %d = %d." %(ccid, nc))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num discordant edges in amplicon %d = %d." %(ccid, nd))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Num source edges in amplicon %d = %d." %(ccid, nsrc))


	def assign_cov(self):
		"""
		Extract the long read coverage from bam file, if missing, for each sequence edge 
		"""
		for seg in self.seq_edges:
			if seg[3] == -1:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Finding LR cov for sequence edge %s." %(seg))
				"""
				For long read, use the total number of nucleotides
				"""
				rl_list = [read for read in self.lr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if read.infer_read_length()]
				seg[3] = [len(rl_list), sum([sum(nc) for nc in self.lr_bamfh.count_coverage(seg[0], seg[1], seg[2] + 1, quality_threshold = 0, read_callback = 'nofilter')])]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "LR cov assigned for sequence edge %s." %(seg))
		"""
		Extract the long read coverage from bam file, if missing, for each concordant edge 
		"""
		for eci in range(len(self.concordant_edges)):
			ec = self.concordant_edges[eci]
			ec_ccid = self.concordant_edge_ccids[eci]
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Finding cov for concordant edge %s." %(ec))
			rls = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[0], start = ec[1], stop = ec[1] + 1)])
			rrs = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[3], start = ec[4], stop = ec[4] + 1)])
			rls1 = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[0], start = ec[1] - self.min_bp_match_cutoff_ - 1, stop = ec[1] - self.min_bp_match_cutoff_)])
			rrs1 = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[3], start = ec[4] + self.min_bp_match_cutoff_, stop = ec[4] + self.min_bp_match_cutoff_ + 1)]) 
			rbps = set([])
			for bpi in self.nodes[ec_ccid][(ec[0], ec[1], ec[2])][2]:
				rbps |= self.new_bp_list[bpi][-1]
			for bpi in self.nodes[ec_ccid][(ec[3], ec[4], ec[5])][2]:
				rbps |= self.new_bp_list[bpi][-1]
			self.concordant_edges[eci][6] = len((rls & rrs & rls1 & rrs1) - rbps)
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tLR cov assigned for concordant edge %s." %(self.concordant_edges[eci]))



	def assign_cn(self):
		"""
		Compute the maximum likelihood CN assignment on each edge
		"""
		amplicon = 1
		for ccid in self.nodes.keys():
			self.max_CN[ccid] = 0.0
			seq_edges_ccid = [seqi for seqi in range(len(self.seq_edge_ccids)) if self.seq_edge_ccids[seqi] == ccid]
			concordant_edges_ccid = [ci for ci in range(len(self.concordant_edges)) if self.concordant_edge_ccids[ci] == ccid]
			discordant_edges_ccid = [bpi for bpi in range(len(self.new_bp_list)) if self.new_bp_ccids[bpi] == ccid]
			src_edges_ccid = [srci for srci in range(len(self.source_edges)) if self.source_edge_ccids[srci] == ccid]
			lseg, lc, ld, lsrc = len(seq_edges_ccid), len(concordant_edges_ccid), len(discordant_edges_ccid), len(src_edges_ccid)

			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Adjacent list for amplicon %d:" %(amplicon))
			for node in self.nodes[ccid].keys():
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Node %s; adjacent list = %s." %(str(node), self.nodes[ccid][node]))

			nconstraints = len([node for node in self.nodes[ccid].keys() if node not in self.endnodes])
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d constraints for cvxopt." %(nconstraints))
			nvariables = lseg + lc + ld + lsrc
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d variables for cvxopt." %(nvariables))

			wcn = []
			wlncn = []
			wlrseg = []
			wcn = [0.5 * self.normal_cov * self.seq_edges[seqi][-1] for seqi in seq_edges_ccid]
			wcn += [self.normal_cov for eci in concordant_edges_ccid]
			wcn += [self.normal_cov for bpi in discordant_edges_ccid]
			wcn += [0.5 * self.normal_cov for srci in src_edges_ccid]
			wlncn = [-0.5 for seqi in seq_edges_ccid]
			wlncn += [self.concordant_edges[eci][6] * 1.0 for eci in concordant_edges_ccid]
			wlncn += [len(self.new_bp_list[edi][-1]) * 1.0 for edi in discordant_edges_ccid]
			wlncn += [-0.5 for srci in src_edges_ccid]
			wlrseg = [(0.5 * self.seq_edges[seqi][3][1] * self.seq_edges[seqi][3][1] / \
				(self.normal_cov * self.seq_edges[seqi][-1])) for seqi in seq_edges_ccid]
			wlrseg += [0.0 for eci in concordant_edges_ccid]
			wlrseg += [0.0 for bpi in discordant_edges_ccid]
			wlrseg += [(0.5 * self.source_edges[srci][-1] * self.source_edges[srci][-1] / self.normal_cov) for srci in src_edges_ccid]
			wcn = cvxopt.matrix(wcn)
			wlncn = cvxopt.matrix(wlncn)
			wlrseg = cvxopt.matrix(wlrseg)
		
			#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "The gradient of CN vector is: wcn * CN - wlncn * CN:")
			#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\twcn = %s;" %(wcn))
			#logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\twlncn = %s." %(wlncn))
		
			ci = 0
			balance_constraints = np.zeros([nconstraints, nvariables])
			for node in self.nodes[ccid].keys():
				if node not in self.endnodes:
					for seqi in self.nodes[ccid][node][0]:
						balance_constraints[ci][seq_edges_ccid.index(seqi)] = 1
					for eci in self.nodes[ccid][node][1]:
						balance_constraints[ci][lseg + concordant_edges_ccid.index(eci)] = -1
					for edi in self.nodes[ccid][node][2]:
						balance_constraints[ci][lseg + lc + discordant_edges_ccid.index(edi)] = -1
					for srci in self.nodes[ccid][node][3]:
						balance_constraints[ci][lseg + lc + ld + src_edges_ccid.index(srci)] = -1
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
				
			options = {'maxiters': 1000, 'show_progress': False}
			sol = ''
			if nconstraints > 0:
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
					for i in range(lseg):
						self.seq_edges[seq_edges_ccid[i]] += [sol['x'][i] * 2]
						if sol['x'][i] * 2 > self.max_CN[ccid]:
							self.max_CN[ccid] = sol['x'][i] * 2
					for i in range(lc):
						self.concordant_edges[concordant_edges_ccid[i]] += [sol['x'][lseg + i] * 2]
						if sol['x'][lseg + i] * 2 > self.max_CN[ccid]:
							self.max_CN[ccid] = sol['x'][lseg + i] * 2
					for i in range(ld):
						self.new_bp_list[discordant_edges_ccid[i]] += [sol['x'][lseg + lc + i] * 2]
						if sol['x'][lseg + lc + i] * 2 > self.max_CN[ccid]:
							self.max_CN[ccid] = sol['x'][lseg + lc + i] * 2
					for i in range(lsrc):
						self.source_edges[src_edges_ccid[i]] += [sol['x'][lseg + lc + ld + i] * 2]
						if sol['x'][lseg + lc + ld + i] * 2 > self.max_CN[ccid]:
							self.max_CN[ccid] = sol['x'][lseg + lc + ld + i] * 2
			else:
				assert lc == 0 and ld == 0 and lsrc == 0
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Skipped convex optimization.")
				for i in range(lseg):
					sseg = self.seq_edges[seq_edges_ccid[i]]
					cn_segi = self.normal_cov * sseg[-1] * 2.0 / sseg[3][1]
					self.seq_edges[seq_edges_ccid[i]] += [cn_segi]
					if cn_segi > self.max_CN[ccid]:
						self.max_CN[ccid] = cn_segi
			self.max_CN[ccid] += 1.0
			amplicon += 1


	
	def output_breakpoint_info(self, obpfile):
		"""
		Write the list of breakpoints to file
		"""
		with open(obpfile, 'w') as fp:
			fp.write("chr1\tpos1\tchr2\tpos2\torientation\tlr_support\tlr_info=[avg1, avg2, std1, std2, mapq1, mapq2]\n")
			for bpi in range(len(self.new_bp_list)):
				bp = self.new_bp_list[bpi]
				bp_stats = self.new_bp_stats[bpi]
				fp.write("%s\t%s\t%s\t%s\t%s%s\t-1\t%d\t%s\n" 
					%(bp[3], bp[4], bp[0], bp[1], bp[5], bp[2], len(bp[-1]), bp_stats))


	def output_breakpoint_graph(self, ogfile_prefix):
		"""
		Write a breakpoint graph to file in AA graph format
		"""
		amplicon = 1
		for ccid in self.nodes.keys():
			with open(ogfile_prefix + "_amplicon" + str(amplicon) + "_graph.txt", 'w') as fp:
				fp.write("SequenceEdge: StartPosition, EndPosition, PredictedCN, NumberOfLongReads, Size\n")
				for seqi in range(len(self.seq_edges)):
					if self.seq_edge_ccids[seqi] == ccid:
						sseg = self.seq_edges[seqi]
						fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%s\t%s\n" %(sseg[0], sseg[1], sseg[0], sseg[2], sseg[-1], sseg[3][0], sseg[4]))
				fp.write("BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads\n")
				#fp.write("HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence\n")
				for srci in range(len(self.source_edges)):
					if self.source_edge_ccids[srci] == ccid:
						es = self.source_edges[srci]
						fp.write("source\t%s:%s%s->%s:%s%s\t%f\t-1\n" %(es[0], es[1], es[2], es[3], 
							es[4], es[5], es[-1]))
				for eci in range(len(self.concordant_edges)):
					if self.concordant_edge_ccids[eci] == ccid:
						ec = self.concordant_edges[eci]
						fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\n" %(ec[0], ec[1], ec[2], ec[3], 
							ec[4], ec[5], ec[-1], ec[6]))
				for bpi in range(len(self.new_bp_list)):
					if self.new_bp_ccids[bpi] == ccid:
						ed = self.new_bp_list[bpi]
						fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t%d\n" %(ed[0], ed[1], ed[2], ed[3], 
							ed[4], ed[5], ed[-1], len(ed[-2])))
			amplicon += 1




	def closebam(self):
		"""
		Close the short read and long read bam file
		"""
		self.lr_bamfh.close()



if __name__ == '__main__':
	global_names.TSTART = start_time
	parser = argparse.ArgumentParser(description = "Examine ")
	parser.add_argument("--lr_bam", help = "Sorted indexed long read bam file.", required = True)
	parser.add_argument("--seed", help = "File including seed intervals.", required = True)
	parser.add_argument("--output_prefix", help = ".", required = True)
	parser.add_argument("--output_bp", help = "If specified, only output the list of breakpoints.",  action = 'store_true')
	parser.add_argument("--lr_cnseg", help = "Long read *.cns file.", required = True)
	args = parser.parse_args()

	logging.basicConfig(filename = 'infer_breakpoint_graph.log', filemode = 'w', level = logging.DEBUG, 
			format = '[%(name)s:%(levelname)s]\t%(message)s')
	logging.info("Python version " + sys.version + "\n")
	commandstring = 'Commandline: '
	for arg in sys.argv:
		if ' ' in arg:
			commandstring += '"{}" '.format(arg)
		else:
			commandstring += "{} ".format(arg)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + commandstring)

	b2bn = bam_to_breakpoint_nanopore(args.lr_bam, args.seed)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Opened LR bam files.")
	b2bn.read_cns(args.lr_cnseg)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed parsing CN segment files.")
	b2bn.fetch()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed fetching reads containing breakpoints.")
	b2bn.hash_alignment_to_seg()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed hashing chimeric reads to CN segments.")
	b2bn.find_amplicon_intervals()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed finding amplicon intervals.")
	b2bn.find_smalldel_breakpoints()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed finding small del breakpoints.")
	b2bn.find_breakpoints()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed finding all breakpoints.")
	b2bn.find_cn_breakpoints()
	b2bn.build_graph()
	b2bn.assign_cov()
	b2bn.assign_cn()
	#b2bn.output_breakpoint_info(args.output_fn)
	
	b2bn.output_breakpoint_graph(args.output_prefix)
	b2bn.closebam()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total runtime.")
	#print (len(b2bn.new_bp_list), b2bn.new_bp_list)



