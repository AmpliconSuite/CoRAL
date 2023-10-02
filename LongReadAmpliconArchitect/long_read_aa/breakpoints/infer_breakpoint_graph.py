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
from scipy import stats
import cvxopt
import cvxopt.modeling
import intervaltree
import gurobipy as gp
from gurobipy import GRB


import cigar_parsing
from breakpoint_utilities import *
from long_read_aa.breakpoints import global_names
global_names.TSTART = start_time


edge_type_to_index = {'s': 0, 'c': 1, 'd': 2}


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
	ccid2id = dict()
	endnodes = [] 
	small_del_indices = [] # Indices of +- breakpoints on the same chr and < small_del_cutoff in discordant_edges
	source_edges = [] # Special AA discordant edges connected to a source node
	source_edge_ccids = []
	max_CN = dict()
	max_seq_repeat = dict()

	concordant_edges_reads = []
	path_constraints = [[], [], []]
	valid_path_constraints = dict()
	cycles = dict() # cycles, paths
	cycle_weights = dict() # cycles, paths
	path_constraints_satisfied = dict() # cycles, paths


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
			#if r == "amplicon-3_18538_aligned_27073_F_36_52065_12":
			#	print (self.chimeric_alignments[r])
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
			# Fix 09/28/23 - added flanking region while searching for amplicon intervals
			self.amplicon_intervals[ai][1] = self.cns_intervals_by_chr[chr][lcni][1]
			if len(list(self.pos2cni(chr, self.cns_intervals_by_chr[chr][lcni][1] - self.interval_delta))) > 0:
				self.amplicon_intervals[ai][1] = self.cns_intervals_by_chr[chr][lcni][1] - self.interval_delta
			self.amplicon_intervals[ai][2] = self.cns_intervals_by_chr[chr][rcni][2]
			if len(list(self.pos2cni(chr, self.cns_intervals_by_chr[chr][rcni][2] + self.interval_delta))) > 0:
				self.amplicon_intervals[ai][2] = self.cns_intervals_by_chr[chr][rcni][2] + self.interval_delta
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tUpdated amplicon interval %s" %self.amplicon_intervals[ai])
 
		ccid = 0
		for ai in range(len(self.amplicon_intervals)):
			if self.amplicon_intervals[ai][3] == -1:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Begin processing amplicon interval %d" %ai)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAmplicon interval %s" %self.amplicon_intervals[ai])
				# Fix 09/28/23 - added flanking region while searching for amplicon intervals
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
			# Fix 09/28/23 - added flanking region while searching for amplicon intervals
			if not (interval_adjacent(amplicon_intervals_sorted[ai + 1], amplicon_intervals_sorted[ai]) or \
				interval_overlap(amplicon_intervals_sorted[ai], amplicon_intervals_sorted[ai + 1])):
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
				if bp[0] == bp[3] and bp[2] == '-' and bp[5] == '+' and abs(bp[1] - bp[4]) < self.small_del_cutoff:
					self.small_del_indices.append(bpi)
				return bpi
		bpi = len(self.new_bp_list)
		self.new_bp_list.append(bp_ + [bpr_])
		self.new_bp_ccids.append(ccid)
		self.new_bp_stats.append(bp_stats_)
		if bp_[0] == bp_[3] and bp_[2] == '-' and bp_[5] == '+' and abs(bp_[1] - bp_[4]) < self.small_del_cutoff:
			self.small_del_indices.append(bpi)
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
								new_bp_list.append(interval2bp(rr_int[ri], rr_int[ri + 1], (r, ri, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri][1])) + [q_[ri], q_[ri + 1]])
								bassigned[ri] = 1
								#bp_flag_r = 1
							elif interval_overlap(rr_int[ri + 1], [nint_[0], ns, ne]) and \
								interval_overlap(rr_int[ri], self.amplicon_intervals[ai_]) and q_[ri] >= 20 and q_[ri + 1] >= 20:
								new_bp_list.append(interval2bp(rr_int[ri], rr_int[ri + 1], (r, ri, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri][1])) + [q_[ri], q_[ri + 1]])
								bassigned[ri] = 1

						# Breakpoint from local alignment i - 1 and i + 1
						for ri in range(1, len(rr_int) - 1):
							if bassigned[ri - 1] == 0 and bassigned[ri] == 0 and q_[ri] < 10 and q_[ri - 1] >= 20 and q_[ri + 1] >= 20 and \
								interval_overlap(rr_int[ri - 1], [nint_[0], ns, ne]) and \
								interval_overlap(rr_int[ri + 1], self.amplicon_intervals[ai_]):
								new_bp_list.append(interval2bp(rr_int[ri - 1], rr_int[ri + 1], (r, ri - 1, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri - 1][1])) + [q_[ri - 1], q_[ri + 1]])
							elif bassigned[ri - 1] == 0 and bassigned[ri] == 0 and q_[ri] < 10 and q_[ri - 1] >= 20 and q_[ri + 1] >= 20 and \
								interval_overlap(rr_int[ri + 1], [nint_[0], ns, ne]) and \
								interval_overlap(rr_int[ri - 1], self.amplicon_intervals[ai_]):
								new_bp_list.append(interval2bp(rr_int[ri - 1], rr_int[ri + 1], (r, ri - 1, ri + 1), int(r_int[ri + 1][0]) - int(r_int[ri - 1][1])) + [q_[ri - 1], q_[ri + 1]])

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
								l = max(self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1] - self.interval_delta, 
									self.cns_intervals_by_chr[chr_][0][1])
								r = min(lir + self.interval_delta, self.cns_intervals_by_chr[chr_][-1][2])
								if nint_segs[lasti][1] - int(self.max_interval_cutoff / 2) > l:
									l = nint_segs[lasti][1] - int(self.max_interval_cutoff / 2)
								if nint_segs[i][1] + int(self.max_interval_cutoff / 2) < r:
									r = nint_segs[i][1] + int(self.max_interval_cutoff / 2)
								if len(list(self.pos2cni(chr_, l))) == 0:
									l = self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1]
								if len(list(self.pos2cni(chr_, r))) == 0:
									r = lir
								new_intervals_refined.append([chr_, l, r, -1])
								new_intervals_connections.append([])
								for i_ in range(lasti, i + 1):
									new_intervals_connections[-1].append(nint_segs[i_][2])
								lasti = i + 1
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFixed new interval: %s." %new_intervals_refined[-1])
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tList of breakpoints connected to the new interval:")
								for bpi in new_intervals_connections[-1]:
									logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\t\t%s" %self.new_bp_list[bpi][:6])
						if len(nint_segs) > 0:
							l = max(self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1] - self.interval_delta,
								self.cns_intervals_by_chr[chr_][0][1])
							r = min(self.cns_intervals_by_chr[chr_][nint_segs[-1][0]][2] + self.interval_delta,
								self.cns_intervals_by_chr[chr_][-1][2])
							if nint_segs[lasti][1] - int(self.max_interval_cutoff / 2) > l:
								l = nint_segs[lasti][1] - int(self.max_interval_cutoff / 2) > l
							if nint_segs[-1][1] + int(self.max_interval_cutoff / 2) < r:
								r = nint_segs[-1][1] + int(self.max_interval_cutoff / 2)
							if len(list(self.pos2cni(chr_, l))) == 0:
								l = self.cns_intervals_by_chr[chr_][nint_segs[lasti][0]][1]
							if len(list(self.pos2cni(chr_, r))) == 0:
								r = self.cns_intervals_by_chr[chr_][nint_segs[-1][0]][2]
							new_intervals_refined.append([chr_, l, r, -1])
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
								l = max(self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1] - self.interval_delta,
									self.cns_intervals_by_chr[nint_segs_[lasti][0]][0][1])
								r = min(lir + self.interval_delta,
									self.cns_intervals_by_chr[nint_segs_[i][0]][-1][2])
								if nint_segs_[lasti][2] - int(self.max_interval_cutoff / 2) > l:
									l = nint_segs_[lasti][2] - int(self.max_interval_cutoff / 2)
								if nint_segs_[i][2] + int(self.max_interval_cutoff / 2) < r:
									r = nint_segs_[i][2] + int(self.max_interval_cutoff / 2)
								if len(list(self.pos2cni(nint_segs_[lasti][0], l))) == 0:
									l = self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1]
								if len(list(self.pos2cni(nint_segs_[i][0], r))) == 0:
									r = lir
								new_intervals_refined.append([nint_segs_[lasti][0], l, r, -1])
								new_intervals_connections.append([])
								lasti = i + 1
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tFixed new interval: %s." %new_intervals_refined[-1])
								logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\t\tSkip breakpoints connected to the new interval.")
						if len(nint_segs_) > 0:
							l = max(self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1] - self.interval_delta,
								self.cns_intervals_by_chr[nint_segs_[lasti][0]][0][1])
							r = min(self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[-1][1]][2] + self.interval_delta,
								self.cns_intervals_by_chr[nint_segs_[lasti][0]][-1][2])
							if nint_segs_[lasti][2] - int(self.max_interval_cutoff / 2) > l:
								l = nint_segs_[lasti][2] - int(self.max_interval_cutoff / 2)
							if nint_segs_[-1][2] + int(self.max_interval_cutoff / 2) < r:
								r = nint_segs_[-1][2] + int(self.max_interval_cutoff / 2)
							if len(list(self.pos2cni(nint_segs_[lasti][0], l))) == 0:
								l = self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[lasti][1]][1]
							if len(list(self.pos2cni(nint_segs_[lasti][0], r))) == 0:
								r = self.cns_intervals_by_chr[nint_segs_[lasti][0]][nint_segs_[-1][1]][2]
							new_intervals_refined.append([nint_segs_[lasti][0], l, r, -1])
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
							new_bp_list_.append(interval2bp(rr_int[i], rr_int[i + 1], (r, i, i + 1), int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '+':
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list_.append(interval2bp(rr_int[i], rr_int[i + 1], (r, i, i + 1), int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
							bassigned[i] = 1
					elif rr_int[i + 1][3] == '-':
						gr =  int(r_int[i + 1][0]) - int(r_int[i][1])
						grr = int(rr_int[i][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)) and q_[i] >= 20 and q_[i + 1] >= 20:
							new_bp_list_.append(interval2bp(rr_int[i], rr_int[i + 1], (r, i, i + 1), int(r_int[i + 1][0]) - int(r_int[i][1])) + [q_[i], q_[i + 1]])
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
						new_bp_list_.append(interval2bp(rr_int[i - 1], rr_int[i + 1], (r, i - 1, i + 1), int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '+':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i + 1][1]) - int(rr_int[i - 1][2])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list_.append(interval2bp(rr_int[i - 1], rr_int[i + 1], (r, i - 1, i + 1), int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
					elif rr_int[i + 1][3] == '-':
						gr =  int(r_int[i + 1][0]) - int(r_int[i - 1][1])
						grr = int(rr_int[i - 1][2]) - int(rr_int[i + 1][1])
						if abs(gr - grr) > max(100, abs(gr * 0.2)):
							new_bp_list_.append(interval2bp(rr_int[i - 1], rr_int[i + 1], (r, i - 1, i + 1), int(r_int[i + 1][0]) - int(r_int[i - 1][1])) + [q_[i - 1], q_[i + 1]])
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
				rq = read.mapping_quality
				if rq < 20:
					continue # Fixed on 04/09/23 - Also introduced quality control on small del breakpoints
				blocks = read.get_blocks()
				for bi in range(len(blocks) - 1):
					if abs(blocks[bi + 1][0] - blocks[bi][1]) > self.min_del_len:
						try:
							self.large_indel_alignments[rn].append([ai[0], blocks[bi + 1][0], blocks[bi][1], blocks[0][0], blocks[-1][1], rq])
						except:
							self.large_indel_alignments[rn] = [[ai[0], blocks[bi + 1][0], blocks[bi][1], blocks[0][0], blocks[-1][1], rq]]
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Fetched %d reads with large indels in CIGAR." %(len(self.large_indel_alignments)))
		
		for r in self.large_indel_alignments.keys():
			for rr_gap_i in range(len(self.large_indel_alignments[r])):
				rr_gap = self.large_indel_alignments[r][rr_gap_i][:3]
				rr_gap_ = rr_gap
				if rr_gap[2] > rr_gap[1]:
					rr_gap_[2] = rr_gap[1]
					rr_gap_[1] = rr_gap[2]
				new_bp_list_.append([rr_gap_[0], rr_gap_[1], '-', rr_gap_[0], rr_gap_[2], '+', (r, rr_gap_i, rr_gap_i), 0, 0, -1, -1])
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
				#print (nl, nr, cb) 
				cov_profile = [sum([sum(nc) for nc in self.lr_bamfh.count_coverage(cb[1], cb[3] - (i + 1) * b + 1, cb[3] - i * b + 1)]) * 1.0 / b for i in range(nl)][::-1] + \
						[sum([sum(nc) for nc in self.lr_bamfh.count_coverage(cb[1], cb[3] + i * b + 1, cb[3] + (i + 1) * b + 1)]) * 1.0 / b for i in range(nr)]
				#print(cov_profile)
				cb_refined = [-1, 0.0]
				for i in range(max(1, nl - 6000 // b), nl + min(nr - 1, 6000 // b)):
					dmu = np.mean(cov_profile[:i]) - np.mean(cov_profile[i:])
					if abs(dmu) > abs(cb_refined[1]):
						cb_refined[0] = i
						cb_refined[1] = dmu
				#print (cb_refined)
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
				#print ("P val = ", pval)
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
			#print (self.amplicon_intervals[ai])
			#print ('split_int:', split_int[ai])
			if self.amplicon_intervals[ai][2] - self.amplicon_intervals[ai][1] > self.max_interval_cutoff:
				if split_int[ai][0][0] - self.amplicon_intervals[ai][1] > 100000:
					self.amplicon_intervals[ai][1] = split_int[ai][0][0] - 100000
				if self.amplicon_intervals[ai][2] - split_int[ai][-1][1] > 100000:
					self.amplicon_intervals[ai][2] = split_int[ai][-1][1] + 100000
				#print ('corrected', self.amplicon_intervals[ai])
			
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
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset the ccid for breakpoint %s at index %d from %d to %d." \
						%(bp[:6], bpi, bp_ccid, self.amplicon_intervals[io1][3]))
				bp_ccid = self.amplicon_intervals[io1][3]
				self.new_bp_ccids[bpi] = self.amplicon_intervals[io1][3]
			self.nodes[bp_ccid][(bp[0], bp[1], bp[2])][2].append(bpi)
			self.nodes[bp_ccid][(bp[3], bp[4], bp[5])][2].append(bpi)

		for srci in range(len(self.source_edges)):
			srce = self.source_edges[srci]
			src_ccid = self.source_edge_ccids[srci]
			self.nodes[src_ccid][(srce[3], srce[4], srce[5])][3].append(srci)

		amplicon_id = 1
		for ccid in self.nodes.keys():
			self.ccid2id[ccid] = amplicon_id
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
			amplicon_id += 1


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
			self.concordant_edges_reads.append(rls | rrs)
			self.concordant_edges[eci][6] = len((rls & rrs & rls1 & rrs1) - rbps)
			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tLR cov assigned for concordant edge %s." %(self.concordant_edges[eci]))



	def assign_cn(self):
		"""
		Compute the maximum likelihood CN assignment on each edge
		"""
		for ccid in self.nodes.keys():
			self.max_CN[ccid] = 0.0
			seq_edges_ccid = [seqi for seqi in range(len(self.seq_edge_ccids)) if self.seq_edge_ccids[seqi] == ccid]
			concordant_edges_ccid = [ci for ci in range(len(self.concordant_edges)) if self.concordant_edge_ccids[ci] == ccid]
			discordant_edges_ccid = [bpi for bpi in range(len(self.new_bp_list)) if self.new_bp_ccids[bpi] == ccid]
			src_edges_ccid = [srci for srci in range(len(self.source_edges)) if self.source_edge_ccids[srci] == ccid]
			lseg, lc, ld, lsrc = len(seq_edges_ccid), len(concordant_edges_ccid), len(discordant_edges_ccid), len(src_edges_ccid)

			logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Adjacent list for amplicon %d, ccid = %d:" %(self.ccid2id[ccid], ccid))
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
				seq_cn_list = [self.seq_edges[segi][-1] for segi in seq_edges_ccid if self.seq_edges[segi][-2] >= 10000]
				if len(seq_cn_list) > 0 and max(seq_cn_list) >= 5.0:
					self.max_seq_repeat[ccid] = int(round(max(seq_cn_list) / \
							np.average([self.seq_edges[segi][-1] for segi in seq_edges_ccid if self.seq_edges[segi][-1] >= 5.0], \
							weights = [self.seq_edges[segi][-2] for segi in seq_edges_ccid if self.seq_edges[segi][-1] >= 5.0]))) + 1
				else:
					self.max_seq_repeat[ccid] = 2

			else:
				assert lc == 0 and ld == 0 and lsrc == 0
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Skipped convex optimization.")
				for i in range(lseg):
					sseg = self.seq_edges[seq_edges_ccid[i]]
					cn_segi = sseg[3][1] * 2.0 / (self.normal_cov * sseg[-1]) 
					self.seq_edges[seq_edges_ccid[i]] += [cn_segi]
					if cn_segi > self.max_CN[ccid]:
						self.max_CN[ccid] = cn_segi
				self.max_seq_repeat[ccid] = 2
			self.max_CN[ccid] += 1.0


	
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
		for ccid in self.nodes.keys():
			with open(ogfile_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "_graph.txt", 'w') as fp:
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


	
	def valid_path(self, adj_list, path):
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
					if e1[1] not in adj_list[path[i]][edge_type_to_index[e1[0]]]:
						return False
					if e2[1] not in adj_list[path[i]][edge_type_to_index[e2[0]]]:
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
		ccid1 = self.seq_edge_ccids[segi0]
		segi0 = seq_edge_list[-1]
		node2 = (self.seq_edges[segi0][0], self.seq_edges[segi0][2], '+')
		ccid2 = self.seq_edge_ccids[segi0]
		assert ccid1 == ccid2
		#print (rint, node1, node2, seq_edge_list)
		path_ = self.traverse_through_sequence_edge(self.nodes[ccid1], node1, node2)[1:-1] 
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


	def chimeric_alignment_to_path_i(self, rints, ai1, ai2, bpi):
		path_ = [('d', bpi)]
		node1 = (self.new_bp_list[bpi][0], self.new_bp_list[bpi][1], self.new_bp_list[bpi][2])
		node2 = (self.new_bp_list[bpi][3], self.new_bp_list[bpi][4], self.new_bp_list[bpi][5])
		if ai1 > ai2:
			#node1, node2 = node2, node1
			path_ = self.chimeric_alignment_to_path_l(rints, ai2, node2) + path_ + self.chimeric_alignment_to_path_r(rints, ai1, node1)
			#print (node1, node2, path_)
		else:
			path_ = self.chimeric_alignment_to_path_l(rints, ai1, node1) + path_ + self.chimeric_alignment_to_path_r(rints, ai2, node2)
			#print (node1, node2, path_)
		return path_


	def traverse_through_sequence_edge(self, adj_list, start_node, end_node):
		assert start_node[2] != end_node[2]
		path_ = [start_node]
		seqi = adj_list[start_node][0][0]
		seq_edge = self.seq_edges[seqi]
		next_end = (seq_edge[0], seq_edge[1], '-')
		if start_node[2] == '-':
			next_end = (seq_edge[0], seq_edge[2], '+')
		path_.append(('s', seqi))
		path_.append(next_end)
		while next_end != end_node:
			try:
				ci = adj_list[next_end][1][0] # 07/20/2023 - read may span two amplicon intervals
			except:
				return path_
			path_.append(('c', ci))
			cedge = self.concordant_edges[ci]
			next_start = (cedge[0], cedge[1], cedge[2])
			if next_start == next_end:
				next_start = (cedge[3], cedge[4], cedge[5])
			path_.append(next_start)
			seqi = adj_list[next_start][0][0]
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
			bpi = bp_list[i]
			node1 = (self.new_bp_list[bpi][0], self.new_bp_list[bpi][1], self.new_bp_list[bpi][2])
			node2 = (self.new_bp_list[bpi][3], self.new_bp_list[bpi][4], self.new_bp_list[bpi][5])
			ccid = self.new_bp_ccids[bpi]
			if ai_list[i][0] > ai_list[i][1]:
				if i == 0:
					path_ = self.chimeric_alignment_to_path_l(rints, ai_list[i][1], node2) + [('d', bp_list[i])]
					lastnode = node1
				else:
					path_ += self.traverse_through_sequence_edge(self.nodes[ccid], lastnode, node2)
					path_.append(('d', bp_list[i]))
					lastnode = node1
					if i == len(bp_list) - 1:
						path_ += self.chimeric_alignment_to_path_r(rints, ai_list[i][0], node1)
			else:
				if i == 0:
					path_ = self.chimeric_alignment_to_path_l(rints, ai_list[i][0], node1) + [('d', bp_list[i])]
					lastnode = node2
				else:
					path_ += self.traverse_through_sequence_edge(self.nodes[ccid], lastnode, node1)
					path_.append(('d', bp_list[i]))
					lastnode = node2
					if i == len(bp_list) - 1:
						path_ += self.chimeric_alignment_to_path_r(rints, ai_list[i][1], node2)
		#print (path_)
		return path_


	def compute_path_constraints(self):
		lseg = len(self.seq_edges)
		lc = len(self.concordant_edges)
		ld = len(self.new_bp_list)

		bp_reads = dict()
		concordant_reads = dict()
		for bpi in range(ld):
			for r_ in self.new_bp_list[bpi][-2]:
				if r_[1] == r_[2]:
					if r_[0] in bp_reads:
						bp_reads[r_[0]][1].append([r_[1], r_[2], bpi])
					else:
						bp_reads[r_[0]] = [[], [[r_[1], r_[2], bpi]]]
				else:
					if r_[0] in bp_reads:
						bp_reads[r_[0]][0].append([r_[1], r_[2], bpi])
					else:
						bp_reads[r_[0]] = [[[r_[1], r_[2], bpi]], []]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d reads in total covering at least one breakpoint." %len(bp_reads))

		for rn in bp_reads.keys():
			bp_reads_rn = bp_reads[rn][0]
			bp_reads_rn_sdel = bp_reads[rn][1]
			path_ccids = []
			paths = []
			if len(bp_reads_rn) == 1 and len(bp_reads_rn_sdel) == 0:
				rints = [aint[:4] for aint in self.chimeric_alignments[rn][1]]
				ai1 = bp_reads_rn[0][0]
				ai2 = bp_reads_rn[0][1]
				bpi = bp_reads_rn[0][2]
				ccid = self.new_bp_ccids[bpi]
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "Read %s covers a single breakpoint." %rn)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bp = (%d, %d, %d)" \
						%(rints, self.chimeric_alignments[rn][2], ai1, ai2, bpi))
				path = self.chimeric_alignment_to_path_i(rints, ai1, ai2, bpi)
				paths.append(path)
				path_ccids.append(ccid)
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
					rints = [aint[:4] for aint in self.chimeric_alignments[rn][1]]
					ai_list = [bp_reads_rn[bi][:2] for bi in ai_block]
					bp_list = [bp_reads_rn[bi][2] for bi in ai_block]
					ccid = self.new_bp_ccids[bp_list[0]]
					if len(bp_list) > 1:
						for bpi in bp_list:
							assert self.new_bp_ccids[bpi] == ccid
					if len(set(bp_list)) < len(bp_list):
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the block due to repeated breakpoints.")
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %ai_block)
						continue
					path = self.chimeric_alignment_to_path(rints, ai_list, bp_list)
					paths.append(path)
					path_ccids.append(ccid)
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
				ccid = self.new_bp_ccids[bpi]
				path_ccids.append(ccid)
				path = []
				if rints[0][3] == '+':
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bp = (1, 0, %d)" \
							%(rints, rq, bpi))
					path = self.chimeric_alignment_to_path_i(rints, 1, 0, bpi)
					paths.append(path)
				else:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bp = (0, 1, %d)" \
							%(rints, rq, bpi))
					path = self.chimeric_alignment_to_path_i(rints, 0, 1, bpi)
					paths.append(path)
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
					ccid = self.new_bp_ccids[bp_list[0]]
					if len(bp_list) > 1:
						for bpi in bp_list:
							assert self.new_bp_ccids[bpi] == ccid
					if len(set(bp_list)) < len(bp_list):
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the block due to repeated breakpoints.")
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %ai_block)
						continue
					path = self.chimeric_alignment_to_path(rints_, ai_list, bp_list)
					paths.append(path)
					path_ccids.append(ccid)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq = %s; bps = %s" \
						%(rints_, rq, bp_reads_rn_sdel))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tResulting subpath = %s" %path)
			else:
				rints = [aint[:4] for aint in self.chimeric_alignments[rn][1]]
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
					ccid = self.new_bp_ccids[bp_list[0]]
					if len(bp_list) > 1:
						for bpi in bp_list:
							assert self.new_bp_ccids[bpi] == ccid
					if len(set(bp_list)) < len(bp_list):
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tDiscarded the block due to repeated breakpoints.")
						logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tBlocks of local alignments: %s" %ai_block)
						continue
					path = self.chimeric_alignment_to_path(rints, ai_list, bp_list)
					paths.append(path)
					path_ccids.append(ccid)
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tAlignment intervals on reference = %s; mapq (unsplit) = %s; bps = %s" \
						%(rints, self.chimeric_alignments[rn][2], bp_reads_rn))
					logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tResulting subpath = %s" %path)
			for pi in range(len(paths)):
				path = paths[pi]
				path_ccid = path_ccids[pi]
				#print (pi, path, path_ccid)
				if len(path) > 5 and self.valid_path(self.nodes[path_ccid], path):
					#print ('aaa')
					if path not in self.path_constraints[0]:
						self.path_constraints[0].append(path)
						self.path_constraints[1].append(1)
						self.path_constraints[2].append(path_ccid)
					else:
						pci = self.path_constraints[0].index(path)
						self.path_constraints[1][pci] += 1
						assert path_ccid == self.path_constraints[2][pci]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d distinct subpaths due to reads involving breakpoints." %len(self.path_constraints[0]))
		#extract reads in concordant_edges_reads
		for ci in range(lc):
			for rn in self.concordant_edges_reads[ci]:
				if rn not in self.large_indel_alignments and rn not in self.chimeric_alignments:
					concordant_reads[rn] = self.concordant_edge_ccids[ci]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d concordant reads within amplicon intervals." %len(concordant_reads))
		for aint in self.amplicon_intervals:
			for read in self.lr_bamfh.fetch(aint[0], aint[1], aint[2] + 1):
				rn = read.query_name
				q = read.mapq
				if q >= 20 and rn in concordant_reads:
					path = self.alignment_to_path([read.reference_name, read.reference_start, read.reference_end])
					#print ("case c:", [read.reference_name, read.reference_start, read.reference_end], path)
					#print ()
					if len(path) > 5 and self.valid_path(self.nodes[concordant_reads[rn]], path):
						if path not in self.path_constraints[0]:
							self.path_constraints[0].append(path)
							self.path_constraints[1].append(1)
							self.path_constraints[2].append(concordant_reads[rn])
						else:
							pci = self.path_constraints[0].index(path)
							self.path_constraints[1][pci] += 1
							assert concordant_reads[rn] == self.path_constraints[2][pci]
		logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "There are %d distinct subpaths in total." %len(self.path_constraints[0]))
		#print (self.small_del_indices)
		#print (self.small_del_indices_)
		#os.abort()


	def minimize_cycles(self, ccid, k, total_weights, node_order, pc_list, pc_indices,
				max_seq_repeat = 2, p_total_weight = 0.9, p_bp_cn = 0.9, num_threads = 16, of_prefix = ""):
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Regular cycle decomposition with at most %d cycles/paths allowed." %k)
		seq_edges_ccid = [seqi for seqi in range(len(self.seq_edge_ccids)) if self.seq_edge_ccids[seqi] == ccid]
		concordant_edges_ccid = [ci for ci in range(len(self.concordant_edges)) if self.concordant_edge_ccids[ci] == ccid]
		discordant_edges_ccid = [bpi for bpi in range(len(self.new_bp_list)) if self.new_bp_ccids[bpi] == ccid]
		src_edges_ccid = [srci for srci in range(len(self.source_edges)) if self.source_edge_ccids[srci] == ccid]
		endnodes_ccid = [endnode for endnode in self.endnodes if endnode in self.nodes[ccid]]
		lseg, lc, ld, lsrc = len(seq_edges_ccid), len(concordant_edges_ccid), len(discordant_edges_ccid), len(src_edges_ccid)

		nnodes = len(self.nodes[ccid])
		nedges = lseg + lc + ld + 2 * lsrc + 2 * len(endnodes_ccid)

		# Gurobi model
		m = gp.Model(of_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "_cycle_decomposition_k=" + str(k))
		
		# z[i]: indicating whether cycle or path i exists
		z = m.addVars(k, vtype = GRB.BINARY, name = ["z" + str(i) for i in range(k)])
		
		# w[i]: the weight of cycle or path i, continuous variable
		w = m.addVars(k, lb = 0.0, ub = self.max_CN[ccid], vtype = GRB.CONTINUOUS, name = ["w" + str(i) for i in range(k)])
		
		# Relationship between w[i] and z[i]
		for i in range(k):
			m.addConstr(w[i] <= z[i] * self.max_CN[ccid])

		# x: the number of times an edge occur in cycle or path i
		x_names = []
		for ei in range(nedges):
			for i in range(k):
				x_names.append("x" + str(ei) + "," + str(i))
		x = m.addVars(k * nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)
	
		# Objective: minimize the total number of cycles
		obj = gp.QuadExpr(1.0)
		for i in range(k):
			obj += z[i]
			for seqi in seq_edges_ccid:
				obj -= (x[seq_edges_ccid.index(seqi) * k + i] * w[i] * self.seq_edges[seqi][-2] / total_weights)
		m.setObjective(obj, GRB.MINIMIZE)

		# Must include at least 0.9 * total CN weights
		total_weights_expr = gp.QuadExpr(0.0)
		for i in range(k):
			for seqi in seq_edges_ccid:
				total_weights_expr += (x[seq_edges_ccid.index(seqi) * k + i] * w[i] * self.seq_edges[seqi][-2])
		m.addConstr(total_weights_expr >= p_total_weight * total_weights)

		# Eulerian constraint
		for node in self.nodes[ccid].keys():
			if node in self.endnodes:
				for i in range(k):
					m.addConstr(x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] + \
							x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] \
							== x[seq_edges_ccid.index(self.nodes[ccid][node][0][0]) * k + i])
			else:
				for i in range(k):
					ec_expr = gp.LinExpr(0.0)
					for seqi in self.nodes[ccid][node][0]:
						ec_expr += x[seq_edges_ccid.index(seqi) * k + i]
					for ci in self.nodes[ccid][node][1]:
						ec_expr -= x[(lseg + concordant_edges_ccid.index(ci)) * k + i]
					for di in self.nodes[ccid][node][2]:
						# Fix 10/01
						#dedge = self.new_bp_list[di]
						#if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]:
						#	ec_expr -= 2 * x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
						#else:
						ec_expr -= x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
					for srci in self.nodes[ccid][node][3]:
						ec_expr -= x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i] # connected to s
						ec_expr -= x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + k + i] # connected to t
					m.addConstr(ec_expr == 0.0)
		for i in range(k):
			path_expr = gp.LinExpr(0.0)
			for enodei in range(len(endnodes_ccid)):
				path_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + i] # (s, v)
				path_expr -= x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + k + i] # (v, t)
			for srci in range(lsrc): 
				path_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
				path_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # (v, t)
			m.addConstr(path_expr == 0.0)

		# CN constraint
		for seqi in range(lseg):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[seqi * k + i]
			m.addQConstr(cn_expr <= self.seq_edges[seq_edges_ccid[seqi]][-1])
		for ci in range(lc):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + ci) * k + i]
			m.addQConstr(cn_expr <= self.concordant_edges[concordant_edges_ccid[ci]][-1])
		for bpi in range(ld):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + lc + bpi) * k + i]
			m.addQConstr(cn_expr <= self.new_bp_list[discordant_edges_ccid[bpi]][-1])
			m.addQConstr(cn_expr >= p_bp_cn * self.new_bp_list[discordant_edges_ccid[bpi]][-1])
		for srci in range(lsrc):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + i]
				cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + k + i]
			m.addQConstr(cn_expr <= self.source_edges[src_edges_ccid[srci]][-1])
			
		# Occurrence of breakpoints in each cycle/path
		for i in range(k):
			"""
			for bpi in range(ld):
				if discordant_edges_ccid[bpi] not in self.small_del_indices:
					m.addConstr(x[(lseg + lc + bpi) * k + i] <= max_bp_repeat)
			"""
			for seqi in range(lseg):
                                m.addConstr(x[seqi * k + i] <= max_seq_repeat)

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
			for eni in range(len(endnodes_ccid)):
				cycle_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * eni) * k + i] # (s, v)
			for srci in range(lsrc): 
				cycle_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
			m.addConstr(cycle_expr <= 1.0)

		# special request for c added for max_seq_repeat >= 2
		for i in range(k):
			for node in self.nodes[ccid].keys():
				seqi = self.nodes[ccid][node][0][0]
				m.addConstr(c[k * node_order[node] + i] * x[seq_edges_ccid.index(seqi) * k + i] <= 1.0)
			
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
			for bpi in range(ld):
				dedge = self.new_bp_list[discordant_edges_ccid[bpi]]
				if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
					m.addConstr(y1[(lseg + lc + bpi) * k + i] == 0)
					m.addConstr(y2[(lseg + lc + bpi) * k + i] == 0)
		for i in range(k):
			t_expr_x = gp.LinExpr(0.0)
			t_expr_y = gp.LinExpr(0.0)
			t_expr_yd = gp.QuadExpr(0.0)
			for node in endnodes_ccid:
				t_expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i]
				t_expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i]
				t_expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] * \
						(d[k * nnodes + k + i] - d[k * node_order[node] + i]) # node -> t
				expr_x = gp.LinExpr(0.0)
				expr_y = gp.LinExpr(0.0)
				expr_xc = gp.QuadExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in self.nodes[ccid][node][0]:
					sseg = self.seq_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seq_edges_ccid.index(seqi) * k + i]
					expr_xc += x[seq_edges_ccid.index(seqi) * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seq_edges_ccid.index(seqi) * k + i]
						expr_yd += y1[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[seq_edges_ccid.index(seqi) * k + i]
						expr_yd += y2[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						
				expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] # from s
				expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] * c[k * node_order[node] + i]
				expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] # to t
				expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] * c[k * node_order[node] + i]
				expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] # from s
				expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] * \
								(d[k * node_order[node] + i] - d[k * nnodes + i])
				m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
				m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)
			
			for srci in range(lsrc):
				srce = self.source_edges[src_edges_ccid[srci]]
				t_expr_x += x[(lseg + lc + ld + 2 * srci) * k + k + i]
				t_expr_y += y1[(lseg + lc + ld + 2 * srci) * k + k + i]
				t_expr_yd += y1[(lseg + lc + ld + 2 * srci) * k + k + i] * \
						(d[k * nnodes + k + i] - d[k * node_order[(srce[3], srce[4], srce[5])] + i])
			m.addConstr(t_expr_x * (nnodes + 2) >= d[k * nnodes + k + i])
			m.addConstr(t_expr_y <= 1.0)
			m.addConstr(t_expr_y * nedges * k >= t_expr_x)
			m.addConstr(t_expr_yd >= t_expr_x)
			
			for node in self.nodes[ccid].keys():
				if node not in self.endnodes:
					expr_x = gp.LinExpr(0.0)
					expr_y = gp.LinExpr(0.0)
					expr_xc = gp.QuadExpr(0.0)
					expr_yd = gp.QuadExpr(0.0)
					for seqi in self.nodes[ccid][node][0]:
						sseg = self.seq_edges[seqi]
						node_ = (sseg[0], sseg[1], '-')
						if node_ == node:
							node_ = (sseg[0], sseg[2], '+')
						expr_x += x[seq_edges_ccid.index(seqi) * k + i]
						expr_xc += x[seq_edges_ccid.index(seqi) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[seq_edges_ccid.index(seqi) * k + i]
							expr_yd += y1[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[seq_edges_ccid.index(seqi) * k + i]
							expr_yd += y2[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for ci in self.nodes[ccid][node][1]:
						cedge = self.concordant_edges[ci]
						node_ = (cedge[0], cedge[1], cedge[2])
						if node_ == node:
							node_ = (cedge[3], cedge[4], cedge[5])
						expr_x += x[(lseg + concordant_edges_ccid.index(ci)) * k + i]
						expr_xc += x[(lseg + concordant_edges_ccid.index(ci)) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[(lseg + concordant_edges_ccid.index(ci)) * k + i]
							expr_yd += y1[(lseg + concordant_edges_ccid.index(ci)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[(lseg + concordant_edges_ccid.index(ci)) * k + i]
							expr_yd += y2[(lseg + concordant_edges_ccid.index(ci)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for di in self.nodes[ccid][node][2]:
						dedge = self.new_bp_list[di]
						node_ = (dedge[0], dedge[1], dedge[2])
						if node_ == node:
							node_ = (dedge[3], dedge[4], dedge[5])
						expr_x += x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
						expr_xc += x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
							expr_yd += y1[(lseg + lc + discordant_edges_ccid.index(di)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
							expr_yd += y2[(lseg + lc + discordant_edges_ccid.index(di)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for srci in self.nodes[ccid][node][3]:
						expr_x += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i]
						expr_x += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + k + i]
						expr_xc += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i] * c[k * node_order[node] + i]
						expr_xc += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + k + i] * c[k * node_order[node] + i]
						expr_y += y1[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i]
						expr_yd += y1[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i] * \
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
						m.addConstr(x[seq_edges_ccid.index(edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])
					elif edge[0] == 'c':
						m.addConstr(x[(lseg + concordant_edges_ccid.index(edge[1])) * k + i] >= r[pi * k + i] * path_constraint_[edge])
					else:
						m.addConstr(x[(lseg + lc + discordant_edges_ccid.index(edge[1])) * k + i] >= r[pi * k + i] * path_constraint_[edge])
			
		m.setParam(GRB.Param.Threads, num_threads)
		m.setParam(GRB.Param.NonConvex, 2)
		m.setParam(GRB.Param.TimeLimit, max(7200, ld * 300)) # each breakpoint edge is assigned 5 minutes 
		m.write(of_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "_model.lp")
		m.optimize()
		print('MS:', m.Status)
		if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
			print('Optimization was stopped with status %d' % m.Status)
			return GRB.INFEASIBLE
		else:
			if ccid not in self.cycles:
				self.cycles[ccid] = [[], []] # cycles, paths
				self.cycle_weights[ccid] = [[], []] # cycles, paths
				self.path_constraints_satisfied[ccid] = [[], []] # cycles, paths

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
									cycle[('e', seq_edges_ccid[xi_])] = x_xi
								elif xi_ < lseg + lc:
									cycle[('c', concordant_edges_ccid[xi_ - lseg])] = x_xi
								elif xi_ < lseg + lc + ld:
									cycle[('d', discordant_edges_ccid[xi_ - lseg - lc])] = x_xi
								elif xi_ < lseg + lc + ld + 2 * lsrc:
									assert x_xi == 1
									if (xi_ - lseg - lc - ld) % 2 == 0:
										cycle[('s', src_edges_ccid[(xi_ - lseg - lc - ld) // 2])] = 1 # source edge connected to s
									else:
										cycle[('t', src_edges_ccid[(xi_ - lseg - lc - ld - 1) // 2])] = 1 # source edge connected to t
								else:
									assert x_xi == 1
									if (xi_ - lseg - lc - ld - 2 * lsrc) % 2 == 0:
										nsi = self.endnodes.index(endnodes_ccid[(xi_ - lseg - lc - ld - 2 * lsrc) // 2])
										cycle[('ns', nsi)] = 1 # source edge connected to s
									else:
										nti = self.endnodes.index(endnodes_ccid[(xi_ - lseg - lc - ld - 2 * lsrc - 1) // 2])
										cycle[('nt', nti)] = 1 # source edge connected to t
						for pi in range(len(pc_list)):
							if sol_r[pi * k + i] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
						self.cycles[ccid][1].append(cycle)
						self.cycle_weights[ccid][1].append(sol_w[i])
						self.path_constraints_satisfied[ccid][1].append(path_constraints_s)
					else:
						cycle = dict()
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if xi % k == i and sol_x[xi] >= 0.9:
								xi_ = xi // k
								x_xi = int(round(sol_x[xi]))
								if xi_ < lseg:
									cycle[('e', seq_edges_ccid[xi_])] = x_xi
								elif xi_ < lseg + lc:
									cycle[('c', concordant_edges_ccid[xi_ - lseg])] = x_xi
								elif xi_ < lseg + lc + ld:
									cycle[('d', discordant_edges_ccid[xi_ - lseg - lc])] = x_xi
								else:
									print ("Cyclic path cannot connect to source nodes.")
									os.abort()
						for pi in range(len(pc_list)):
							if sol_r[pi * k + i] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
						self.cycles[ccid][0].append(cycle)
						self.cycle_weights[ccid][0].append(sol_w[i])
						self.path_constraints_satisfied[ccid][0].append(path_constraints_s)
					for seqi in range(lseg):
						total_weights_included += (sol_x[seqi * k + i] * sol_w[i] * self.seq_edges[seq_edges_ccid[seqi]][-2])
			print ("Total weights = ", total_weights_included, total_weights)
			return m.Status
			

	def minimize_cycles_post(self, ccid, total_weights, node_order, pc_list, pc_indices,
				max_seq_repeat = 2, p_total_weight = 0.9, resolution = 0.1, num_threads = 16, of_prefix = ""):
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Cycle decomposition with initial solution from the greedy strategy.")

		k = len(self.cycles[ccid][0]) + len(self.cycles[ccid][1])
		print ("k=", k)
		path_constraint_indices_ = []
		for paths in (self.path_constraints_satisfied[ccid][0] + self.path_constraints_satisfied[ccid][1]):
			for pathi in paths:
				if pathi not in path_constraint_indices_:
					path_constraint_indices_.append(pathi)
		p_path_constraints = 0
		if len(pc_list) > 0:
			p_path_constraints = len(path_constraint_indices_) * 0.9999 / len(pc_indices)

		seq_edges_ccid = [seqi for seqi in range(len(self.seq_edge_ccids)) if self.seq_edge_ccids[seqi] == ccid]
		concordant_edges_ccid = [ci for ci in range(len(self.concordant_edges)) if self.concordant_edge_ccids[ci] == ccid]
		discordant_edges_ccid = [bpi for bpi in range(len(self.new_bp_list)) if self.new_bp_ccids[bpi] == ccid]
		src_edges_ccid = [srci for srci in range(len(self.source_edges)) if self.source_edge_ccids[srci] == ccid]
		endnodes_ccid = [endnode for endnode in self.endnodes if endnode in self.nodes[ccid]]
		lseg, lc, ld, lsrc = len(seq_edges_ccid), len(concordant_edges_ccid), len(discordant_edges_ccid), len(src_edges_ccid)

		nnodes = len(self.nodes[ccid])
		nedges = lseg + lc + ld + 2 * lsrc + 2 * len(endnodes_ccid)

		# Gurobi model
		m = gp.Model(of_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "_cycle_postprocessing_k=" + str(k))
		
		# z[i]: indicating whether cycle or path i exists
		z = m.addVars(k, vtype = GRB.BINARY, name = ["z" + str(i) for i in range(k)])
		
		# w[i]: the weight of cycle or path i, continuous variable
		w = m.addVars(k, lb = 0.0, ub = self.max_CN[ccid], vtype = GRB.CONTINUOUS, name = ["w" + str(i) for i in range(k)])
		
		# Relationship between w[i] and z[i]
		for i in range(k):
			m.addConstr(w[i] <= z[i] * self.max_CN[ccid])

		# x: the number of times an edge occur in cycle or path i
		x_names = []
		for ei in range(nedges):
			for i in range(k):
				x_names.append("x" + str(ei) + "," + str(i))
		x = m.addVars(k * nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)

		# r and R: path constraints
		r = []
		R = []
		if len(pc_list) > 0:
			r_names = []
			R_names = []
			for pi in range(len(pc_list)):
				R_names.append("R" + str(pi))
				for i in range(k):
					r_names.append("r" + str(pi) + "," + str(i))
			r = m.addVars(k * len(pc_list), vtype = GRB.BINARY, name = r_names)
			R = m.addVars(len(pc_list), vtype = GRB.BINARY, name = R_names)

		# Objective: minimize the total number of cycles
		obj = gp.QuadExpr(1.0)
		if len(pc_list) > 0:
                        obj = gp.QuadExpr(2.0)
		for i in range(k):
			obj += z[i]
			for seqi in seq_edges_ccid:
				obj -= (x[seq_edges_ccid.index(seqi) * k + i] * w[i] * self.seq_edges[seqi][-2] / total_weights)
		for pi in range(len(pc_list)):
                        obj -= (R[pi] / len(pc_list))
		m.setObjective(obj, GRB.MINIMIZE)

		# Must include at least 0.9 * total CN weights
		total_weights_expr = gp.QuadExpr(0.0)
		for i in range(k):
			for seqi in seq_edges_ccid:
				total_weights_expr += (x[seq_edges_ccid.index(seqi) * k + i] * w[i] * self.seq_edges[seqi][-2])
		m.addConstr(total_weights_expr >= p_total_weight * total_weights)

		# Eulerian constraint
		for node in self.nodes[ccid].keys():
			if node in self.endnodes:
				for i in range(k):
					m.addConstr(x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] + \
							x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] \
							== x[seq_edges_ccid.index(self.nodes[ccid][node][0][0]) * k + i])
			else:
				for i in range(k):
					ec_expr = gp.LinExpr(0.0)
					for seqi in self.nodes[ccid][node][0]:
						ec_expr += x[seq_edges_ccid.index(seqi) * k + i]
					for ci in self.nodes[ccid][node][1]:
						ec_expr -= x[(lseg + concordant_edges_ccid.index(ci)) * k + i]
					for di in self.nodes[ccid][node][2]:
						# Fix 10/01
						#dedge = self.new_bp_list[di]
						#if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]:
						#	ec_expr -= 2 * x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
						#else:
						ec_expr -= x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
					for srci in self.nodes[ccid][node][3]:
						ec_expr -= x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i] # connected to s
						ec_expr -= x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + k + i] # connected to t
					m.addConstr(ec_expr == 0.0)
		for i in range(k):
			path_expr = gp.LinExpr(0.0)
			for enodei in range(len(endnodes_ccid)):
				path_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + i] # (s, v)
				path_expr -= x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + k + i] # (v, t)
			for srci in range(lsrc): 
				path_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
				path_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # (v, t)
			m.addConstr(path_expr == 0.0)

		# CN constraint
		for seqi in range(lseg):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[seqi * k + i]
			m.addQConstr(cn_expr <= self.seq_edges[seq_edges_ccid[seqi]][-1])
		for ci in range(lc):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + ci) * k + i]
			m.addQConstr(cn_expr <= self.concordant_edges[concordant_edges_ccid[ci]][-1])
		for bpi in range(ld):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + lc + bpi) * k + i]
			m.addQConstr(cn_expr <= self.new_bp_list[discordant_edges_ccid[bpi]][-1])
			#m.addQConstr(cn_expr >= p_bp_cn * self.new_bp_list[discordant_edges_ccid[bpi]][-1])
		for srci in range(lsrc):
			cn_expr = gp.QuadExpr(0.0)
			for i in range(k):
				cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + i]
				cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + k + i]
			m.addQConstr(cn_expr <= self.source_edges[src_edges_ccid[srci]][-1])
			
		# Occurrence of breakpoints in each cycle/path
		for i in range(k):
			"""
			for bpi in range(ld):
				if discordant_edges_ccid[bpi] not in self.small_del_indices:
					m.addConstr(x[(lseg + lc + bpi) * k + i] <= max_bp_repeat)
			"""
			for seqi in range(lseg):
                                m.addConstr(x[seqi * k + i] <= max_seq_repeat)

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
			for eni in range(len(endnodes_ccid)):
				cycle_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * eni) * k + i] # (s, v)
			for srci in range(lsrc): 
				cycle_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
			m.addConstr(cycle_expr <= 1.0)

		# special request for c added for max_seq_repeat >= 2
		for i in range(k):
			for node in self.nodes[ccid].keys():
				seqi = self.nodes[ccid][node][0][0]
				m.addConstr(c[k * node_order[node] + i] * x[seq_edges_ccid.index(seqi) * k + i] <= 1.0)
			
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
			for bpi in range(ld):
				dedge = self.new_bp_list[discordant_edges_ccid[bpi]]
				if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
					m.addConstr(y1[(lseg + lc + bpi) * k + i] == 0)
					m.addConstr(y2[(lseg + lc + bpi) * k + i] == 0)
		for i in range(k):
			t_expr_x = gp.LinExpr(0.0)
			t_expr_y = gp.LinExpr(0.0)
			t_expr_yd = gp.QuadExpr(0.0)
			for node in endnodes_ccid:
				t_expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i]
				t_expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i]
				t_expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] * \
						(d[k * nnodes + k + i] - d[k * node_order[node] + i]) # node -> t
				expr_x = gp.LinExpr(0.0)
				expr_y = gp.LinExpr(0.0)
				expr_xc = gp.QuadExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in self.nodes[ccid][node][0]:
					sseg = self.seq_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seq_edges_ccid.index(seqi) * k + i]
					expr_xc += x[seq_edges_ccid.index(seqi) * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seq_edges_ccid.index(seqi) * k + i]
						expr_yd += y1[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[seq_edges_ccid.index(seqi) * k + i]
						expr_yd += y2[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						
				expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] # from s
				expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] * c[k * node_order[node] + i]
				expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] # to t
				expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + k + i] * c[k * node_order[node] + i]
				expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] # from s
				expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) * k + i] * \
								(d[k * node_order[node] + i] - d[k * nnodes + i])
				m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
				m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)
			
			for srci in range(lsrc):
				srce = self.source_edges[src_edges_ccid[srci]]
				t_expr_x += x[(lseg + lc + ld + 2 * srci) * k + k + i]
				t_expr_y += y1[(lseg + lc + ld + 2 * srci) * k + k + i]
				t_expr_yd += y1[(lseg + lc + ld + 2 * srci) * k + k + i] * \
						(d[k * nnodes + k + i] - d[k * node_order[(srce[3], srce[4], srce[5])] + i])
			m.addConstr(t_expr_x * (nnodes + 2) >= d[k * nnodes + k + i])
			m.addConstr(t_expr_y <= 1.0)
			m.addConstr(t_expr_y * nedges * k >= t_expr_x)
			m.addConstr(t_expr_yd >= t_expr_x)
			
			for node in self.nodes[ccid].keys():
				if node not in self.endnodes:
					expr_x = gp.LinExpr(0.0)
					expr_y = gp.LinExpr(0.0)
					expr_xc = gp.QuadExpr(0.0)
					expr_yd = gp.QuadExpr(0.0)
					for seqi in self.nodes[ccid][node][0]:
						sseg = self.seq_edges[seqi]
						node_ = (sseg[0], sseg[1], '-')
						if node_ == node:
							node_ = (sseg[0], sseg[2], '+')
						expr_x += x[seq_edges_ccid.index(seqi) * k + i]
						expr_xc += x[seq_edges_ccid.index(seqi) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[seq_edges_ccid.index(seqi) * k + i]
							expr_yd += y1[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[seq_edges_ccid.index(seqi) * k + i]
							expr_yd += y2[seq_edges_ccid.index(seqi) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for ci in self.nodes[ccid][node][1]:
						cedge = self.concordant_edges[ci]
						node_ = (cedge[0], cedge[1], cedge[2])
						if node_ == node:
							node_ = (cedge[3], cedge[4], cedge[5])
						expr_x += x[(lseg + concordant_edges_ccid.index(ci)) * k + i]
						expr_xc += x[(lseg + concordant_edges_ccid.index(ci)) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[(lseg + concordant_edges_ccid.index(ci)) * k + i]
							expr_yd += y1[(lseg + concordant_edges_ccid.index(ci)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[(lseg + concordant_edges_ccid.index(ci)) * k + i]
							expr_yd += y2[(lseg + concordant_edges_ccid.index(ci)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for di in self.nodes[ccid][node][2]:
						dedge = self.new_bp_list[di]
						node_ = (dedge[0], dedge[1], dedge[2])
						if node_ == node:
							node_ = (dedge[3], dedge[4], dedge[5])
						expr_x += x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
						expr_xc += x[(lseg + lc + discordant_edges_ccid.index(di)) * k + i] * c[k * node_order[node] + i]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
							expr_yd += y1[(lseg + lc + discordant_edges_ccid.index(di)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						else:
							expr_y += y2[(lseg + lc + discordant_edges_ccid.index(di)) * k + i]
							expr_yd += y2[(lseg + lc + discordant_edges_ccid.index(di)) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					for srci in self.nodes[ccid][node][3]:
						expr_x += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i]
						expr_x += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + k + i]
						expr_xc += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i] * c[k * node_order[node] + i]
						expr_xc += x[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + k + i] * c[k * node_order[node] + i]
						expr_y += y1[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i]
						expr_yd += y1[(lseg + lc + ld + 2 * src_edges_ccid.index(srci)) * k + i] * \
								(d[k * node_order[node] + i] - d[k * nnodes + i])
					m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
					m.addConstr(expr_y <= 1.0)
					m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
					m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)

		# Subpath constraints
		if len(pc_list) > 0:
			sum_R = gp.LinExpr(0.0)
			for pi in range(len(pc_list)):
				sum_R += R[pi]
				path_constraint_ = pc_list[pi]
				sum_r = gp.LinExpr(0.0)
				for ri in range(pi * k, (pi + 1) * k):
					sum_r += r[ri]
					m.addConstr(R[pi] >= r[ri])
				m.addConstr(sum_r >= R[pi])	
				for edge in path_constraint_.keys():
					for i in range(k):
						if edge[0] == 's':
							m.addConstr(x[seq_edges_ccid.index(edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])
						elif edge[0] == 'c':
							m.addConstr(x[(lseg + concordant_edges_ccid.index(edge[1])) * k + i] >= r[pi * k + i] * path_constraint_[edge])
						else:
							m.addConstr(x[(lseg + lc + discordant_edges_ccid.index(edge[1])) * k + i] >= r[pi * k + i] * path_constraint_[edge])
			m.addConstr(sum_R >= p_path_constraints * len(pc_list))

		# Initialize variables
		for i in range(len(self.cycles[ccid][0])):
			z[i].start = 1
			w[i].start = self.cycle_weights[ccid][0][i]
			for (v, vi) in self.cycles[ccid][0][i].keys():
				if v == 'x':
					x[vi * k + i].start = self.cycles[ccid][0][i][(v, vi)]
				elif v == 'c':
					c[vi * k + i].start = self.cycles[ccid][0][i][(v, vi)]
				elif v == 'd':
					d[vi * k + i].start = self.cycles[ccid][0][i][(v, vi)]
				elif v == 'y1':
					y1[vi * k + i].start = self.cycles[ccid][0][i][(v, vi)]
				elif v == 'y2':
					y2[vi * k + i].start = self.cycles[ccid][0][i][(v, vi)]
		for i in range(len(self.cycles[ccid][1])):
			i_ = i + len(self.cycles[ccid][0])
			z[i_].start = 1
			w[i_].start = self.cycle_weights[ccid][1][i]
			for (v, vi) in self.cycles[ccid][1][i].keys():
				if v == 'x':
					x[vi * k + i_].start = self.cycles[ccid][1][i][(v, vi)]
				elif v == 'c':
					c[vi * k + i_].start = self.cycles[ccid][1][i][(v, vi)]
				elif v == 'd':
					d[vi * k + i_].start = self.cycles[ccid][1][i][(v, vi)]
				elif v == 'y1':
					y1[vi * k + i_].start = self.cycles[ccid][1][i][(v, vi)]
				elif v == 'y2':
					y2[vi * k + i_].start = self.cycles[ccid][1][i][(v, vi)]
		for i in range(len(self.path_constraints_satisfied[ccid][0])):
			for pathi in self.path_constraints_satisfied[ccid][0][i]:
				pi = pc_indices.index(pathi)
				r[pi * k + i].start = 1
				R[pi].start = 1
		for i in range(len(self.path_constraints_satisfied[ccid][1])):
			i_ = i + len(self.path_constraints_satisfied[ccid][0])
			for pathi in self.path_constraints_satisfied[ccid][1][i]:
				pi = pc_indices.index(pathi)
				r[pi * k + i_].start = 1
				R[pi].start = 1
		m.update()
		self.cycles[ccid][0].clear()
		self.cycle_weights[ccid][0].clear()
		self.path_constraints_satisfied[ccid][0].clear()
		self.cycles[ccid][1].clear()
		self.cycle_weights[ccid][1].clear()
		self.path_constraints_satisfied[ccid][1].clear()

		m.setParam(GRB.Param.Threads, num_threads)
		m.setParam(GRB.Param.NonConvex, 2)
		m.setParam(GRB.Param.TimeLimit, max(7200, ld * 300)) # each breakpoint edge is assigned 5 minutes 
		m.write(of_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "_postprocessing_model.lp")
		m.optimize()
		print('MS:', m.Status)
		if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
			print('Optimization was stopped with status %d' % m.Status)
			return GRB.INFEASIBLE
		else:
			"""
			if ccid not in self.cycles:
				self.cycles[ccid] = [[], []] # cycles, paths
				self.cycle_weights[ccid] = [[], []] # cycles, paths
				self.path_constraints_satisfied[ccid] = [[], []] # cycles, paths
			"""

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
									cycle[('e', seq_edges_ccid[xi_])] = x_xi
								elif xi_ < lseg + lc:
									cycle[('c', concordant_edges_ccid[xi_ - lseg])] = x_xi
								elif xi_ < lseg + lc + ld:
									cycle[('d', discordant_edges_ccid[xi_ - lseg - lc])] = x_xi
								elif xi_ < lseg + lc + ld + 2 * lsrc:
									assert x_xi == 1
									if (xi_ - lseg - lc - ld) % 2 == 0:
										cycle[('s', src_edges_ccid[(xi_ - lseg - lc - ld) // 2])] = 1 # source edge connected to s
									else:
										cycle[('t', src_edges_ccid[(xi_ - lseg - lc - ld - 1) // 2])] = 1 # source edge connected to t
								else:
									assert x_xi == 1
									if (xi_ - lseg - lc - ld - 2 * lsrc) % 2 == 0:
										nsi = self.endnodes.index(endnodes_ccid[(xi_ - lseg - lc - ld - 2 * lsrc) // 2])
										cycle[('ns', nsi)] = 1 # source edge connected to s
									else:
										nti = self.endnodes.index(endnodes_ccid[(xi_ - lseg - lc - ld - 2 * lsrc - 1) // 2])
										cycle[('nt', nti)] = 1 # source edge connected to t
						for pi in range(len(pc_list)):
							if sol_r[pi * k + i] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
						self.cycles[ccid][1].append(cycle)
						self.cycle_weights[ccid][1].append(sol_w[i])
						self.path_constraints_satisfied[ccid][1].append(path_constraints_s)
					else:
						cycle = dict()
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if xi % k == i and sol_x[xi] >= 0.9:
								xi_ = xi // k
								x_xi = int(round(sol_x[xi]))
								if xi_ < lseg:
									cycle[('e', seq_edges_ccid[xi_])] = x_xi
								elif xi_ < lseg + lc:
									cycle[('c', concordant_edges_ccid[xi_ - lseg])] = x_xi
								elif xi_ < lseg + lc + ld:
									cycle[('d', discordant_edges_ccid[xi_ - lseg - lc])] = x_xi
								else:
									print ("Cyclic path cannot connect to source nodes.")
									os.abort()
						for pi in range(len(pc_list)):
							if sol_r[pi * k + i] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
						self.cycles[ccid][0].append(cycle)
						self.cycle_weights[ccid][0].append(sol_w[i])
						self.path_constraints_satisfied[ccid][0].append(path_constraints_s)
					for seqi in range(lseg):
						total_weights_included += (sol_x[seqi * k + i] * sol_w[i] * self.seq_edges[seq_edges_ccid[seqi]][-2])
			print ("Total weights = ", total_weights_included, total_weights)
			return m.Status


	def maximize_weights_greedy(self, ccid, total_weights, node_order, pc_list, pc_indices, alpha = 0.01,
				max_seq_repeat = 2, p_total_weight = 0.9, resolution = 0.1, num_threads = 16, of_prefix = ""):
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Integer program too large, perform greedy cycle decomposition.")

		seq_edges_ccid = [seqi for seqi in range(len(self.seq_edge_ccids)) if self.seq_edge_ccids[seqi] == ccid]
		concordant_edges_ccid = [ci for ci in range(len(self.concordant_edges)) if self.concordant_edge_ccids[ci] == ccid]
		discordant_edges_ccid = [bpi for bpi in range(len(self.new_bp_list)) if self.new_bp_ccids[bpi] == ccid]
		src_edges_ccid = [srci for srci in range(len(self.source_edges)) if self.source_edge_ccids[srci] == ccid]
		endnodes_ccid = [endnode for endnode in self.endnodes if endnode in self.nodes[ccid]]
		lseg, lc, ld, lsrc = len(seq_edges_ccid), len(concordant_edges_ccid), len(discordant_edges_ccid), len(src_edges_ccid)

		nnodes = len(self.nodes[ccid])
		nedges = lseg + lc + ld + 2 * lsrc + 2 * len(endnodes_ccid)

		remaining_weights = total_weights
		unsatisfied_pc = [pc_indices[i] for i in range(len(pc_indices))]
		remaining_CN = dict()
		for segi in range(lseg):
			remaining_CN[('s', seq_edges_ccid[segi])] = self.seq_edges[seq_edges_ccid[segi]][-1]
		for ci in range(lc):
			remaining_CN[('c', concordant_edges_ccid[ci])] = self.concordant_edges[concordant_edges_ccid[ci]][-1]
		for bpi in range(ld):
			remaining_CN[('d', discordant_edges_ccid[bpi])] = self.new_bp_list[discordant_edges_ccid[bpi]][-1]
		for srci in range(lsrc):
			remaining_CN[('src', src_edges_ccid[srci])] = self.source_edges[src_edges_ccid[srci]][-1]
		#print (remaining_CN)

		next_w = resolution * 1.1
		cycle_id = 0
		num_unsatisfied_pc = len(pc_indices)
		while next_w >= resolution and (remaining_weights > (1.0 - p_total_weight) * total_weights or num_unsatisfied_pc > math.floor(0.1 * len(pc_indices))):
			pp = 1.0
			if alpha > 0 and num_unsatisfied_pc > 0:
				pp = alpha * remaining_weights / num_unsatisfied_pc # multi - objective optimization parameter
			print ("Cycle id = ", cycle_id)
			print ("Remaining weights = ", remaining_weights, total_weights)
			print ("Num unsatisfied path constraints = ", num_unsatisfied_pc, len(pc_indices))
			print ("Path constraints factor = ", pp)

			# Gurobi model
			m = gp.Model(of_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "_cycle_decomposition_greedy_" + str(cycle_id + 1))

			# z[i]: indicating whether cycle or path i exists
			z = m.addVars(1, vtype = GRB.BINARY, name = ["z0"])

			# w[i]: the weight of cycle or path i, continuous variable
			w = m.addVars(1, lb = 0.0, ub = self.max_CN, vtype = GRB.CONTINUOUS, name = ["w0"])

			# Relationship between w[i] and z[i]
			m.addConstr(w[0] <= z[0] * self.max_CN[ccid])
			m.addConstr(w[0] >= z[0] * resolution)

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
				obj += (x[seqi] * w[0] * self.seq_edges[seq_edges_ccid[seqi]][-2])
			for pi in range(len(pc_list)):
				if unsatisfied_pc[pi] >= 0: 
					obj += (r[pi] * max(pp, 1.0))
			m.setObjective(obj, GRB.MAXIMIZE)

			# Eulerian constraint
			for node in self.nodes[ccid].keys():
				if node in self.endnodes:
					m.addConstr(x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node))] + \
							x[(lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)) + 1] \
							== x[seq_edges_ccid.index(self.nodes[ccid][node][0][0])])
				else:
					ec_expr = gp.LinExpr(0.0)
					for seqi in self.nodes[ccid][node][0]:
						ec_expr += x[seq_edges_ccid.index(seqi)]
					for ci in self.nodes[ccid][node][1]:
						ec_expr -= x[lseg + concordant_edges_ccid.index(ci)]
					for di in self.nodes[ccid][node][2]:
						# Fix 10/01
						#dedge = self.new_bp_list[di]
						#if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]:
						#	ec_expr -= 2 * x[lseg + lc + discordant_edges_ccid.index(di)]
						#else:
						ec_expr -= x[lseg + lc + discordant_edges_ccid.index(di)]
					for srci in self.nodes[ccid][node][3]:
						ec_expr -= x[lseg + lc + ld + 2 * src_edges_ccid.index(srci)] # connected to s
						ec_expr -= x[lseg + lc + ld + 2 * src_edges_ccid.index(srci) + 1] # connected to t
					m.addConstr(ec_expr == 0.0)
			path_expr = gp.LinExpr(0.0)
			for enodei in range(len(endnodes_ccid)):
				path_expr += x[lseg + lc + ld + 2 * lsrc + 2 * enodei] # (s, v)
				path_expr -= x[lseg + lc + ld + 2 * lsrc + 2 * enodei + 1] # (v, t)
			for srci in range(lsrc):
				path_expr += x[lseg + lc + ld + 2 * srci] # (s, v)
				path_expr -= x[lseg + lc + ld + 2 * srci + 1] # (v, t)
			m.addConstr(path_expr == 0.0)

			# CN constraint
			for seqi in range(lseg):
				m.addQConstr(w[0] * x[seqi] <= remaining_CN[('s', seq_edges_ccid[seqi])])
			for ci in range(lc):
				m.addQConstr(w[0] * x[lseg + ci] <= remaining_CN[('c', concordant_edges_ccid[ci])])
			for di in range(ld):
				m.addQConstr(w[0] * x[lseg + lc + di] <= remaining_CN[('d', discordant_edges_ccid[di])])
				if self.new_bp_list[discordant_edges_ccid[di]][-1] < resolution:
					m.addConstr(x[lseg + lc + di] == 0.0)
					print ("Set coverage of bp edge at index %d to 0." %(di))
			for srci in range(lsrc):
				cn_expr = gp.QuadExpr(0.0)
				cn_expr += w[0] * x[lseg + lc + ld + 2 * srci]
				cn_expr += w[0] * x[lseg + lc + ld + 2 * srci + 1]
				m.addQConstr(cn_expr <= remaining_CN[('src', src_edges_ccid[srci])])
			
			# Occurrence of breakpoints in each cycle/path
			"""
			for edi in range(ld):
				if discordant_edges_ccid[edi] not in self.small_del_indices:
					m.addConstr(x[lseg + lc + edi] <= max_bp_repeat)
			"""
			for seqi in range(lseg):
				m.addConstr(x[seqi] <= max_seq_repeat)
			
			# c: decomposition i is a cycle, and start at particular node
			c_names = []
			for ni in range(nnodes):
				c_names.append("c" + str(ni))
			c = m.addVars(nnodes, vtype = GRB.BINARY, name = c_names)

			# Relationship between c and x
			cycle_expr = gp.LinExpr(0.0)
			for ni in range(nnodes):
				cycle_expr += c[ni]
			for eni in range(len(endnodes_ccid)):
				cycle_expr += x[lseg + lc + ld + 2 * lsrc + 2 * eni] # (s, v)
			for srci in range(lsrc): 
				cycle_expr += x[lseg + lc + ld + 2 * srci] # (s, v)
			m.addConstr(cycle_expr <= 1.0)
			
			# special request for c added for max_seq_repeat >= 2
			for node in self.nodes[ccid].keys():
				seqi = self.nodes[ccid][node][0][0]
				m.addConstr(c[node_order[node]] * x[seq_edges_ccid.index(seqi)] <= 1.0)

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
			for bpi in range(ld):
				dedge = self.new_bp_list[discordant_edges_ccid[bpi]]
				if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
					m.addConstr(y1[lseg + lc + bpi] == 0)
					m.addConstr(y2[lseg + lc + bpi] == 0)
			t_expr_x = gp.LinExpr(0.0)
			t_expr_y = gp.LinExpr(0.0)
			t_expr_yd = gp.QuadExpr(0.0)
			for node in endnodes_ccid:
				t_expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node) + 1]
				t_expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node) + 1]
				t_expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node) + 1] * \
						(d[nnodes + 1] - d[node_order[node]]) # node -> t
				expr_x = gp.LinExpr(0.0)
				expr_y = gp.LinExpr(0.0)
				expr_xc = gp.QuadExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in self.nodes[ccid][node][0]:
					sseg = self.seq_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seq_edges_ccid.index(seqi)]
					expr_xc += x[seq_edges_ccid.index(seqi)] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seq_edges_ccid.index(seqi)]
						expr_yd += y1[seq_edges_ccid.index(seqi)] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[seq_edges_ccid.index(seqi)]
						expr_yd += y2[seq_edges_ccid.index(seqi)] * (d[node_order[node]] - d[node_order[node_]])
						
				expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)] # from s
				expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)] * c[node_order[node]]
				expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node) + 1] # to t
				expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node) + 1] * c[node_order[node]]
				expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)] # from s
				expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * endnodes_ccid.index(node)] * \
								(d[node_order[node]] - d[nnodes])
				m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges + expr_xc >= expr_x)
				m.addConstr(expr_yd * nedges + expr_xc >= expr_x)
			for srci in range(lsrc):
				srce = self.source_edges[src_edges_ccid[srci]]
				t_expr_x += x[lseg + lc + ld + 2 * srci + 1]
				t_expr_y += y1[lseg + lc + ld + 2 * srci + 1]
				t_expr_yd += y1[lseg + lc + ld + 2 * srci + 1] * \
						(d[nnodes + 1] - d[node_order[(srce[3], srce[4], srce[5])]])
			m.addConstr(t_expr_x * (nnodes + 2) >= d[nnodes + 1])
			m.addConstr(t_expr_y <= 1.0)
			m.addConstr(t_expr_y * nedges >= t_expr_x)
			m.addConstr(t_expr_yd >= t_expr_x)
			
			for node in self.nodes[ccid].keys():
				if node not in self.endnodes:
					expr_x = gp.LinExpr(0.0)
					expr_y = gp.LinExpr(0.0)
					expr_xc = gp.QuadExpr(0.0)
					expr_yd = gp.QuadExpr(0.0)
					for seqi in self.nodes[ccid][node][0]:
						sseg = self.seq_edges[seqi]
						node_ = (sseg[0], sseg[1], '-')
						if node_ == node:
							node_ = (sseg[0], sseg[2], '+')
						expr_x += x[seq_edges_ccid.index(seqi)]
						expr_xc += x[seq_edges_ccid.index(seqi)] * c[node_order[node]]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[seq_edges_ccid.index(seqi)]
							expr_yd += y1[seq_edges_ccid.index(seqi)] * (d[node_order[node]] - d[node_order[node_]])
						else:
							expr_y += y2[seq_edges_ccid.index(seqi)]
							expr_yd += y2[seq_edges_ccid.index(seqi)] * (d[node_order[node]] - d[node_order[node_]])
					for ci in self.nodes[ccid][node][1]:
						cedge = self.concordant_edges[ci]
						node_ = (cedge[0], cedge[1], cedge[2])
						if node_ == node:
							node_ = (cedge[3], cedge[4], cedge[5])
						expr_x += x[lseg + concordant_edges_ccid.index(ci)]
						expr_xc += x[lseg + concordant_edges_ccid.index(ci)] * c[node_order[node]]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[lseg + concordant_edges_ccid.index(ci)]
							expr_yd += y1[lseg + concordant_edges_ccid.index(ci)] * (d[node_order[node]] - d[node_order[node_]])
						else:
							expr_y += y2[lseg + concordant_edges_ccid.index(ci)]
							expr_yd += y2[lseg + concordant_edges_ccid.index(ci)] * (d[node_order[node]] - d[node_order[node_]])
					for di in self.nodes[ccid][node][2]:
						dedge = self.new_bp_list[di]
						node_ = (dedge[0], dedge[1], dedge[2])
						if node_ == node:
							node_ = (dedge[3], dedge[4], dedge[5])
						expr_x += x[lseg + lc + discordant_edges_ccid.index(di)]
						expr_xc += x[lseg + lc + discordant_edges_ccid.index(di)] * c[node_order[node]]
						if node_order[node_] <= node_order[node]:
							expr_y += y1[lseg + lc + discordant_edges_ccid.index(di)]
							expr_yd += y1[lseg + lc + discordant_edges_ccid.index(di)] * (d[node_order[node]] - d[node_order[node_]])
						else:
							expr_y += y2[lseg + lc + discordant_edges_ccid.index(di)]
							expr_yd += y2[lseg + lc + discordant_edges_ccid.index(di)] * (d[node_order[node]] - d[node_order[node_]])
					for srci in self.nodes[ccid][node][3]:
						expr_x += x[lseg + lc + ld + 2 * src_edges_ccid.index(srci)]
						expr_x += x[lseg + lc + ld + 2 * src_edges_ccid.index(srci) + 1]
						expr_xc += x[lseg + lc + ld + 2 * src_edges_ccid.index(srci)] * c[node_order[node]]
						expr_xc += x[lseg + lc + ld + 2 * src_edges_ccid.index(srci) + 1] * c[node_order[node]]
						expr_y += y1[lseg + lc + ld + 2 * src_edges_ccid.index(srci)]
						expr_yd += y1[lseg + lc + ld + 2 * src_edges_ccid.index(srci)] * (d[node_order[node]] - d[nnodes])
					m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
					m.addConstr(expr_y <= 1.0)
					m.addConstr(expr_y * nedges + expr_xc >= expr_x)
					m.addConstr(expr_yd * nedges + expr_xc >= expr_x)

			# Subpath constraints
			for pi in range(len(pc_list)):
				path_constraint_ = pc_list[pi]
				for edge in path_constraint_.keys():
					if edge[0] == 's':
						m.addConstr(x[seq_edges_ccid.index(edge[1])] >= r[pi] * path_constraint_[edge])
					elif edge[0] == 'c':
						m.addConstr(x[lseg + concordant_edges_ccid.index(edge[1])] >= r[pi] * path_constraint_[edge])
					else:
						m.addConstr(x[lseg + lc + discordant_edges_ccid.index(edge[1])] >= r[pi] * path_constraint_[edge])
			
			m.setParam(GRB.Param.Threads, num_threads)
			m.setParam(GRB.Param.NonConvex, 2)
			m.setParam(GRB.Param.TimeLimit, 7200) # each breakpoint edge is assigned 5 minutes
			m.write(of_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "greedy_model_" + str(cycle_id + 1) + "_alpha=" + str(alpha) + ".lp") 
			m.optimize()
			print('MS:', m.Status)
			if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
				print('Optimization was stopped with status %d' % m.Status)
				break
			else:
				if ccid not in self.cycles:
					self.cycles[ccid] = [[], []] # cycles, paths
					self.cycle_weights[ccid] = [[], []] # cycles, paths
					self.path_constraints_satisfied[ccid] = [[], []] # cycles, paths

				cycle_id += 1
				sol_z = m.getAttr('X', z)
				sol_w = m.getAttr('X', w)
				sol_d = m.getAttr('X', d)
				sol_r = m.getAttr('X', r)
				total_weights_included = 0.0
				if sol_z[0] >= 0.9:
					print ("Next cycle/Path exists; weight = %f" %(sol_w[0]))
					next_w = sol_w[0]
					if next_w < resolution:
						break
					sol_x = m.getAttr('X', x)
					sol_y1 = m.getAttr('X', y1)
					sol_y2 = m.getAttr('X', y2)
					sol_c = m.getAttr('X', c)
					cycle_flag = -1
					for ci in range(len(sol_c)):
						if sol_c[ci] >= 0.9:
							cycle_flag = ci
							break
					cycle = dict()
					for xi in range(len(sol_x)):
						cycle[('x', xi)] = sol_x[xi]
					for ci in range(len(sol_c)):
						cycle[('c', ci)] = sol_c[ci]
					for di in range(len(sol_d)):
						cycle[('d', di)] = sol_d[di]
					for yi in range(len(sol_y1)):
						cycle[('y1', yi)] = sol_y1[yi]
					for yi in range(len(sol_y2)):
						cycle[('y2', yi)] = sol_y2[yi]
					if cycle_flag == -1:
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if sol_x[xi] >= 0.9:
								x_xi = int(round(sol_x[xi]))
								if xi < lseg:
									print (cycle_id, 'path', 'seq', x_xi, self.seq_edges[seq_edges_ccid[xi]][:3])
									remaining_CN[('s', seq_edges_ccid[xi])] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('s', seq_edges_ccid[xi])] < resolution:
										remaining_CN[('s', seq_edges_ccid[xi])] = 0.0
								elif xi < lseg + lc:
									print (cycle_id, 'path', 'concordant', x_xi, self.concordant_edges[concordant_edges_ccid[xi - lseg]][:6])
									remaining_CN[('c', concordant_edges_ccid[xi - lseg])] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('c', concordant_edges_ccid[xi - lseg])] < resolution:
										remaining_CN[('c', concordant_edges_ccid[xi - lseg])] = 0.0
								elif xi < lseg + lc + ld:
									print (cycle_id, 'path', 'discordant', x_xi, self.new_bp_list[discordant_edges_ccid[xi - lseg - lc]][:6])
									remaining_CN[('d', discordant_edges_ccid[xi - lseg - lc])] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('d', discordant_edges_ccid[xi - lseg - lc])] < resolution:
										remaining_CN[('d', discordant_edges_ccid[xi - lseg - lc])] = 0.0
								elif xi < lseg + lc + ld + 2 * lsrc:
									assert x_xi == 1
									if (xi - lseg - lc - ld) % 2 == 0:
										#cycle[('s', src_edges_ccid[(xi - lseg - lc - ld) // 2])] = 1 # source edge connected to s
										print (cycle_id, 'path', 'source', x_xi, self.source_edges[src_edges_ccid[(xi - lseg - lc - ld) // 2]][:6])
										remaining_CN[('src', src_edges_ccid[(xi - lseg - lc - ld) // 2])] -= sol_w[0]
										if remaining_CN[('src', src_edges_ccid[(xi - lseg - lc - ld) // 2])] < resolution:
											remaining_CN[('src', src_edges_ccid[(xi - lseg - lc - ld) // 2])] = 0.0
									else:
										#cycle[('t', src_edges_ccid[(xi - lseg - lc - ld - 1) // 2])] = 1 # source edge connected to t
										print (cycle_id, 'path', 'source', x_xi, self.source_edges[src_edges_ccid[(xi - lseg - lc - ld - 1) // 2]][:6])
										remaining_CN[('src', src_edges_ccid[(xi - lseg - lc - ld - 1) // 2])] -= sol_w[0]
										if remaining_CN[('src', src_edges_ccid[(xi - lseg - lc - ld - 1) // 2])] < resolution:
											remaining_CN[('src', src_edges_ccid[(xi - lseg - lc - ld - 1) // 2])] = 0.0
								else:
									assert x_xi == 1
									if (xi - lseg - lc - ld - 2 * lsrc) % 2 == 0:
										nsi = self.endnodes.index(endnodes_ccid[(xi - lseg - lc - ld - 2 * lsrc) // 2])
										#cycle[('ns', nsi)] = 1 # source edge connected to s
										print (cycle_id, 'path', 'source', x_xi, endnodes_ccid[(xi - lseg - lc - ld - 2 * lsrc) // 2])
									else:
										nti = self.endnodes.index(endnodes_ccid[(xi - lseg - lc - ld - 2 * lsrc - 1) // 2])
										#cycle[('nt', nti)] = 1 # source edge connected to t
										print (cycle_id, 'path', 'source', x_xi, endnodes_ccid[(xi - lseg - lc - ld - 2 * lsrc - 1) // 2])
						for pi in range(len(pc_list)):
							if sol_r[pi] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
								unsatisfied_pc[pi] = -1
						self.cycles[ccid][1].append(cycle)
						self.cycle_weights[ccid][1].append(sol_w[0])
						self.path_constraints_satisfied[ccid][1].append(path_constraints_s)
					else:
						path_constraints_s = []
						for xi in range(len(sol_x)):
							if sol_x[xi] >= 0.9:
								x_xi = int(round(sol_x[xi]))
								if xi < lseg:
									print (cycle_id, 'cycle', 'seq', x_xi, self.seq_edges[seq_edges_ccid[xi]][:3])
									remaining_CN[('s', seq_edges_ccid[xi])] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('s', seq_edges_ccid[xi])] < resolution:
										remaining_CN[('s', seq_edges_ccid[xi])] = 0.0
								elif xi < lseg + lc:
									print (cycle_id, 'cycle', 'concordant', x_xi, self.concordant_edges[concordant_edges_ccid[xi - lseg]][:6])
									remaining_CN[('c', concordant_edges_ccid[xi - lseg])] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('c', concordant_edges_ccid[xi - lseg])] < resolution:
										remaining_CN[('c', concordant_edges_ccid[xi - lseg])] = 0.0
								elif xi < lseg + lc + ld:
									print (cycle_id, 'cycle', 'discordant', x_xi, self.new_bp_list[discordant_edges_ccid[xi - lseg - lc]][:6])
									remaining_CN[('d', discordant_edges_ccid[xi - lseg - lc])] -= (sol_x[xi] * sol_w[0])
									if remaining_CN[('d', discordant_edges_ccid[xi - lseg - lc])] < resolution:
										remaining_CN[('d', discordant_edges_ccid[xi - lseg - lc])] = 0.0
								else:
									print ("Cyclic path cannot connect to source nodes.")
									os.abort()
						for pi in range(len(pc_list)):
							if sol_r[pi] >= 0.9:
								path_constraints_s.append(pc_indices[pi])
								unsatisfied_pc[pi] = -1
						self.cycles[ccid][0].append(cycle)
						self.cycle_weights[ccid][0].append(sol_w[0])
						self.path_constraints_satisfied[ccid][0].append(path_constraints_s)
					for seqi in range(lseg):
						total_weights_included += (sol_x[seqi] * sol_w[0] * self.seq_edges[seq_edges_ccid[seqi]][-2])
					print ("Total weights = ", total_weights_included, total_weights)
					remaining_weights -= total_weights_included
					if total_weights_included < 0.005 * total_weights:
						break
				else:
					break
			num_unsatisfied_pc = 0
			for i in range(len(pc_indices)):
				if unsatisfied_pc[i] >= 0:
					num_unsatisfied_pc += 1
		return remaining_weights


	def cycle_decomposition(self, alpha = 0.01, max_seq_repeat = 2, p_total_weight = 0.9, resolution = 0.1, of_prefix = ""):
		for ccid in self.nodes.keys():
			seq_edges_ccid = [seqi for seqi in range(len(self.seq_edge_ccids)) if self.seq_edge_ccids[seqi] == ccid]
			concordant_edges_ccid = [ci for ci in range(len(self.concordant_edges)) if self.concordant_edge_ccids[ci] == ccid]
			discordant_edges_ccid = [bpi for bpi in range(len(self.new_bp_list)) if self.new_bp_ccids[bpi] == ccid]
			src_edges_ccid = [srci for srci in range(len(self.source_edges)) if self.source_edge_ccids[srci] == ccid]
			endnodes_ccid = [endnode for endnode in self.endnodes if endnode in self.nodes[ccid]]
			lseg, lc, ld, lsrc = len(seq_edges_ccid), len(concordant_edges_ccid), len(discordant_edges_ccid), len(src_edges_ccid)

			total_weights = 0.0
			for seqi in seq_edges_ccid:
				sseg = self.seq_edges[seqi]
				total_weights += sseg[-2] * sseg[-1]
			logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Begin cycle decomposition.")
			logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total CN weights = %f." %total_weights)

			path_constraints_ccid = [pathi for pathi in range(len(self.path_constraints[0])) if self.path_constraints[2][pathi] == ccid]
			lpc_ccid = len(path_constraints_ccid)
			self.valid_path_constraints[ccid] = [[], [], []]
			for pathi in path_constraints_ccid:
				path = self.path_constraints[0][pathi]
				path_constraint = dict()
				for ei in range(len(path)):
					if ei % 2 == 0:
						try:
							path_constraint[path[ei]] += 1
						except:
							path_constraint[path[ei]] = 1
				self.valid_path_constraints[ccid][0].append(path_constraint)
				self.valid_path_constraints[ccid][1].append(pathi)
				self.valid_path_constraints[ccid][2].append(self.path_constraints[1][pathi])
			for pathi in range(lpc_ccid)[::-1]:
				path_constraint = self.valid_path_constraints[ccid][0][pathi]
				subpath_flag = -1
				for pathi_ in range(len(self.valid_path_constraints[ccid][0])):
					path_constraint_ = self.valid_path_constraints[ccid][0][pathi_]
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
					del self.valid_path_constraints[ccid][0][pathi]
					del self.valid_path_constraints[ccid][1][pathi]
					self.valid_path_constraints[ccid][2][subpath_flag] = max(self.valid_path_constraints[ccid][2][subpath_flag], self.valid_path_constraints[ccid][2][pathi])
					del self.valid_path_constraints[ccid][2][pathi]

			logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total num subpath constraints = %d." %len(self.valid_path_constraints[ccid][0]))
			for pathi in self.valid_path_constraints[ccid][1]:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - start_time) + "\tSubpath constraints %d = %s" %(pathi, self.path_constraints[0][pathi]))

			k = max(10, ld // 2) # Initial num cycles/paths
			logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total num initial cycles / paths = %d." %k)
			nnodes = len(self.nodes[ccid]) # Does not include s and t
			node_order = dict()
			ni_ = 0
			for node in self.nodes[ccid].keys():
				node_order[node] = ni_
				ni_ += 1
			nedges = lseg + lc + ld + 2 * lsrc + 2 * len(endnodes_ccid)
			if nedges < k:
				k = nedges
				logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Reset num cycles/paths to %d." %k)
			sol_flag = 0
			while k <= nedges:
				if nedges > 100 or (3 * k + 3 * k * nedges + 2 * k * nnodes + k * len(self.valid_path_constraints[ccid][0])) >= 10000:
					remaining_weights = self.maximize_weights_greedy(ccid, total_weights, node_order, self.valid_path_constraints[ccid][0], \
									self.valid_path_constraints[ccid][1], alpha, self.max_seq_repeat[ccid], p_total_weight, resolution, 16, of_prefix)
					print (ccid, self.max_seq_repeat[ccid], 1.0 - remaining_weights / total_weights)
					self.minimize_cycles_post(ccid, total_weights, node_order, self.valid_path_constraints[ccid][0], \
						self.valid_path_constraints[ccid][1], self.max_seq_repeat[ccid], min((1.0 - remaining_weights / total_weights) * 0.9999, 0.8), resolution, 16, of_prefix)
					sol_flag = 1
					break 
				else:
					if self.minimize_cycles(ccid, k, total_weights, node_order, self.valid_path_constraints[ccid][0], self.valid_path_constraints[ccid][1], \
								self.max_seq_repeat[ccid], p_total_weight, 0.9, 16, of_prefix) == GRB.INFEASIBLE:
						k *= 2
					else:
						sol_flag = 1
						break
			if sol_flag == 0:
				remaining_weights = self.maximize_weights_greedy(ccid, total_weights, node_order, self.valid_path_constraints[ccid][0], \
							self.valid_path_constraints[ccid][1], alpha, self.max_seq_repeat[ccid], p_total_weight, resolution, 16, of_prefix)
				self.minimize_cycles_post(ccid, total_weights, node_order, self.valid_path_constraints[ccid][0], \
							self.valid_path_constraints[ccid][1], self.max_seq_repeat[ccid], min((1.0 - remaining_weights / total_weights) * 0.9999, p_total_weight), resolution, 16, of_prefix)



	def eulerian_cycle_t(self, edges_next_cycle, path_constraints_next_cycle):
		lseg = len(self.seq_edges)

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
				ccid = self.seq_edge_ccids[last_seq_edge]
				node = (seq_edge[0], seq_edge[2], '+')
				if last_edge_dir == '-':
					node = (seq_edge[0], seq_edge[1], '-')
				eulerian_cycle.append(node)
				next_bp_edges = [] # Since cycle, only consider discordant edges and concordant edges
				for ci in self.nodes[ccid][node][1]:
					next_bp_edges.append(('c', ci))
				for di in self.nodes[ccid][node][2]:
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
					else:
						bp_edge = self.new_bp_list[next_bp_edges[0][1]][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_cycle.append(node_)
					last_seq_edge = self.nodes[ccid][node_][0][0]
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
					else:
						bp_edge = self.new_bp_list[next_bp_edges[r][1]][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_cycle.append(node_)
					last_seq_edge = self.nodes[ccid][node_][0][0]
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
			#print ('EC', edges_cur)
			path_constraints_cur = path_constraints_next_path.copy()
			src_edge = ()
			last_seq_edge = lseg
			last_edge_dir = '+'
			for edge in edges_cur.keys(): #Start with the edge with smallest index
				if (edge[0] == 's' or edge[0] == 't'):
					src_edge = edge
					node = (self.source_edges[edge[1]][3], self.source_edges[edge[1]][4], self.source_edges[edge[1]][5])
					ccid = self.source_edge_ccids[edge[1]]
					if len(eulerian_path) == 0:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path.append(('$', -1))
						eulerian_path.append(node)
						last_seq_edge = self.nodes[ccid][node][0][0]
					elif self.nodes[ccid][node][0][0] < last_seq_edge:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path[-1] = node
						last_seq_edge = self.nodes[ccid][node][0][0]
				elif (edge[0] == 'ns' or edge[0] == 'nt'):
					src_edge = edge
					node = self.endnodes[edge[1]]
					ccid = -1
					for ccid_ in self.nodes.keys():
						if node in self.nodes[ccid_]:
							ccid = ccid_
							break
					if len(eulerian_path) == 0:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path.append(('$', -1))
						eulerian_path.append(node)
						last_seq_edge = self.nodes[ccid][node][0][0]
					elif self.nodes[ccid][node][0][0] < last_seq_edge:
						last_edge_dir = global_names.neg_plus_minus[node[2]]
						eulerian_path[-1] = node
						last_seq_edge = self.nodes[ccid][node][0][0]					
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
				ccid = self.seq_edge_ccids[last_seq_edge]
				node = (seq_edge[0], seq_edge[2], '+')
				if last_edge_dir == '-':
					node = (seq_edge[0], seq_edge[1], '-')
				eulerian_path.append(node)
				if len(edges_cur) == 1 and (list(edges_cur.keys())[0][0] == 's' or list(edges_cur.keys())[0][0] == 'ns' or \
					list(edges_cur.keys())[0][0] == 't' or list(edges_cur.keys())[0][0] == 'nt'):
					eulerian_path.append(('$', -1))
					break
				next_bp_edges = [] # Since cycle, only consider discordant edges and concordant edges
				for ci in self.nodes[ccid][node][1]:
					next_bp_edges.append(('c', ci))
				for di in self.nodes[ccid][node][2]:
					next_bp_edges.append(('d', di))
				del_list = [i for i in range(len(next_bp_edges)) if next_bp_edges[i] not in edges_cur]
				for i in del_list[::-1]:
					del next_bp_edges[i]
				if len(next_bp_edges) == 0:
					valid = 0
					break
				#print ('EC', edges_cur, next_bp_edges)
				if len(next_bp_edges) == 1: # No branching on the path
					eulerian_path.append(next_bp_edges[0])
					edges_cur[next_bp_edges[0]] = int(edges_cur[next_bp_edges[0]]) - 1
					if edges_cur[next_bp_edges[0]] == 0:
						del edges_cur[next_bp_edges[0]]
					bp_edge = []
					if next_bp_edges[0][0] == 'c':
						bp_edge = self.concordant_edges[next_bp_edges[0][1]][:6]
					else:
						bp_edge = self.new_bp_list[next_bp_edges[0][1]][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_path.append(node_)
					last_seq_edge = self.nodes[ccid][node_][0][0]
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
					else:
						bp_edge = self.new_bp_list[next_bp_edges[r][1]][:6]
					node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
					if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
						node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
					eulerian_path.append(node_)
					last_seq_edge = self.nodes[ccid][node_][0][0]
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


	def output_cycles(self, cycle_file_prefix):
		for ccid in self.nodes.keys():
			fp = open(cycle_file_prefix + "_amplicon" + str(self.ccid2id[ccid]) + "_cycles.txt", 'w')
			interval_num = 1
			ai_ccid = [ai for ai in self.amplicon_intervals if ai[3] == ccid]
			ai_ccid = sorted(ai_ccid, key = lambda ai: (global_names.chr_idx[ai[0]], ai[1]))
			for ai in ai_ccid:
				fp.write("Interval\t%d\t%s\t%d\t%d\n" %(interval_num, ai[0], ai[1], ai[2]))
				interval_num += 1
			fp.write("List of cycle segments\n")
			seq_edges_ccid = [seqi for seqi in range(len(self.seq_edge_ccids)) if self.seq_edge_ccids[seqi] == ccid]
			for segi in range(len(seq_edges_ccid)):
				sseg = self.seq_edges[seq_edges_ccid[segi]]
				fp.write("Segment\t%d\t%s\t%d\t%d\n" %(segi + 1, sseg[0], sseg[1], sseg[2]))
			fp.write("List of subpath constraints\n")
			path_constraint_indices_ = []
			for paths in (self.path_constraints_satisfied[ccid][0] + self.path_constraints_satisfied[ccid][1]):
				for pathi in paths:
					if pathi not in path_constraint_indices_:
						path_constraint_indices_.append(pathi)
			for constraint_i in range(len(self.valid_path_constraints[ccid][1])):
				fp.write("Path constraint\t%d\t" %(constraint_i + 1))
				pathi = self.valid_path_constraints[ccid][1][constraint_i]
				path_ = self.path_constraints[0][pathi]
				assert ccid == self.path_constraints[2][pathi]
				if path_[0][1] > path_[-1][1]:
					path_ = path_[::-1]
				for i in range(len(path_)):
					if i % 4 == 0:
						if i < len(path_) - 1:
							if path_[i + 1][2] == '+':
								fp.write("%d+," %(seq_edges_ccid.index(path_[i][1]) + 1))
							else:
								fp.write("%d-," %(seq_edges_ccid.index(path_[i][1]) + 1))
						else:
							if path_[i - 1][2] == '+':
								fp.write("%d-\t" %(seq_edges_ccid.index(path_[i][1]) + 1))
							else:
								fp.write("%d+\t" %(seq_edges_ccid.index(path_[i][1]) + 1))
				fp.write("Support=%d\t" %(self.valid_path_constraints[ccid][2][constraint_i]))
				if pathi in path_constraint_indices_:
					fp.write("Satisfied\n")
				else:
					fp.write("Unsatisfied\n")
			# sort cycles according to weights
			cycle_indices = sorted([(0, i) for i in range(len(self.cycle_weights[ccid][0]))] + [(1, i) for i in range(len(self.cycle_weights[ccid][1]))], 
						key = lambda item: self.cycle_weights[ccid][item[0]][item[1]], reverse = True)
			#print (ccid, cycle_indices, self.cycle_weights)
			for cycle_i in cycle_indices: 
				cycle_edge_list = []
				if cycle_i[0] == 0: # cycles
					cycle_seg_list = self.eulerian_cycle_t(self.cycles[ccid][cycle_i[0]][cycle_i[1]], self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]])
					print (self.cycles[ccid][cycle_i[0]][cycle_i[1]], cycle_seg_list)
					assert cycle_seg_list[0] == cycle_seg_list[-1]
					fp.write("Cycle=%d;" %(cycle_indices.index(cycle_i) + 1))
					fp.write("Copy_count=%s;" %str(self.cycle_weights[ccid][cycle_i[0]][cycle_i[1]]))
					fp.write("Segments=")
					for segi in range(len(cycle_seg_list) - 2):
						fp.write("%d%s," %(seq_edges_ccid.index(int(cycle_seg_list[segi][:-1]) - 1) + 1, cycle_seg_list[segi][-1]))
					fp.write("%d%s;" %(seq_edges_ccid.index(int(cycle_seg_list[-2][:-1]) - 1) + 1, cycle_seg_list[-2][-1]))
					fp.write("Path_constraints_satisfied=")
					for pathi in range(len(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]]) - 1):
						fp.write("%d," %(self.valid_path_constraints[ccid][1].index(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]][pathi]) + 1))
					if len(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]]) > 0:
						fp.write("%d\n" %(self.valid_path_constraints[ccid][1].index(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]][-1]) + 1))
					else:
						fp.write("\n")
				else: # paths
					#print (self.cycles[ccid][cycle_i[0]][cycle_i[1]])
					cycle_seg_list = self.eulerian_path_t(self.cycles[ccid][cycle_i[0]][cycle_i[1]], self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]])
					fp.write("Cycle=%d;" %(cycle_indices.index(cycle_i) + 1))
					fp.write("Copy_count=%s;" %str(self.cycle_weights[ccid][cycle_i[0]][cycle_i[1]]))
					fp.write("Segments=0+,")
					for segi in range(len(cycle_seg_list) - 1):
						fp.write("%d%s," %(seq_edges_ccid.index(int(cycle_seg_list[segi][:-1]) - 1) + 1, cycle_seg_list[segi][-1]))
					fp.write("%d%s,0-;" %(seq_edges_ccid.index(int(cycle_seg_list[-1][:-1]) - 1) + 1, cycle_seg_list[-1][-1]))
					fp.write("Path_constraints_satisfied=")
					for pathi in range(len(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]]) - 1):
						fp.write("%d," %(self.valid_path_constraints[ccid][1].index(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]][pathi]) + 1))
					if len(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]]) > 0:
						fp.write("%d\n" %(self.valid_path_constraints[ccid][1].index(self.path_constraints_satisfied[ccid][cycle_i[0]][cycle_i[1]][-1]) + 1))
					else:
						fp.write("\n")
			fp.close()


	def closebam(self):
		"""
		Close the short read and long read bam file
		"""
		self.lr_bamfh.close()



if __name__ == '__main__':
	global_names.TSTART = start_time
	parser = argparse.ArgumentParser(description = "Long read only amplicon reconstruction pipeline.")
	parser.add_argument("--lr_bam", help = "Sorted indexed (long read) bam file.", required = True)
	parser.add_argument("--seed", help = "File including seed intervals.", required = True)
	parser.add_argument("--output_prefix", help = "Prefix of output files.", required = True)
	parser.add_argument("--output_bp", help = "If specified, only output the list of breakpoints.",  action = 'store_true')
	parser.add_argument("--cnseg", help = "Long read CNV segmentation file.", required = True)
	parser.add_argument("--ilp_alpha", help = "Parameter used to balance CN weight and path constraints in greedy cycle extraction.", type = float)
	parser.add_argument("--max_repeat", help = "Maximum allowed number of times a sequence edge can be traversed in any cycle/path.", type = int)
	parser.add_argument("--log_fn", help = "Name of log file.")
	args = parser.parse_args()

	log_fn = "infer_breakpoint_graph.log"
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

	b2bn = bam_to_breakpoint_nanopore(args.lr_bam, args.seed)
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Opened LR bam files.")
	b2bn.read_cns(args.cnseg)
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
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed finding all discordant breakpoints.")
	if args.output_bp:
		b2bn.output_breakpoint_info(args.output_prefix + '_breakpoints.tsv')
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Wrote breakpoint information, for all amplicons, to %s." %(args.output_prefix + '_breakpoints.tsv'))
	else:
		b2bn.find_cn_breakpoints()
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed finding breakpoints corresponding to CN changes.")
		b2bn.build_graph()
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Breakpoint graph built for all amplicons.")
		b2bn.assign_cov()
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Fetched read coverage for all sequence and concordant edges.")
		b2bn.assign_cn()
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Computed CN for all edges.")
		b2bn.output_breakpoint_graph(args.output_prefix)
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Wrote breakpoint graph for all complicons to %s." %(args.output_prefix + '_amplicon*_graph.txt'))
		b2bn.compute_path_constraints()
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Computed all subpath constraints.")
		max_repeat = 2
		ilp_alpha = 0.01
		if args.max_repeat:
			max_repeat = args.max_repeat
		if args.ilp_alpha:
			ilp_alpha = args.ilp_alpha
		b2bn.cycle_decomposition(alpha = args.ilp_alpha, max_seq_repeat = max_repeat, of_prefix = args.output_prefix) # This needs to be corrected for parsing max_repeat
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Completed cycle decomposition for all amplicons.")
		b2bn.output_cycles(args.output_prefix)
		logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Wrote cycles for all complicons to %s." %(args.output_prefix + '_amplicon*_cycles.txt'))
		
	b2bn.closebam()
	logging.info("#TIME " + '%.4f\t' %(time.time() - start_time) + "Total runtime.")
	#print (len(b2bn.new_bp_list), b2bn.new_bp_list)



