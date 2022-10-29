import time
import math
import pysam
import argparse
import sys
import os
import numpy as np
import cvxopt
import cvxopt.modeling


def cigar2posSM(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		rs = int(cigar[: cigar.index('S')])
		re = rl - 1
	else:
		rs = 0
		re = rrl - 1
	return rs, re, rrl


def cigar2posMS(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[: cigar.index('M')])
	if strand == '+':
		rs = 0
		re = rrl - 1
	else:
		rs = int(cigar[cigar.index('M') + 1 : cigar.index('S')])
		re = rl - 1
	return rs, re, rrl


def cigar2posSMS(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		rs = int(cigar[: cigar.index('S')])
		re = rs + rrl - 1
	else:
		rs = int(cigar[cigar.index('M') + 1 : -1])
		re = rs + rrl - 1
	return rs, re, rrl


def cigar2posSMD(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[cigar.index('S') + 1 : cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		rs = int(cigar[: cigar.index('S')])
		re = rl - 1
	else:
		rs = 0
		re = int(cigar[cigar.index('S') + 1 : cigar.index('M')]) - 1
	return rs, re, rrl


def cigar2posMDS(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[: cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		rs = 0
		re = int(cigar[: cigar.index('M')]) - 1
	else:
		rs = int(cigar[cigar.index('D') + 1 : cigar.index('S')])
		re = rl - 1
	return rs, re, rrl


def cigar2posSMDS(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[cigar.index('S') + 1: cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		rs = int(cigar[: cigar.index('S')])
		re = rl - int(cigar[cigar.index('D') + 1 : -1]) -1
	else:
		rs = int(cigar[cigar.index('D') + 1 : -1])
		re = rl - int(cigar[: cigar.index('S')]) -1
	return rs, re, rrl


def cigar2posSMI(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		rs = int(cigar[: cigar.index('S')])
		re = rl - 1
	else:
		rs = 0
		re = rl - rrl - 1
	return rs, re, rrl


def cigar2posMIS(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[: cigar.index('M')])
	if strand == '+':
		rs = 0
		re = rl - int(cigar[cigar.index('I') + 1 : cigar.index('S')]) - 1
	else:
		rs = int(cigar[cigar.index('I') + 1 : cigar.index('S')])
		re = rl - 1
	return rs, re, rrl


def cigar2posSMIS(cigar, strand, rl):
	rs, re, rrl = 0, 0, 0
	rrl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		rs = int(cigar[: cigar.index('S')])
		re = rl - int(cigar[cigar.index('I') + 1 : -1]) - 1
	else:
		rs = int(cigar[cigar.index('I') + 1 : -1])
		re = rl - int(cigar[: cigar.index('S')]) -1
	return rs, re, rrl


cigar2pos_ops = {"SM": cigar2posSM,
	"MS": cigar2posMS,
	"SMS": cigar2posSMS,
	"SMD": cigar2posSMD,	
	"MDS": cigar2posMDS,
	"SMDS": cigar2posSMDS,
	"SMI": cigar2posSMI,
	"MIS": cigar2posMIS,
	"SMIS": cigar2posSMIS
}


neg_plus_minus = {"+": "-", "-": "+"}


chr_idx = {'chr1': 0, 'chr2': 1, 'chr3': 2, 'chr4': 3,
		'chr5': 4, 'chr6': 5, 'chr7': 6, 'chr8': 7,
		'chr9': 8, 'chr10': 9, 'chr11': 10, 'chr12': 11,
		'chr13': 12, 'chr14': 13, 'chr15': 14, 'chr16': 15,
		'chr17': 16, 'chr18': 17, 'chr19': 18, 'chr20': 19,
		'chr21': 20, 'chr22': 21, 'chrX': 22, 'chrY': 23, 'chrM': 24}



def interval_overlap(int1, int2):
	return (int1[0] == int2[0] and int(int1[1]) <= int(int2[2]) and int(int2[1]) <= int(int1[2]))


def interval_overlap_l(int1, intl):
	for int2i in range(len(intl)):
		if interval_overlap(int1, intl[int2i]):
			return int2i
	return -1


# A breakpoint (chr1, e1, chr2, s2) must either satisfy chr1 > chr2 or chr1 == chr2 and e1 >= s2 
def bp_match(bp1, bp2, rgap, bp_distance_cutoff):
	if bp1[0] == bp2[0] and bp1[3] == bp2[3] and bp1[2] == bp2[2] and bp1[5] == bp2[5]:
		if rgap <= 0:
			return (abs(int(bp1[1]) - int(bp2[1])) < bp_distance_cutoff and \
				abs(int(bp1[4]) - int(bp2[4])) < bp_distance_cutoff)
		rgap_ = rgap
		consume_rgap = [0, 0]
		if bp1[2] == '+' and int(bp1[1]) <= int(bp2[1]) - bp_distance_cutoff:
			rgap_ -= (int(bp2[1]) - bp_distance_cutoff - int(bp1[1]) + 1)
			consume_rgap[0] = 1
		if bp1[2] == '-' and int(bp1[1]) >= int(bp2[1]) + bp_distance_cutoff:
			rgap_ -= (int(bp1[1]) - int(bp2[1]) - bp_distance_cutoff + 1)
			consume_rgap[0] = 1
		if bp1[5] == '+' and int(bp1[4]) <= int(bp2[4]) - bp_distance_cutoff:
			rgap_ -= (int(bp2[4]) - bp_distance_cutoff - int(bp1[4]) + 1)
			consume_rgap[1] = 1
		if bp1[5] == '-' and int(bp1[4]) >= int(bp2[4]) + bp_distance_cutoff:
			rgap_ -= (int(bp1[4]) - int(bp2[4]) - bp_distance_cutoff + 1)
			consume_rgap[1] = 1
		return (((consume_rgap[0] == 1 and rgap_ >= 0) or (abs(int(bp1[1]) - int(bp2[1])) < bp_distance_cutoff)) and \
			((consume_rgap[1] == 1 and rgap_ >= 0) or (abs(int(bp1[4]) - int(bp2[4])) < bp_distance_cutoff)))
	return False


# Convert split alignment to breakpoint
def interval2bp(R1, R2, rn = '', rgap = 0):
	if (chr_idx[R2[0]] < chr_idx[R1[0]]) or (chr_idx[R2[0]] == chr_idx[R1[0]] and R2[1] < R1[2]):
		return [R1[0], R1[2], R1[3], R2[0], R2[1], neg_plus_minus[R2[3]], rn, rgap, 0]
	return [R2[0], R2[1], neg_plus_minus[R2[3]], R1[0], R1[2], R1[3], rn, rgap, 1]


# Call exact breakpoint from a breakpoint cluster
def bpc2bp(bp_cluster, bp_distance_cutoff):
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
	for bp_ in bp_cluster:
		if bp_match(bp_, bp, bp_[7] * 1.2, bp_distance_cutoff):
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


class bam_to_breakpoint_nanopore():

	sr_bamfh = "" # Short read bam file
	lr_bamfh = "" # Long read bam file

	min_sequence_size = -1
	min_bp_distance_cutoff = 50
	min_clip_cutoff = 1000
	max_interval_cutoff = 2000000
	min_cluster_cutoff = 3
	max_breakpoint_cutoff = 2000
	small_del_cutoff = 10000
	min_bp_len = 10000

	read_length = dict()
	chimeric_alignments = dict()
	large_clip_alignments = dict()
	large_indel_alignments = dict()

	chimeric_alignments_bin = dict()
	
	#median_cov_sr = 0.3658362218228398
	#median_cov_lr = 0.3011307861950682
	amplicon_intervals = []
	amplicon_interval_connections = dict()
	aa_segs_list = []
	"""
	discordant_edges_pos: AA incoming/outgoing discordant edge at a sequence edge junction
	[[], []] = No discordant edges, i.e., CN boundary breakpoint
	[[], [edge indices]] = Only incoming breakpoint edge 
	[[edge indices], []] = Only outgoing breakpoint edge
	[[edge indices], [edge indices]] = Both incoming and outgoing breakpoint edges
	"""
	discordant_edges_pos = dict()
	discordant_edges = []
	small_del_indices = []
	concordant_edges = []
	source_edges = []
	new_bp_list_ = []
	new_bp_stats_ = []
	sr_length = 0.0
	min_sr_alignment_length = 30.0
	normal_cov_sr = 0.3658362218228398
	#normal_cov_sr = 60.70809673628406
	normal_cov_lr = 0.3011307861950682
	#normal_cov_lr = 10.477484296982896
	

	def __init__(self, sr_bamfile, lr_bamfile, aacfile):
		self.sr_bamfh = pysam.AlignmentFile(sr_bamfile, 'rb')
		self.lr_bamfh = pysam.AlignmentFile(lr_bamfile, 'rb')
		with open(aacfile, 'r') as fp:
			for line in fp:
				s = line.strip().split()
				if s[0] == 'Interval':
					self.amplicon_intervals.append([s[2], int(s[3]), int(s[4])])
					if self.min_sequence_size == -1:
						self.min_sequence_size = int(s[4]) - int(s[3])
					else:
						self.min_sequence_size = min(self.min_sequence_size, int(s[4]) - int(s[3]))

	
	def read_graph(self, aagfile):
		lastchr = ''
		lastr = -1
		with open(aagfile, 'r') as fp:
			for line in fp:
				s = line.strip().split()
				if s[0] == 'sequence':
					chr = s[1][: s[1].find(':')]
					l = int(s[1][s[1].find(':') + 1: -1])
					r = int(s[2][s[2].find(':') + 1: -1])
					if chr == lastchr and l == lastr + 1:
						self.discordant_edges_pos[(chr, lastr, l)] = [[], []]
					lastchr = chr
					lastr = r
					self.aa_segs_list.append([chr, l, r, int(s[-1]), -1, 0.0, r - l + 1])
					self.min_sequence_size = min(self.min_sequence_size, r - l)
				if s[0] == 'discordant':
					t = s[1].split('->')
					t0_ = t[0].split(':')
					t1_ = t[1].split(':')
					#print t0_, t1_
					if t0_[1][-1] == '+':
						if (t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1)] = [[], []]
						self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]), int(t0_[1][:-1]) + 1)][0].append(len(self.discordant_edges))
					else:
						if (t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1])) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1]))] = [[], []]
						self.discordant_edges_pos[(t0_[0], int(t0_[1][:-1]) - 1, int(t0_[1][:-1]))][1].append(len(self.discordant_edges))
					if t1_[1][-1] == '+':
						if (t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1)] = [[], []]
						self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]), int(t1_[1][:-1]) + 1)][0].append(len(self.discordant_edges))
					else:
						if (t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1])) not in self.discordant_edges_pos:
							self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1]))] = [[], []]
						self.discordant_edges_pos[(t1_[0], int(t1_[1][:-1]) - 1, int(t1_[1][:-1]))][1].append(len(self.discordant_edges))
					self.discordant_edges.append([t0_[0], int(t0_[1][:-1]), t0_[1][-1], t1_[0], int(t1_[1][:-1]), t1_[1][-1],
									int(s[3]), 0, s[-2], s[-1], set([])])
					if t0_[0] == t1_[0] and t0_[1][-1] == '-' and t1_[1][-1] == '+' and \
						abs(int(t0_[1][:-1]) - int(t1_[1][:-1])) < self.small_del_cutoff:
						#print [t0_[0], int(t0_[1][:-1]), t0_[1][-1], t1_[0], int(t1_[1][:-1]), t1_[1][-1],
						#			int(s[3]), 0, s[-2], s[-1]]
						self.small_del_indices.append(len(self.discordant_edges) - 1)
						self.min_bp_len = min(self.min_bp_len, abs(int(t0_[1][:-1]) - int(t1_[1][:-1])) - 2 * self.min_bp_distance_cutoff)
				if s[0] == 'concordant':
					#print s
					t = s[1].split('->')
					t0_ = t[0].split(':')
					t1_ = t[1].split(':')
					assert t0_[1][-1] == '+' and t1_[1][-1] == '-'
					self.concordant_edges.append([t0_[0], int(t0_[1][:-1]), t0_[1][-1], t1_[0], int(t1_[1][:-1]), t1_[1][-1],
									int(s[3]), -1, s[-2], s[-1]])
				if s[0] == 'source':
					#print s
					t = s[1].split('->')
					t0_ = t[0].split(':')
					t1_ = t[1].split(':')
					#print t0_, t1_
					assert t0_[1][:-1] == '-1'
					self.source_edges.append(['source', -1, '-', t1_[0], int(t1_[1][:-1]), t1_[1][-1], float(s[3]), s[-2], s[-1]])
		#print len(self.aa_segs_list), len(self.concordant_edges)
		#print len(self.discordant_edges_pos)
		self.min_bp_len = max(100, self.min_bp_len)
		#print self.min_bp_len
		#print self.aa_segs_list, self.discordant_edges
		#print self.min_sequence_size
		#print self.discordant_edges_pos
		#self.min_bp_distance_cutoff = min(self.min_bp_distance_cutoff, self.min_sequence_size / 2 - 1)
		#print self.min_bp_distance_cutoff
					


	def fetch(self):
		start_time = time.time()
		for ai in self.amplicon_intervals:
			for read in self.lr_bamfh.fetch(ai[0], ai[1], ai[2] + 1):
				rn = read.query_name
				blocks = read.get_blocks()
				for bi in range(len(blocks) - 1):
					if abs(blocks[bi + 1][0] - blocks[bi][1]) > self.min_bp_len:
						try:
							self.large_indel_alignments[rn].append([ai[0], blocks[bi + 1][0], blocks[bi][1]])
						except:
							self.large_indel_alignments[rn] = [[ai[0], blocks[bi + 1][0], blocks[bi][1]]]
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
			if rn not in self.chimeric_alignments and (rcs[4] > self.min_clip_cutoff or rcs[5] > self.min_clip_cutoff):
				#print read.cigartuples
				self.large_clip_alignments[rn] = (read.cigartuples[0], read.cigartuples[-1])
		print("--- %s seconds ---" % (time.time() - start_time))
		for r in self.chimeric_alignments.keys():
			rl = self.read_length[r]
			r_int = []
			rr_int = []
			q_ = []
			for sa in self.chimeric_alignments[r]:
				t = sa.split(',')
				assert 'S' in t[3] and 'M' in t[3]
				op = ''.join(c for c in t[3] if not c.isdigit())
				rs, re, rrl = cigar2pos_ops[op](t[3], t[2], rl)
				r_int.append([rs, re])
				if t[2] == '+':
					rr_int.append([t[0], int(t[1]) - 1, int(t[1]) + rrl - 2, '+']) # converted to 0 based coordinates
				else:
					rr_int.append([t[0], int(t[1]) + rrl - 2, int(t[1]) - 1, '-']) # converted to 0 based coordinates
				for i in range(int(t[1]) // 10000, (int(t[1]) + rrl) // 10000 + 1):
					try:
						self.chimeric_alignments_bin[t[0]][i].append(r)
					except:
						if t[0] not in self.chimeric_alignments_bin:
							self.chimeric_alignments_bin[t[0]] = dict()
						self.chimeric_alignments_bin[t[0]][i] = [r]
				q_.append(int(t[4]))
			r_int_ind = sorted(range(len(r_int)), key = lambda i: (r_int[i][0], r_int[i][1]))
			r_int = [r_int[i] for i in r_int_ind]
			rr_int = [rr_int[i] for i in r_int_ind]
			q_ = [q_[i] for i in r_int_ind]
			self.chimeric_alignments[r] = [r_int, rr_int, q_]
			"""
			if r == 'SRR12880625.55938':
				print self.chimeric_alignments[r]
			if r == 'SRR12880625.38298':
				print self.chimeric_alignments[r]
			"""
			"""
			if r == 'SRR12880625.73057':
				print self.chimeric_alignments[r]
			if r == 'SRR12880625.4268':
				print self.chimeric_alignments[r]
			"""
		print("--- %s seconds ---" % (time.time() - start_time))


	
	
	def find_breakpoints(self):
		start_time = time.time()
		new_bp_list = []
		for r in self.chimeric_alignments.keys():
			r_int = self.chimeric_alignments[r][0]
			rr_int = self.chimeric_alignments[r][1]
			q_ = self.chimeric_alignments[r][2]
			bassigned = [0 for i in range(len(rr_int) - 1)]
			for i in range(len(rr_int) - 1):
				overlap_i = 0
				for di in range(len(self.discordant_edges)):
					d = self.discordant_edges[di]
					if bp_match([rr_int[i][0], rr_int[i][2], rr_int[i][3], rr_int[i + 1][0], 
							rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]]], d, 
							int(r_int[i + 1][0]) - int(r_int[i][1]), self.min_bp_distance_cutoff):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add(r)
						overlap_i = 1
						bassigned[i] = 1
					elif bp_match([rr_int[i + 1][0], rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]], 
							rr_int[i][0], rr_int[i][2], rr_int[i][3]], d, 
							int(r_int[i + 1][0]) - int(r_int[i][1]), self.min_bp_distance_cutoff):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add(r)
						overlap_i = 1
						bassigned[i] = 1
				if overlap_i == 0 and interval_overlap_l(rr_int[i], self.amplicon_intervals) > 0 and \
					interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) > 0:
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
			for i in range(1, len(rr_int) - 1):
				overlap_i = 0
				if bassigned[i - 1] == 0 and bassigned[i] == 0 and q_[i] < 10 and q_[i - 1] >= 20 and q_[i + 1] >= 20 and \
					interval_overlap_l(rr_int[i - 1], self.amplicon_intervals) > 0 and interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) > 0:
					for di in range(len(self.discordant_edges)):
						d = self.discordant_edges[di]
						if bp_match([rr_int[i - 1][0], rr_int[i - 1][2], rr_int[i - 1][3], rr_int[i + 1][0], 
								rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]]], d, 
								int(r_int[i + 1][0]) - int(r_int[i - 1][1]), self.min_bp_distance_cutoff):
							self.discordant_edges[di][7] += 1
							self.discordant_edges[di][-1].add(r)
							overlap_i = 1
						elif bp_match([rr_int[i + 1][0], rr_int[i + 1][1], neg_plus_minus[rr_int[i + 1][3]], 
								rr_int[i - 1][0], rr_int[i - 1][2], rr_int[i - 1][3]], d, 
								int(r_int[i + 1][0]) - int(r_int[i - 1][1]), self.min_bp_distance_cutoff):
							self.discordant_edges[di][7] += 1
							self.discordant_edges[di][-1].add(r)
							overlap_i = 1
				if overlap_i == 0 and bassigned[i - 1] == 0 and bassigned[i] == 0 and q_[i] < 10 and q_[i - 1] >= 20 and q_[i + 1] >= 20 and \
					interval_overlap_l(rr_int[i - 1], self.amplicon_intervals) > 0 and interval_overlap_l(rr_int[i + 1], self.amplicon_intervals) > 0:
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
				
		#for dei in range(len(self.discordant_edges)):
		#	print (self.discordant_edges[dei])
		new_bp_clusters = []
		new_bp_list.sort(key = lambda item: (chr_idx[item[0]], chr_idx[item[3]], item[1], item[4]))
		#print new_bp_list
		for bpi in range(len(new_bp_list)):
			if bpi == 0:
				new_bp_clusters.append([new_bp_list[bpi]])
				continue
			else:
				bp = new_bp_list[bpi]
				lbp = new_bp_list[bpi - 1]
				if bp[0] == lbp[0] and bp[3] == lbp[3] and bp[2] == lbp[2] and bp[5] == lbp[5] and abs(int(bp[1]) - int(lbp[1])) < 2000 and abs(int(bp[4]) - int(lbp[4])) < 2000:
					new_bp_clusters[-1].append(bp)
				else:
					new_bp_clusters.append([bp])
		print ('New breakpoints ---')
		for c in new_bp_clusters:
			if len(c) >= 3:
				bp, bpr, bp_stats_ = bpc2bp(c, self.min_bp_distance_cutoff)
				if len(set(bpr)) >= 3:
					print (bp, len(set(bpr)))
					if interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals) > 0 and \
						interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals) > 0:
						self.new_bp_list_.append(bp[:-3] + [set(bpr)])
						self.new_bp_stats_.append(bp_stats_)
		print("--- %s seconds ---" % (time.time() - start_time))



	def find_smalldel_breakpoints(self):
		start_time = time.time()
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
							self.min_bp_distance_cutoff):
						self.discordant_edges[di][7] += 1
						self.discordant_edges[di][-1].add(r)
						overlap_ = 1
				if overlap_ == 0:
					new_bp_list.append([rr_gap_[0], rr_gap_[1], '-', rr_gap_[0], rr_gap_[2], '+', r, 0, 0, -1, -1])
		#print new_bp_list
		new_bp_clusters = []
		new_bp_list.sort(key = lambda item: (chr_idx[item[0]], chr_idx[item[3]], item[1], item[4]))
		#print new_bp_list
		for bpi in range(len(new_bp_list)):
			if bpi == 0:
				new_bp_clusters.append([new_bp_list[bpi]])
				continue
			else:
				bp = new_bp_list[bpi]
				lbp = new_bp_list[bpi - 1]
				if bp[0] == lbp[0] and bp[3] == lbp[3] and bp[2] == lbp[2] and bp[5] == lbp[5] and abs(int(bp[1]) - int(lbp[1])) < 2000 and abs(int(bp[4]) - int(lbp[4])) < 2000:
					new_bp_clusters[-1].append(bp)
				else:
					new_bp_clusters.append([bp])
		print ('New breakpoints ---')
		#print new_bp_clusters
		for c in new_bp_clusters:
			if len(c) >= 3:
				bp, bpr, bp_stats_ = bpc2bp(c, self.min_bp_distance_cutoff)
				#print bp, bpr
				if len(set(bpr)) >= 3:
					print (bp, len(set(bpr)))
					if interval_overlap_l([bp[0], bp[1], bp[1]], self.amplicon_intervals) > 0 and \
						interval_overlap_l([bp[3], bp[4], bp[4]], self.amplicon_intervals) > 0:
						self.new_bp_list_.append(bp[:-3] + [set(bpr)])
						self.new_bp_stats_.append(bp_stats_)	
		print("--- %s seconds ---" % (time.time() - start_time))



	def split_seg_bp(self):
		start_time = time.time()
		split_seg = dict()
		for bpi in range(len(self.new_bp_list_)):
			bp = self.new_bp_list_[bpi]
			nsplit = [1, 1]
			for seg in self.aa_segs_list:
				ep, em = -1, -1
				try:
					ep = 1 if len(self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0]) > 0 else 0
					em = 1 if len(self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1]) > 0 else 0
				except:
					pass
				if bp[0] == seg[0]:
					if bp[2] == '+':
						if seg[2] - self.min_bp_distance_cutoff < bp[1] <= seg[2]:
							self.new_bp_list_[bpi][1] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
						elif seg[2] < bp[1] < seg[2] + self.min_bp_distance_cutoff and em == 0:
							self.new_bp_list_[bpi][1] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
					if bp[2] == '-':
						if seg[2] < bp[1] <= seg[2] + self.min_bp_distance_cutoff and ep + em >= 0:
							self.new_bp_list_[bpi][1] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
						elif seg[2] - self.min_bp_distance_cutoff < bp[1] < seg[2] and ep == 0:
							self.new_bp_list_[bpi][1] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[0] = 0
				if bp[3] == seg[0]:
					if bp[5] == '+':
						if seg[2] - self.min_bp_distance_cutoff < bp[4] <= seg[2]:
							self.new_bp_list_[bpi][4] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
						elif seg[2] < bp[4] < seg[2] + self.min_bp_distance_cutoff and em == 0:
							self.new_bp_list_[bpi][4] = seg[2]
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][0].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
					if bp[5] == '-':
						if seg[2] < bp[4] <= seg[2] + self.min_bp_distance_cutoff and ep + em >= 0:
							self.new_bp_list_[bpi][4] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
						elif seg[2] - self.min_bp_distance_cutoff < bp[4] < seg[2] and ep == 0:
							self.new_bp_list_[bpi][4] = seg[2] + 1
							self.discordant_edges_pos[(seg[0], int(seg[2]), int(seg[2]) + 1)][1].append(len(self.discordant_edges) + bpi)
							nsplit[1] = 0
			#print bp, nsplit
			for segi in range(len(self.aa_segs_list)):	
				seg = self.aa_segs_list[segi]
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
		#print split_seg
		new_segs = dict()
		for segi in split_seg:
			#print split_seg[segi]
			split_seg[segi].sort(key = lambda item: item[0])
			sseg = self.aa_segs_list[segi]	
			new_segs[segi] = []
			lssi = 0
			del_list = []
			for ssi in range(len(split_seg[segi])):
				spp = split_seg[segi][ssi]
				if ssi > 0 and int(spp[0]) - int(split_seg[segi][lssi][1]) <= self.min_bp_distance_cutoff and \
					spp[3] == split_seg[segi][lssi][3]:
					if self.new_bp_list_[spp[2]][spp[3] + 1] == '+':
						self.new_bp_list_[spp[2]][spp[3]] = split_seg[segi][lssi][0]
					else:
						self.new_bp_list_[spp[2]][spp[3]] = split_seg[segi][lssi][1]
					del_list.append(ssi)
				else:
					lssi = ssi
			#print del_list[::-1]
			for dssi in del_list[::-1]:
				del split_seg[segi][dssi]
			#print split_seg[segi]
			for ssi in range(len(split_seg[segi])):
				if ssi == 0:
					new_segs[segi].append([sseg[0], sseg[1], split_seg[segi][ssi][0], -1, -1, 0.0, 
								int(split_seg[segi][ssi][0]) - int(sseg[1]) + 1])
				else:
					new_segs[segi].append([sseg[0], split_seg[segi][ssi - 1][1], split_seg[segi][ssi][0], 
								-1, -1, 0.0, int(split_seg[segi][ssi][0]) - \
								int(split_seg[segi][ssi - 1][1]) + 1])
			new_segs[segi].append([sseg[0], split_seg[segi][-1][1], sseg[2], -1, -1, 0.0, 
					int(sseg[2]) - int(split_seg[segi][-1][1]) + 1])
		for segi in sorted(split_seg.keys())[::-1]:
			del self.aa_segs_list[segi]
			for nsi in range(len(new_segs[segi])):
				self.aa_segs_list.insert(segi + nsi, new_segs[segi][nsi]) # Add sequence edges
			for nsi in range(len(new_segs[segi]) - 1):
				self.concordant_edges.append([new_segs[segi][nsi][0], new_segs[segi][nsi][2], '+', 
					new_segs[segi][nsi + 1][0], new_segs[segi][nsi + 1][1], '-', -1, -1, 'None', 'None'])
		print("--- %s seconds ---" % (time.time() - start_time))



	def assign_cov_sequence(self):
		start_time = time.time()
		max_seg = -1
		avg_rl = 0.0
		for seg in self.aa_segs_list:
			#print seg
			if seg[3] == -1:
				rl_list = [read.infer_read_length() for read in self.sr_bamfh.fetch(seg[0], seg[1], seg[2] + 1) if read.infer_read_length()]
				#seg[3] = sum([sum(nc) for nc in self.sr_bamfh.fetch(seg[0], seg[1], seg[2] + 1)]) \
				#		/ max(1.0, float(seg[2] - seg[1] + 1))
				len_rl_list = len(rl_list)
				seg[3] = len_rl_list
				if len_rl_list > max_seg:
					max_seg = len_rl_list
					avg_rl = np.average(rl_list)
			if seg[4] == -1:
				#seg[4] = sum([sum(nc) for nc in self.lr_bamfh.count_coverage(seg[0], seg[1], seg[2] + 1)]) \
				#		/ max(1.0, float(seg[2] - seg[1] + 1))
				rl_list = [read.infer_read_length() for read in self.lr_bamfh.fetch(seg[0], seg[1], seg[2] + 1)]
				#print (rl_list)
				seg[4] = len(rl_list)
				seg[5] = np.average(rl_list)
		if self.sr_length == 0.0:
			print (avg_rl)
			self.sr_length = avg_rl
		print("assign_cov --- %s seconds ---" % (time.time() - start_time))


	
	def output_breakpoint_graph(self, ogfile):
		with open(ogfile, 'w') as fp:
			fp.write("SequenceEdge: StartPosition, EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads, Size\n")
			for segi in range(len(self.aa_segs_list)):
				sseg = self.aa_segs_list[segi]
				fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%s\t%s\t%s\n" %(sseg[0], sseg[1], sseg[0], sseg[2], sseg[-1], sseg[3], sseg[4], str(sseg[6])))
			fp.write("BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads, ")
			fp.write("HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence\n")
			for es in self.source_edges:
				fp.write("source\t%s:%s%s->%s:%s%s\t%f\t%d\t-1\t%s\t%s\n" %(es[0], es[1], es[2], es[3], 
					es[4], es[5], es[-1], int(es[6]), es[7], es[8]))
			for ec in self.concordant_edges:
				#print ec
				#print([type(item) for item in ec])
				rls = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[0], start = ec[1], stop = ec[1] + 1)])
				rrs = set([read.query_name for read in self.lr_bamfh.fetch(contig = ec[3], start = ec[4], stop = ec[4] + 1)])
				rbps = set([])
				for bpi in self.discordant_edges_pos[(ec[0], ec[1], ec[4])][0]:
					if bpi >= len(self.discordant_edges):
						rbps |= self.new_bp_list_[bpi - len(self.discordant_edges)][-2]
					else:
						rbps |= self.discordant_edges[bpi][-2]
				for bpi in self.discordant_edges_pos[(ec[0], ec[1], ec[4])][1]:
					if bpi >= len(self.discordant_edges):
						rbps |= self.new_bp_list_[bpi - len(self.discordant_edges)][-2]
					else:
						rbps |= self.discordant_edges[bpi][-2]
				if ec[6] == -1:
					rls_ = set([read.query_name for read in self.sr_bamfh.fetch(contig = ec[0], start = ec[1], stop = ec[1] + 1)])
					rrs_ = set([read.query_name for read in self.sr_bamfh.fetch(contig = ec[3], start = ec[4], stop = ec[4] + 1)])
					fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\t%s\t%s\n" %(ec[0], ec[1], ec[2], ec[3], 
						ec[4], ec[5], ec[-1], len(rls_ & rrs_), len((rls & rrs) - rbps), ec[8], ec[9]))
				else:
					fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\t%s\t%s\n" %(ec[0], ec[1], ec[2], ec[3], 
						ec[4], ec[5], ec[-1], ec[6], len((rls & rrs) - rbps), ec[8], ec[9]))
				#print ec, rls, rrs, len(rls & rrs)
			for ed in self.discordant_edges:
				fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\t%s\t%s\n" %(ed[0], ed[1], ed[2], ed[3], 
					ed[4], ed[5], ed[-1], ed[6], ed[7], ed[8], ed[9]))
			for bp in self.new_bp_list_:
				fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t-1\t%d\tNone\tNone\n" %(bp[0], bp[1], bp[2], bp[3], bp[4], bp[5], bp[-1], len(bp[6])))
	

	def output_breakpoint_info(self, obpfile):
		with open(obpfile, 'w') as fp:
			fp.write("chr1\tpos1\tchr2\tpos2\torientation\tsr_support\tlr_support\tlr_info=[avg1, avg2, std1, std2, mapq1, mapq2]\n")
			for ed in self.discordant_edges:
				fp.write("%s\t%s\t%s\t%s\t%s%s\t%d\t%d\tN/A\n" %(ed[0], ed[1], ed[3], ed[4], ed[2], ed[5], 
					ed[6], ed[7]))
			for bpi in range(len(self.new_bp_list_)):
				bp = self.new_bp_list_[bpi]
				bp_stats = self.new_bp_stats_[bpi]
				fp.write("%s\t%s\t%s\t%s\t%s%s\t-1\t%d\t%s\n" 
					%(bp[0], bp[1], bp[3], bp[4], bp[2], bp[5], len(bp[-1]), bp_stats))
	

	def closebam(self):
		self.sr_bamfh.close()
		self.lr_bamfh.close()


	def assign_cn(self, mode = "both"):
		lseg = len(self.aa_segs_list)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lnbp = len(self.new_bp_list_)
		lsrc = len(self.source_edges)

		# Prepare adjacency list
		endnodes = []
		nodes = dict()
		for segi in range(lseg):
			sseg = self.aa_segs_list[segi]
			nodes[(sseg[0], sseg[1])] = [[segi], [], [], []] # sequence edges, breakpoint edges
			nodes[(sseg[0], sseg[2])] = [[segi], [], [], []]
			if segi == 0:
				endnodes.append((sseg[0], sseg[1]))
			else:
				lastseg = self.aa_segs_list[segi - 1]
				if lastseg[0] != sseg[0] or lastseg[2] + 1 != sseg[1]:
					endnodes.append((lastseg[0], lastseg[2]))
					endnodes.append((sseg[0], sseg[1]))
				if segi == lseg - 1:
					endnodes.append((sseg[0], sseg[2]))
		for eci in range(lc):
			ec = self.concordant_edges[eci]
			nodes[(ec[0], ec[1])][1].append(eci)
			nodes[(ec[3], ec[4])][1].append(eci)
		for edi in range(ld):
			ed = self.discordant_edges[edi]
			nodes[(ed[0], ed[1])][2].append(edi)
			nodes[(ed[3], ed[4])][2].append(edi)
		for bpi in range(lnbp):
			bp = self.new_bp_list_[bpi]
			#print (bp)
			nodes[(bp[0], bp[1])][2].append(ld + bpi)
			nodes[(bp[3], bp[4])][2].append(ld + bpi)
		for esi in range(lsrc)[::-1]: # Clean up AA source edges
			es = self.source_edges[esi]
			if (es[3], es[4]) in endnodes:
				del self.source_edges[esi]
		lsrc = len(self.source_edges)
		for esi in range(lsrc):
			es = self.source_edges[esi]
			nodes[(es[3], es[4])][3].append(esi)

		nconstraints = len([node for node in nodes.keys() if node not in endnodes])
		nvariables = lseg + lc + ld + lnbp + lsrc
		wcn = []
		wlncn = []
		if mode == "both":
			wcn = [(self.normal_cov_sr * sseg[-1] / self.sr_length) + \
				(self.normal_cov_lr * sseg[-1] / sseg[5]) for sseg in self.aa_segs_list]
			wcn += [self.normal_cov_sr * (self.sr_length - 1) / self.sr_length + self.normal_cov_lr \
				for eci in range(lc)]
			wcn += [self.normal_cov_sr * (self.sr_length - 2 * self.min_sr_alignment_length) / \
				self.sr_length + self.normal_cov_lr for edi in range(ld)]
			wcn += [self.normal_cov_lr for bpi in range(lnbp)]
			wcn += [self.normal_cov_sr * (self.sr_length - 2 * self.min_sr_alignment_length) / self.sr_length \
				for es in self.source_edges]
			wlncn = [sseg[3] + sseg[4] for sseg in self.aa_segs_list]
			wlncn += [ec[6] + ec[7] for ec in self.concordant_edges]
			wlncn += [ed[6] + ed[7] for ed in self.discordant_edges]
			wlncn += [len(bp[6]) for bp in self.new_bp_list_]
			wlncn += [es[6] if es[6] >= 1 else 0.1 for es in self.source_edges]
			wcn = cvxopt.matrix(wcn)
			wlncn = cvxopt.matrix(wlncn)
		elif mode == "sr_only":
			wcn = [(self.normal_cov_sr * sseg[-1] / self.sr_length) for sseg in self.aa_segs_list]
			wcn += [self.normal_cov_sr * (self.sr_length - 1) / self.sr_length \
				for eci in range(lc)]
			wcn += [self.normal_cov_sr * (self.sr_length - 2 * self.min_sr_alignment_length) / \
				self.sr_length for edi in range(ld)]
			wcn += [0.1 for bpi in range(lnbp)]
			wcn += [self.normal_cov_sr * (self.sr_length - 2 * self.min_sr_alignment_length) / self.sr_length \
				for es in self.source_edges]
			wlncn = [sseg[3] for sseg in self.aa_segs_list]
			wlncn += [ec[6] for ec in self.concordant_edges]
			wlncn += [ed[6] for ed in self.discordant_edges]
			wlncn += [0.1 for bp in self.new_bp_list_]
			wlncn += [es[6] if es[6] >= 1 else 0.1 for es in self.source_edges]
			wcn = cvxopt.matrix(wcn)
			wlncn = cvxopt.matrix(wlncn)
		else:
			wcn = [(self.normal_cov_lr * sseg[-1] / sseg[5]) for sseg in self.aa_segs_list]
			wcn += [self.normal_cov_lr for eci in range(lc)]
			wcn += [self.normal_cov_lr for edi in range(ld)]
			wcn += [self.normal_cov_lr for bpi in range(lnbp)]
			wcn += [0.1 for es in self.source_edges]
			wlncn = [sseg[4] for sseg in self.aa_segs_list]
			wlncn += [ec[7] for ec in self.concordant_edges]
			wlncn += [ed[7] for ed in self.discordant_edges]
			wlncn += [len(bp[6]) for bp in self.new_bp_list_]
			wlncn += [0.1 for es in self.source_edges]
			wcn = cvxopt.matrix(wcn)
			wlncn = cvxopt.matrix(wlncn)
		ci = 0
		balance_constraints = np.zeros([nconstraints, nvariables])
		for node in nodes.keys():
			if node not in endnodes:
				for segi in nodes[node][0]:
					balance_constraints[ci][segi] = 1
				for eci in nodes[node][1]:
					balance_constraints[ci][lseg + eci] = -1
				for edi in nodes[node][2]:
					balance_constraints[ci][lseg + lc + edi] = -1
				for esi in nodes[node][3]:
					balance_constraints[ci][lseg + lc + ld + lnbp + esi] = -1
				ci += 1
		print (np.linalg.matrix_rank(balance_constraints))
		balance_constraints = cvxopt.matrix(balance_constraints)
		# Convex optimization function required by cvxopt
		def F(x = None, z = None):
			if x is None: 
				return 0, cvxopt.matrix(1.0, (nvariables, 1))
			if min(x) <= 0.0001: 
				return None
			f = cvxopt.modeling.dot(wcn, x) - cvxopt.modeling.dot(wlncn, cvxopt.log(x))
			Df = (wcn - cvxopt.mul(wlncn, (x**-1))).T
			if z is None: 
				return f, Df
			#print ((z[0] * cvxopt.mul(wlncn, x**-2)).size)
			#print (z[0])
			H = cvxopt.spdiag(z[0] * cvxopt.mul(wlncn, x**-2))
			#print ('x', x)
			#print('R', np.linalg.matrix_rank(H))
			#print ('H', H)
			return f, Df, H
		print(wcn.size, wlncn.size, balance_constraints.size)
		options = {'maxiters': 1000}
		sol = cvxopt.solvers.cp(F, A = balance_constraints, b = cvxopt.matrix([0.0 for i in range(nconstraints)]), kktsolver = 'ldl', options = options)
		if sol['status'] == "optimal" or sol['status'] == "unknown":
			for segi in range(lseg):
				self.aa_segs_list[segi] += [sol['x'][segi]]
			for eci in range(lc):
				self.concordant_edges[eci] += [sol['x'][lseg + eci]]
			for edi in range(ld):
				self.discordant_edges[edi] += [sol['x'][lseg + lc + edi]]
			for bpi in range(lnbp):
				self.new_bp_list_[bpi] += [sol['x'][lseg + lc + ld + bpi]]
			for esi in range(len(self.source_edges)):
				self.source_edges[esi] += [sol['x'][lseg + lc + ld + lnbp + esi]]
	# colo320 normal cov: 0.3011307861950682, rl = 
	# gbm39 normal cov: 10.419352080766544


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Examine ")
	parser.add_argument("--sr_bam", help = "Sorted indexed bam file.", required = True)
	parser.add_argument("--lr_bam", help = "Sorted indexed bam file.", required = True)
	parser.add_argument("--aa_graph", help = "AA-formatted graph file.", required = True)
	parser.add_argument("--aa_cycle", help = "AA-formatted cycle file.", required = True)
	args = parser.parse_args()

	b2bn = bam_to_breakpoint_nanopore(args.sr_bam, args.lr_bam, args.aa_cycle)
	b2bn.read_graph(args.aa_graph)
	
	b2bn.fetch()
	b2bn.find_smalldel_breakpoints()
	b2bn.find_breakpoints()
	b2bn.split_seg_bp()
	#b2bn.assign_cov_sequence()
	#b2bn.assign_cn()

	#b2bn.output_breakpoint_graph(args.aa_graph.split('/')[-1][:-4] + '_.txt')
	b2bn.output_breakpoint_info(args.aa_graph.split('/')[-1][:-9] + 'breakpoints.txt')
	b2bn.closebam()

