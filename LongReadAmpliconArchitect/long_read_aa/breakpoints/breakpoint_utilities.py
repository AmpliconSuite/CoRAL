"""
Utilities for breakpoint graph inference.
"""
import math
import numpy as np
from collections import Counter

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

