import os
import sys

import global_names

GAIN = 6.0
CNSIZE_MIN = 99999
CNSIZE_MAX = 5000001  # Not parameterized
CNGAP_MAX = 300001

def interval_overlap(int1, int2):
	"""
	Check if two intervals in the form of [chr, s, e] overlap
	"""
	return (int1[0] == int2[0] and int(int1[1]) <= int(int2[2]) and int(int2[1]) <= int(int1[2]))


def run_seeding(args):

	blocked_intervals = []
	chr_arms = dict()

	if args.gain:
		GAIN = args.gain
	if args.min_seed_size:
		CNSIZE_MIN = args.min_seed_size
	if args.max_seg_gap:
		CNGAP_MAX = args.max_seg_gap

	__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
	with open(os.path.join(__location__, "annotations", "GRCh38_centromere.bed"), 'r') as fp:
		for line in fp:
			s = line.strip().split()
			if len(s[0]) <= 5: # chr1 - chrM
				blocked_intervals.append([s[0], int(s[1]), int(s[2])])
				if 'p' in s[3]:
					chr_arms[s[0]] = [[int(s[1])], [[], []], [int(s[1])]]
				if 'q' in s[3]:
					chr_arms[s[0]][0].append(int(s[2]))
					chr_arms[s[0]][2].append(global_names.chr_sizes[s[0]] - int(s[2]))

	cnv_seeds = []
	cur_seed = []
	with open(args.cn_seg, 'r') as fp:
		for line in fp:
			s = line.strip().split()
			if s[0] != "chromosome":
				if args.cn_seg.endswith(".cns"):
					cn = 2*(2**float(s[4]))
				elif args.cn_seg.endswith(".bed"):
					cn = float(s[3])
				else:
					sys.stderr.write(args.cn_seg + "\n")
					sys.stderr.write("Invalid cn_seg file format!\n")
				# Require absolute CN >= max(GAIN, cn_cutoff_chrarm)
				if cn >= GAIN and (int(s[2]) <= chr_arms[s[0]][0][0] or int(s[1]) >= chr_arms[s[0]][0][1]): 
					# assume input CN segments sorted by chr and pos
					if len(cur_seed) > 0 and s[0] == cur_seed[-1][0] and int(s[1]) - cur_seed[-1][2] <= CNGAP_MAX:
						cur_seed.append((s[0], int(s[1]), int(s[2]), cn))
					else:
						if len(cur_seed) == 0:
							cur_seed = [(s[0], int(s[1]), int(s[2]), cn)]
						else:
							cnv_seeds.append(cur_seed)
							cur_seed = [(s[0], int(s[1]), int(s[2]), cn)]
				if int(s[2]) <= chr_arms[s[0]][0][0]:
					chr_arms[s[0]][1][0].append((s[0], int(s[1]), int(s[2]), cn))
				if int(s[1]) >= chr_arms[s[0]][0][1]:
					chr_arms[s[0]][1][1].append((s[0], int(s[1]), int(s[2]), cn))
	cnv_seeds.append(cur_seed)

	for chr in chr_arms:
		sum_cns_len_parm = sum([cns[2] - cns[1] for cns in chr_arms[chr][1][0]])
		sum_cns_len_qarm = sum([cns[2] - cns[1] for cns in chr_arms[chr][1][1]])
		ccn_p, ccn_q = 2.0, 2.0
		if sum_cns_len_parm >= 0.5 * chr_arms[chr][2][0]:
			scns = sorted(chr_arms[chr][1][0], key=lambda cns: cns[3])
			sum_cns_len_ = 0
			for cns in scns:
				ccn_p = cns[3]
				sum_cns_len_ += (cns[2] - cns[1])
				if sum_cns_len_ >= 0.49 * sum_cns_len_parm:
					break
		if sum_cns_len_qarm >= 0.5 * chr_arms[chr][2][1]:
			scns = sorted(chr_arms[chr][1][1], key=lambda cns: cns[3])
			sum_cns_len_ = 0
			for cns in scns:
				ccn_q = cns[3]
				sum_cns_len_ += (cns[2] - cns[1])
				if sum_cns_len_ >= 0.49 * sum_cns_len_qarm:
					break
		chr_arms[chr].append([ccn_p, ccn_q])
		
	OUTPUT_FN = args.cn_seg.replace(".cns", "_CNV_SEEDS.bed")
	if args.out:
		OUTPUT_FN = args.out
	with open(OUTPUT_FN, 'w') as fp:
		for seed in cnv_seeds:
			sum_seed_len = sum([cns[2] - cns[1] for cns in seed])
			cn_cutoff_chrarm = GAIN
			if sum_seed_len > CNSIZE_MAX:
				cn_cutoff_chrarm = 1.2 * GAIN
			if seed[-1][2] <= chr_arms[seed[-1][0]][0][0]: # p arm
				cn_cutoff_chrarm = cn_cutoff_chrarm + chr_arms[seed[-1][0]][-1][0] - 2.0
			elif seed[0][1] >= chr_arms[seed[-1][0]][0][1]: # q arm
				cn_cutoff_chrarm = cn_cutoff_chrarm + chr_arms[seed[-1][0]][-1][1] - 2.0
			else:
				os.abort()
			for ci in range(len(seed))[::-1]:
				if seed[ci][3] < cn_cutoff_chrarm:
					del seed[ci]
			if len(seed) > 0:
				lastseg = []
				sum_seed_len = 0
				for cns in seed:
					if len(lastseg) > 0 and cns[1] - lastseg[2] <= CNGAP_MAX:
						sum_seed_len += (cns[2] - cns[1])
						lastseg[2] = cns[2]
					else:
						if len(lastseg) == 0:
							lastseg = list(cns)
							sum_seed_len += (cns[2] - cns[1])
						elif sum_seed_len >= CNSIZE_MIN:
							fp.write("%s\t%d\t%d\n" %(lastseg[0], lastseg[1], lastseg[2] - 1))
							sum_seed_len = 0
							lastseg = list(cns)
				if sum_seed_len >= CNSIZE_MIN:
					fp.write("%s\t%d\t%d\n" %(lastseg[0], lastseg[1], lastseg[2] - 1))

	print("Created " + OUTPUT_FN)


