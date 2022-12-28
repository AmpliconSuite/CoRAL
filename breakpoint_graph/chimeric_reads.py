

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


"""
Infer alignment intervals ON THE READ based on the CIGAR string in "SA:Z:" tag 
Require at least one (soft) clip and one match for each canonical alignment record in a chimeric alignment
"""
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
