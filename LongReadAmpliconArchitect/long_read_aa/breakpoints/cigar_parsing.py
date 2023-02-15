import logging
import global_names


def cigar2posSM(cigar, strand, rl):
	"""
	convert cigar string in format *S*M into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = rl - 1
	else:
		qs = 0
		qe = rl - 1
	return qs, qe, rl


def cigar2posMS(cigar, strand, rl):
	"""
	convert cigar string in format *M*S into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[: cigar.index('M')])
	if strand == '+':
		qs = 0
		qe = rl - 1
	else:
		qs = int(cigar[cigar.index('M') + 1 : cigar.index('S')])
		qe = rl - 1
	return qs, qe, rl


def cigar2posSMS(cigar, strand, rl):
	"""
	convert cigar string in format *S*M*S into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = qs + rl - 1
	else:
		qs = int(cigar[cigar.index('M') + 1 : -1])
		qe = qs + rl - 1
	return qs, qe, rl


def cigar2posSMD(cigar, strand, rl):
	"""
	convert cigar string in format *S*M*D into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[cigar.index('S') + 1 : cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = rl - 1
	else:
		qs = 0
		qe = int(cigar[cigar.index('S') + 1 : cigar.index('M')]) - 1
	return qs, qe, rl


def cigar2posMDS(cigar, strand, rl):
	"""
	convert cigar string in format *M*D*S into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[: cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		qs = 0
		qe = int(cigar[: cigar.index('M')]) - 1
	else:
		qs = int(cigar[cigar.index('D') + 1 : cigar.index('S')])
		qe = rl - 1
	return qs, qe, rl


def cigar2posSMDS(cigar, strand, rl):
	"""
	convert cigar string in format *S*M*D*S into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[cigar.index('S') + 1: cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = rl - int(cigar[cigar.index('D') + 1 : -1]) -1
	else:
		qs = int(cigar[cigar.index('D') + 1 : -1])
		qe = rl - int(cigar[: cigar.index('S')]) -1
	return qs, qe, rl


def cigar2posSMI(cigar, strand, rl):
	"""
	convert cigar string in format *S*M*I into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = rl - 1
	else:
		qs = 0
		qe = rl - rl - 1
	return qs, qe, rl


def cigar2posMIS(cigar, strand, rl):
	"""
	convert cigar string in format *M*I*S into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[: cigar.index('M')])
	if strand == '+':
		qs = 0
		qe = rl - int(cigar[cigar.index('I') + 1 : cigar.index('S')]) - 1
	else:
		qs = int(cigar[cigar.index('I') + 1 : cigar.index('S')])
		qe = rl - 1
	return qs, qe, rl


def cigar2posSMIS(cigar, strand, rl):
	"""
	convert cigar string in format *S*M*I*S into alignment
	return the start and end position on the read, and the alignment length on the reference genome 
	"""
	qs, qe, rl = 0, 0, 0
	rl = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = rl - int(cigar[cigar.index('I') + 1 : -1]) - 1
	else:
		qs = int(cigar[cigar.index('I') + 1 : -1])
		qe = rl - int(cigar[: cigar.index('S')]) -1
	return qs, qe, rl


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


"""
Chimeric alignments of long reads can be represented by the following 3 lists
qint: 0-based intervals on query (reads) for each local alignment
rint: 0-based intervals on reference genome for each local alignment
qual: mapping qualities for each local alignment
"""
def alignment_from_blocks(chr, blocks, cutoff, mapq = None):
	"""
	Convert blocks from pysam get_blocks into a new chimeric alignment.
	Merge blocks within distance cutoff.
	Note that the blocks must come from a single alignment record.
	Skip intervals on read and only fill in intervals on the reference. 
	If a single mapping quality value is given, set qual to a list of this value. 
	"""
	qint, rint, qual = [], [], []
	bl = blocks[0][0]
	for bi in range(len(blocks) - 1):
		if abs(blocks[bi + 1][0] - blocks[bi][1]) >= cutoff:
			if bl > blocks[bi][1]:
				rint.append([chr, bl, blocks[bi][1], '-'])
			else:
				rint.append([chr, bl, blocks[bi][1], '+'])
			bl = blocks[bi + 1][0]
	if bl > blocks[-1][1]:
		rint.append([chr, bl, blocks[-1][1], '-'])
	else:
		rint.append([chr, bl, blocks[-1][1], '+'])
	if mapq:
		qual = [mapq for i_ in self.rint]
	return (qint, rint, qual)


def alignment_from_satags(sa_list, rl):
	"""
	Convert "SA:Z" a list of strings into a new chimeric alignment.
	"""
	qint, rint, qual = [], [], []
	for sa in sa_list:
		t = sa.split(',')
		if 'S' not in t[3] or 'M' not in t[3]:
			"""
			Require a chimeric alignment record having at least some (soft)clips and matches 
			"""
			logging.warning("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Found chimeric alignment without match or soft clips.")
			logging.warning("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tRead name: %s; Read length: %d." %(r, rl))
			logging.warning("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tAll CIGAR strings: %s." %(self.chimeric_alignments[r]))
			continue
		op = ''.join(c for c in t[3] if not c.isdigit())
		qs, qe, rl = cigar2pos_ops[op](t[3], t[2], rl)
		qint.append([qs, qe])
		if t[2] == '+':
			rint.append([t[0], int(t[1]) - 1, int(t[1]) + rl - 2, '+']) # converted to 0 based coordinates
		else:
			rint.append([t[0], int(t[1]) + rl - 2, int(t[1]) - 1, '-']) # converted to 0 based coordinates
		qual.append(int(t[4]))
	qint_ind = sorted(range(len(qint)), key = lambda i: (qint[i][0], qint[i][1]))
	qint = [qint[i] for i in qint_ind]
	rint = [rint[i] for i in qint_ind]
	qual = [qual[i] for i in qint_ind]
	return (qint, rint, qual)


#def merge_alignment
