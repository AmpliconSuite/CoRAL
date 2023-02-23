"""
Convert cigar strings from SA:Z: field into chimeric alignments
A chimeric alignment is represented by the following 3 lists
qint: 0-based intervals on query (reads) for each local alignment
rint: 0-based intervals on reference genome for each local alignment
qual: mapping qualities for each local alignment

TBD:
implement merge_alignment
"""
import logging
import time
from typing import Tuple, List

from long_read_aa.breakpoints import global_names


def cigar2posSM(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *S*M into chimeric alignment

	Args:
		cigar: CIGAR string in format *S*M
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = read_length - 1
	else:
		qs = 0
		qe = al - 1
	return qs, qe, al


def cigar2posMS(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *M*S into chimeric alignment

	Args:
		cigar: CIGAR string in format *M*S
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[: cigar.index('M')])
	if strand == '+':
		qs = 0
		qe = al - 1
	else:
		qs = int(cigar[cigar.index('M') + 1 : cigar.index('S')])
		qe = read_length - 1
	return qs, qe, al


def cigar2posSMS(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *S*M*S into chimeric alignment

	Args:
		cigar: CIGAR string in format *S*M*S
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = qs + al - 1
	else:
		qs = int(cigar[cigar.index('M') + 1 : -1])
		qe = qs + al - 1
	return qs, qe, al


def cigar2posSMD(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *S*M*D into chimeric alignment

	Args:
		cigar: CIGAR string in format *S*M*D
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[cigar.index('S') + 1 : cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = read_length - 1
	else:
		qs = 0
		qe = int(cigar[cigar.index('S') + 1 : cigar.index('M')]) - 1
	return qs, qe, al


def cigar2posMDS(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *M*D*S into chimeric alignment

	Args:
		cigar: CIGAR string in format *M*D*S
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[: cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		qs = 0
		qe = int(cigar[: cigar.index('M')]) - 1
	else:
		qs = int(cigar[cigar.index('D') + 1 : cigar.index('S')])
		qe = read_length - 1
	return qs, qe, al


def cigar2posSMDS(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *S*M*D*S into chimeric alignment

	Args:
		cigar: CIGAR string in format *S*M*D*S
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[cigar.index('S') + 1: cigar.index('M')]) + int(cigar[cigar.index('M') + 1 : cigar.index('D')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = read_length - int(cigar[cigar.index('D') + 1 : -1]) -1
	else:
		qs = int(cigar[cigar.index('D') + 1 : -1])
		qe = read_length - int(cigar[: cigar.index('S')]) -1
	return qs, qe, al


def cigar2posSMI(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *S*M*I into chimeric alignment

	Args:
		cigar: CIGAR string in format *S*M*I
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = read_length - 1
	else:
		qs = 0
		qe = read_length - int(cigar[: cigar.index('S')]) - 1
	return qs, qe, al


def cigar2posMIS(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *M*I*S into chimeric alignment

	Args:
		cigar: CIGAR string in format *M*I*S
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[: cigar.index('M')])
	if strand == '+':
		qs = 0
		qe = read_length - int(cigar[cigar.index('I') + 1 : cigar.index('S')]) - 1
	else:
		qs = int(cigar[cigar.index('I') + 1 : cigar.index('S')])
		qe = read_length - 1
	return qs, qe, al


def cigar2posSMIS(cigar, strand, read_length) -> Tuple[int, int, int]:
	"""
	Convert cigar string in format *S*M*I*S into chimeric alignment

	Args:
		cigar: CIGAR string in format *S*M*I*S
		strand: Strand of the read ("+" or "-")
		read_length: Read length
	Returns:
		The start and end position on the read on positive strand, and the alignment length on
			the reference genome.
	"""
	al = int(cigar[cigar.index('S') + 1 : cigar.index('M')])
	if strand == '+':
		qs = int(cigar[: cigar.index('S')])
		qe = read_length - int(cigar[cigar.index('I') + 1 : -1]) - 1
	else:
		qs = int(cigar[cigar.index('I') + 1 : -1])
		qe = read_length - int(cigar[: cigar.index('S')]) -1
	return qs, qe, al

def alignment_from_satags(sa_list: List[str], read_length) -> Tuple[list, list, list]:
	"""
	Convert "SA:Z" a list of strings into a new chimeric alignment. 
	Require at least one (soft) clip and one match for each canonical alignment record in a chimeric alignment
		If not, trigger a warning message in logger

	Args:
		sa_list: A list of "SA:Z" tags from bam
		read_length: Read length
	Returns:
		chimeric alignment in the form of qint, rint and qual list
		Alignments sorted according to the starting positions on the read on positive strand
	"""

	#Dict indicating which cigar2pos operation will be called 
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

	qint, rint, qual = [], [], []
	for sa in sa_list:
		t = sa.split(',')
		if 'S' not in t[3] or 'M' not in t[3]:
			# Require a chimeric alignment record having at least some (soft)clips and matches 
			logging.warning("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Found chimeric alignment without match or soft clips.")
			logging.warning("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tRead name: %s; Read length: %d." %(r, read_length))
			logging.warning("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tAll CIGAR strings: %s." %(self.chimeric_alignments[r]))
			continue
		op = ''.join(c for c in t[3] if not c.isdigit())
		qs, qe, al = cigar2pos_ops[op](t[3], t[2], read_length)
		qint.append([qs, qe])
		if t[2] == '+':
			rint.append([t[0], int(t[1]) - 1, int(t[1]) + al - 2, '+']) # converted to 0 based coordinates
		else:
			rint.append([t[0], int(t[1]) + al - 2, int(t[1]) - 1, '-']) # converted to 0 based coordinates
		qual.append(int(t[4]))
	qint_ind = sorted(range(len(qint)), key = lambda i: (qint[i][0], qint[i][1]))
	qint = [qint[i] for i in qint_ind]
	rint = [rint[i] for i in qint_ind]
	qual = [qual[i] for i in qint_ind]
	return (qint, rint, qual)

