import sys
import os
from Bio.Seq import Seq	

if __name__ == '__main__':
	"""
	Read in chromosome strings
	"""
	genomes = dict()
	chr_string = ""
	chr = ""
	i = 0
	fp = open(sys.argv[1], 'r')
	for line in fp:
		s = line.strip()
		if i % 2 == 0:
			chr = s
		else:
			chr_string = s
			genomes[chr[1:]] = chr_string
		i += 1
	fp.close()


	sample = sys.argv[2]
	dir = sys.argv[3]
	sl = 10000
	gi = 0
	fp_w = open("hg38_" + sample + ".fa", 'w')

	for fn in os.listdir(dir):
		if sample in fn and 'png' not in fn:
			i = 0
			fp = open(dir + fn, 'r')
			for line in fp:
				s = line.strip().split('\t')
				if i > 0:
					if s[-1] == 'N/A' and int(s[6]) == 0 and float(s[5]) > float(sys.argv[4]):
						fp_w.write(">g%d %s,%s,%s,%s,%s\n" %(gi, s[0], s[1], s[2], s[3], s[4]))
						if s[4][0] == '+' and s[4][1] == '-':
							g1 = genomes[s[0]][int(s[1]) - sl : int(s[1]) + 1]
							g2 = genomes[s[2]][int(s[3]) : int(s[3]) + sl]
							fp_w.write("%s%s\n" %(g1, g2))
						if s[4][0] == '+' and s[4][1] == '+':
							g1 = genomes[s[0]][int(s[1]) - sl : int(s[1]) + 1]
							g2 = Seq(genomes[s[2]][int(s[3]) - sl : int(s[3]) + 1]).reverse_complement()
							fp_w.write("%s%s\n" %(g1, g2))
						if s[4][0] == '-' and s[4][1] == '+':
							g1 = Seq(genomes[s[0]][int(s[1]): int(s[1]) + sl]).reverse_complement()
							g2 = Seq(genomes[s[2]][int(s[3]) - sl : int(s[3]) + 1]).reverse_complement()
							fp_w.write("%s%s\n" %(g1, g2))
						if s[4][0] == '-' and s[4][1] == '-':
							g1 = Seq(genomes[s[0]][int(s[1]): int(s[1]) + sl]).reverse_complement()
							g2 = genomes[s[2]][int(s[3]) : int(s[3]) + sl]
							fp_w.write("%s%s\n" %(g1, g2))
						#print s[4][0], s[4][1]
						gi += 1
				i += 1
			fp.close()


	fp_w.close()



	#/ribosome/projects/Nanopore/scripts/breakpoints



