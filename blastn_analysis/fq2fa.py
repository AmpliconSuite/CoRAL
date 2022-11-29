import sys

if __name__ == '__main__':
	fp = open(sys.argv[1], 'r')
	fp_w = open(sys.argv[1].replace("fastq", "fasta"), 'w')

	i = 0
	for line in fp:
		line = line.strip()
		if i % 4 == 0:
			s = line.split()
			fp_w.write(">%s\n" %s[0][1:])
			#print line
		if i % 4 == 1:
			fp_w.write("%s\n" %line)
		i += 1
	fp.close()
	fp_w.close()