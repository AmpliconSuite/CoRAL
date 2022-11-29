import sys

if __name__ == '__main__':
	fp = open(sys.argv[1], 'r')
	
	blast_hit = dict()
	for line in fp:
		s = line.strip().split('\t')
		if ((int(s[8]) <= 9900 and int(s[9]) >= 10100) \
			or (int(s[9]) <= 9900 and int(s[8]) >= 10100)) \
			and (float(s[10]) <= 1e-30):
			try:
				blast_hit[s[1]].add(s[0])
			except:
				blast_hit[s[1]] = set([s[0]])
		
	for g in blast_hit.keys():
		blast_hit[g] = len(blast_hit[g])
	print blast_hit
	#print ("Total # unique reads - %d." % i)



