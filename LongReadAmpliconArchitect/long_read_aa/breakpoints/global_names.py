"""
Miscellaneous global variables used within the module
"""
# Time since main function starts. Used for logging output 
TSTART = 0


# Map a strand to its opposite strand
neg_plus_minus = {"+": "-", "-": "+"}


# Sorted chromosome names
chr_idx = {'chr1': 0, 'chr2': 1, 'chr3': 2, 'chr4': 3,
		'chr5': 4, 'chr6': 5, 'chr7': 6, 'chr8': 7,
		'chr9': 8, 'chr10': 9, 'chr11': 10, 'chr12': 11,
		'chr13': 12, 'chr14': 13, 'chr15': 14, 'chr16': 15,
		'chr17': 16, 'chr18': 17, 'chr19': 18, 'chr20': 19,
		'chr21': 20, 'chr22': 21, 'chrX': 22, 'chrY': 23, 'chrM': 24}
