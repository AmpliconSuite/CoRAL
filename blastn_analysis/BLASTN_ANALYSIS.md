Step 1: ```python(3) refine_breakpoint_graph.py 
        --sr_bam short_reads.bam 
        --lr_bam long_reads.bam 
        --aa_graph AA_graph.txt 
        --aa_cycle AA_cycles.txt```
        which will output a tsv file including all short read and long read breakpoints.

Step 2: ```python extract_substring.py <reference_genome> <sample_name> <a directory containing all tsvs from one sample> <coverage threshold>```
