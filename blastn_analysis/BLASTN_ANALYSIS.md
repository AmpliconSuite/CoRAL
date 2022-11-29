Step 1: ```python(3) refine_breakpoint_graph.py 
        --sr_bam short_reads.bam 
        --lr_bam long_reads.bam 
        --aa_graph AA_graph.txt 
        --aa_cycle AA_cycles.txt```
        which will output a tsv file including all short read and long read breakpoints.

Step 2: Run step1 for all amplicons from a sample and collect all resulting tsv files into a single directory. Then run ```python extract_substring.py <reference_genome> <sample_name> <the directory containing all tsvs from one sample> <coverage threshold>```

The coverage threshold should a floating point number so that every breakpoint with the number of reads below this threshold will not be considered. 

The output from this step is a fasta file containing all junction genomes having 0 long reads coverage from step 1, which shall be used the as the subject in blastn queries. 

Step 3: Run ```fq2fa.py *.fastq``` to prepare a query fasta file containing all long reads. (Make sure the file has extension ```.fastq``` rather than ```.fq```, otherwise the fq file may be overwritten.)

Run ```blastn -query query.fasta -subject subject.fasta -out <sample name>_alignment.txt -outfmt 6``` to compute the local alignments.

Step 4: Run ```analyze_blast.py <sample name>_alignment.txt``` for a very crude analysis (to reproduce the numbers in the spreadsheet).  
