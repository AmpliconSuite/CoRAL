#!/bin/bash

# usage: map_nanopore_reads.sh READS_DIRECTORY OUTPUT_DIRECTORY TRHEADS

FASTQ=$1
OUTPUT_DIR=$2
THREADS=$3

echo $FASTQ, $OUTPUT_DIR, $THREADS

# Specify outputs
MAPPED_BAM=${OUTPUT_DIR}/winnowmap_mapping.bam
MAPPED_BAM_FILTERED=${OUTPUT_DIR}/winnowmap_mapping.filtered.bam
MAPPED_BAM_SORTED=${OUTPUT_DIR}/winnowmap_mapping.sorted.bam

# specify parameters
GENOME_REF="/oak/stanford/groups/howchang/users/mgjones/reference/hg38/hg38.fa"
MIN_MAPQ=25
MIN_LENGTH=1000


mkdir -p ${OUTPUT_DIR}

eval "$(conda shell.bash hook)"
conda activate ecdna-sv

# create kmer reference
echo "Creating kmer reference..."
meryl count k=15 output merylDB ${GENOME_REF} threads=${THREADS}
meryl print greater-than distinct=0.9998 merylDB > ${OUTPUT_DIR}/repetitive_k15.txt

# align sequences
echo "Aligning with Winnowmap2..."
winnowmap -W ${OUTPUT_DIR}/repetitive_k15.txt \
    -ax map-ont \
    ${GENOME_REF} \
    ${FASTQ} | \
    samtools view -bS -@ ${THREADS} > $MAPPED_BAM

# filter by sequence length and mapping quality
echo "Filtering alignments by length and mapping quality..."
samtools view -h $MAPPED_BAM | \
    awk 'length($10) > 1000 || $1 ~ /^@/'|
    samtools view -bSq $MIN_MAPQ -@ ${THREADS} > $MAPPED_BAM_FILTERED


# convert sam to bam
echo "Sorting..."
samtools sort -@ ${THREADS} $MAPPED_BAM_FILTERED > $MAPPED_BAM_SORTED
samtools index $MAPPED_BAM_SORTED

echo "DONE"
