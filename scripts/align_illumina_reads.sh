#!/bin/bash

# usage: map_illumina_reads.sh R1 R2 ID OUTPUT_DIRECTORY TRHEADS

R1=$1
R2=$2
ID=$3
OUTPUT_DIR=$4
THREADS=$5

PICARD="/home/users/mgjones/software/picard.jar"
BWA="/home/users/mgjones/software/bwa/bwa"
TRIMMOMATIC="/home/users/mgjones/software/Trimmomatic-0.39/trimmomatic-0.39.jar"

echo $R1, $R2, $ID, $OUTPUT_DIR, $THREADS

GENOME_REF="/oak/stanford/groups/howchang/users/mgjones/reference/hg38/hg38.fa"

R1_PAIRED=${OUTPUT_DIR}/trimmomatic/paired/${ID}_R1.trim.fastq.gz
R2_PAIRED=${OUTPUT_DIR}/trimmomatic/paired/${ID}_R2.trim.fastq.gz
R1_UNPAIRED=${OUTPUT_DIR}/trimmomatic/unpaired/${ID}_R1.unpaired.trim.fastq.gz
R2_UNPAIRED=${OUTPUT_DIR}/trimmomatic/unpaired/${ID}_R2.unpaired.trim.fastq.gz
SORTBAM=${OUTPUT_DIR}/${ID}.sort.bam
RMDUP_BAM=${OUTPUT_DIR}/${ID}.sort.rmdup.bam
DEDUP_METRICS=${OUTPUT_DIR}/${ID}.markduplicates_metrics.txt
ALIGN_METRICS=${OUTPUT_DIR}/${ID}.alignment_metrics.txt
INSERT_METRICS=${OUTPUT_DIR}/${ID}.insert_metrics.txt
INSERT_HISTOGRAM=${OUTPUT_DIR}/${ID}.insert_size_histogram.pdf

# set up directory
mkdir -p ${OUTPUT_DIR}

LOG=${OUTPUT_DIR}/alignment.log

eval "$(conda shell.bash hook)"
conda activate ecdna-sv

echo "Trimming reads..." >> $LOG
mkdir -p ${OUTPUT_DIR}/trimmomatic/paired
mkdir -p ${OUTPUT_DIR}/trimmomatic/unpaired

java -jar $TRIMMOMATIC PE $R1 $R2 \
	$R1_PAIRED $R1_UNPAIRED \
	$R2_PAIRED $R2_UNPAIRED \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
	MINLEN:20 ILLUMINACLIP:/oak/stanford/groups/howchang/users/mgjones/reference/adapters.fa:2:30:10 \
	-threads $THREADS &> $LOG

echo "Aligning reads..."  >> $LOG
$BWA mem -t $THREADS \
	-R "@RG\tID:${ID}\tLB:${ID}\tPL:ILLUMINA\tSM:${ID}" \
	$GENOME_REF $R1_PAIRED $R2_PAIRED | \
	samtools view -Sb - | \
	samtools sort -@ $THREADS -o $SORTBAM \
	&> $LOG

# remove duplicatepicard MarkDuplicates INPUT=$SORTBAM OUTPUT=$RMDUP_BAM \
echo "Mark and remove duplicates..."  >> $LOG
java -jar $PICARD MarkDuplicates INPUT=$SORTBAM OUTPUT=$RMDUP_BAM \
    VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true \
    METRICS_FILE=$DEDUP_METRICS \
	&> $LOG
samtools index $RMDUP_BAM

#alignment metrics
echo "Record alignmnet statistics..." >> $LOG
java -jar $PICARD CollectAlignmentSummaryMetrics \
	R=$GENOME_REF \
	I=$RMDUP_BAM \
	O=$ALIGN_METRICS

java -jar $PICARD CollectInsertSizeMetrics \
	INPUT=$RMDUP_BAM \
	OUTPUT=$INSERT_METRICS \
	HISTOGRAM_FILE=$INSERT_HISTOGRAM

echo 'DONE' >> $LOG

conda deactivate
