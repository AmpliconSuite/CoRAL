#!/usr/bin/env bash
## Basic script employing cnvkit.py to call CNVs from alignments.

## USAGE
## call_cnvs.sh $BAM $REFERENCE_CNN $OUTPUT_DIR

BAM=$1
REFERENCE_CNN=$2
OUTPUT_DIR=$3


cnvkit.py batch $BAM \
    --seq-method wgs \
    --drop-low-coverage \
    --reference $REFERENCE_CNN \
    --scatter --diagram \
    -d $OUTPUT_DIR \
