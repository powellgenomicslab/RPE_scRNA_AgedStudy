#!/bin/bash
# Sample name
SAMPLE=$1

if [[  $SAMPLE == "H9_RPE_YOUNG" ]]; then
    SAMPLE=${SAMPLE}_S1
elif [[ $SAMPLE == "H9_RPE_AGED" ]]; then
    SAMPLE=${SAMPLE}_S2
fi

# LOOP OVER EACH LANE
for LANE in 1..4
do
    zcat ${SAMPLE}_I1_001.fastq.gz | grep -A 3 \:1\: | gzip > ${SAMPLE}_L00${LANE}_I1_001.fastq.gz
    zcat ${SAMPLE}_R1_001.fastq.gz | grep -A 3 \:1\: | gzip > ${SAMPLE}_L00${LANE}_R1_001.fastq.gz
    zcat ${SAMPLE}_R2_001.fastq.gz | grep -A 3 \:1\: | gzip > ${SAMPLE}_L00${LANE}_R2_001.fastq.gz
done