#!/bin/bash

# FOR YOUNG SAMPLES
# Split Index
zcat H9_RPE_YOUNG_S1_I1_001.fastq.gz | grep -A 3 \:1\: | gzip > H9_RPE_YOUNG_S1_L001_I1_001.fastq.gz
zcat H9_RPE_YOUNG_S1_I1_001.fastq.gz | grep -A 3 \:2\: | gzip > H9_RPE_YOUNG_S1_L002_I1_001.fastq.gz
zcat H9_RPE_YOUNG_S1_I1_001.fastq.gz | grep -A 3 \:3\: | gzip > H9_RPE_YOUNG_S1_L003_I1_001.fastq.gz
zcat H9_RPE_YOUNG_S1_I1_001.fastq.gz | grep -A 3 \:4\: | gzip > H9_RPE_YOUNG_S1_L004_I1_001.fastq.gz

# Split R1
zcat H9_RPE_YOUNG_S1_R1_001.fastq.gz | grep -A 3 \:1\: | gzip > H9_RPE_YOUNG_S1_L001_R1_001.fastq.gz
zcat H9_RPE_YOUNG_S1_R1_001.fastq.gz | grep -A 3 \:2\: | gzip > H9_RPE_YOUNG_S1_L002_R1_001.fastq.gz
zcat H9_RPE_YOUNG_S1_R1_001.fastq.gz | grep -A 3 \:3\: | gzip > H9_RPE_YOUNG_S1_L003_R1_001.fastq.gz
zcat H9_RPE_YOUNG_S1_R1_001.fastq.gz | grep -A 3 \:4\: | gzip > H9_RPE_YOUNG_S1_L004_R1_001.fastq.gz

# Split R2
zcat H9_RPE_YOUNG_S1_R2_001.fastq.gz | grep -A 3 \:1\: | gzip > H9_RPE_YOUNG_S1_L001_R2_001.fastq.gz
zcat H9_RPE_YOUNG_S1_R2_001.fastq.gz | grep -A 3 \:2\: | gzip > H9_RPE_YOUNG_S1_L002_R2_001.fastq.gz
zcat H9_RPE_YOUNG_S1_R2_001.fastq.gz | grep -A 3 \:3\: | gzip > H9_RPE_YOUNG_S1_L003_R2_001.fastq.gz
zcat H9_RPE_YOUNG_S1_R2_001.fastq.gz | grep -A 3 \:4\: | gzip > H9_RPE_YOUNG_S1_L004_R2_001.fastq.gz

# FOR AGED SAMPLES
# Split Index
zcat H9_RPE_AGED_S2_I1_001.fastq.gz | grep -A 3 \:1\: | gzip > H9_RPE_AGED_S2_L001_I1_001.fastq.gz
zcat H9_RPE_AGED_S2_I1_001.fastq.gz | grep -A 3 \:2\: | gzip > H9_RPE_AGED_S2_L002_I1_001.fastq.gz
zcat H9_RPE_AGED_S2_I1_001.fastq.gz | grep -A 3 \:3\: | gzip > H9_RPE_AGED_S2_L003_I1_001.fastq.gz
zcat H9_RPE_AGED_S2_I1_001.fastq.gz | grep -A 3 \:4\: | gzip > H9_RPE_AGED_S2_L004_I1_001.fastq.gz

# Split R1
zcat H9_RPE_AGED_S2_R1_001.fastq.gz | grep -A 3 \:1\: | gzip > H9_RPE_AGED_S2_L001_R1_001.fastq.gz
zcat H9_RPE_AGED_S2_R1_001.fastq.gz | grep -A 3 \:2\: | gzip > H9_RPE_AGED_S2_L002_R1_001.fastq.gz
zcat H9_RPE_AGED_S2_R1_001.fastq.gz | grep -A 3 \:3\: | gzip > H9_RPE_AGED_S2_L003_R1_001.fastq.gz
zcat H9_RPE_AGED_S2_R1_001.fastq.gz | grep -A 3 \:4\: | gzip > H9_RPE_AGED_S2_L004_R1_001.fastq.gz

# Split R2
zcat H9_RPE_AGED_S2_R2_001.fastq.gz | grep -A 3 \:1\: | gzip > H9_RPE_AGED_S2_L001_R2_001.fastq.gz
zcat H9_RPE_AGED_S2_R2_001.fastq.gz | grep -A 3 \:2\: | gzip > H9_RPE_AGED_S2_L002_R2_001.fastq.gz
zcat H9_RPE_AGED_S2_R2_001.fastq.gz | grep -A 3 \:3\: | gzip > H9_RPE_AGED_S2_L003_R2_001.fastq.gz
zcat H9_RPE_AGED_S2_R2_001.fastq.gz | grep -A 3 \:4\: | gzip > H9_RPE_AGED_S2_L004_R2_001.fastq.gz