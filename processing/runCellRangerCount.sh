#!/bin/bash

# Set up Cell Ranger environment
CELLRANGER_PATH=/share/ScratchGeneral/annsen/pipelines/cellranger-3.0.2
export PATH=$CELLRANGER_PATH:$PATH
source $CELLRANGER_PATH/sourceme.bash

# Set up arguments for run
SAMPLE=H9_RPE_YOUNG_S1
INPUT_DIR=/share/ScratchGeneral/annsen/data/experimental_data/PUBLISHING/RPE_scRNA/RAW
REF_DIR=/share/ClusterShare/thingamajigs/SCCG/data/reference_data/refdata-cellranger-GRCh38-3.0.0
OUTPUT_DIR=/share/ScratchGeneral/annsen/data/experimental_data/CLEAN/RPE_scRNA_AgedStudy

# Run Cell Ranger
cellranger count --id=${SAMPLE} --sample=${SAMPLE} --fastqs=${INPUT_DIR} \
    --localmem=128 --localcores=16 --transcriptome=${REF_DIR} --expect-cells=10000 --nosecondary || exit 1


