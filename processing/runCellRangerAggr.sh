#!/bin/bash

# Set up Cell Ranger environment
CELLRANGER_PATH=/share/ScratchGeneral/annsen/pipelines/cellranger-3.0.2
export PATH=$CELLRANGER_PATH:$PATH
source $CELLRANGER_PATH/sourceme.bash

# Set up arguments for run
INPUT_CSV=aggr.csv
ID=H9_RPE_AGGR
OUTPUT_DIR=/share/ScratchGeneral/annsen/data/experimental_data/CLEAN/RPE_scRNA_AgedStudy

cd $OUTPUT_DIR

# Run Cell Ranger
cellranger aggr --id=${ID} --csv=$INPUT_CSV --normalize=mapped
