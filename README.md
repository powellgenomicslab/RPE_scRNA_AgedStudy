# RPE_scRNA_AgedStudy
Public repository for code related to the generation and analysis of iPSC-derived retina pigmented epithelium cells by Lidgerwood et al.

## Sample Information

## Processing
Data was processed using the 10x Genomics Cell Ranger 3.0.2 pipeline.

### Requirements
#### System Requirements
- 8-core Intel or AMD processor (16 cores recommended)
- 64GB RAM (128GB recommended)
- 1TB free disk space
- 64-bit CentOS/RedHat 6.0 or Ubuntu 12.04

#### Software Requirements
- [Cell Ranger version 3.0.2](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0)
- conda
- Snakemake

#### Additional files
- [*Homo sapiens* GRCh38 Cell Ranger reference](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

### Pre-processing of FASTQ.GZ files
Sequence from four sequencing lanes were combined for submission to ArrayExpress. They need to be separated again in order to be run through the Cell Ranger pipeline. This can be done via the [separateFastqs.sh](processing/separateFastqs.sh) script.





