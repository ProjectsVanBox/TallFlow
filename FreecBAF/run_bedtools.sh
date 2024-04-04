#!/bin/bash

#SBATCH --time=2:0:0

module load bedtools

TCR_GENES_BED=/hpc/pmc_vanboxtel/projects/TallFlow/1_Input/FreecBAF/TCR_genes.bed
BAM=$1
OUT=${BAM/.bam/.cov}
OUT=$( basename ${OUT} )

bedtools coverage -a ${TCR_GENES_BED} -b ${BAM} -mean > ${OUT}