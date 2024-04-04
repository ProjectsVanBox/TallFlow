#!/bin/bash

BAM=$1
FASTA=/hpc/pmc_vanboxtel/projects/TallFlow/1_Input/FreecBAF/freec/Homo_sapiens_assembly38.fasta
OUT=${BAM/.bam/.mpileup}
OUT=$( basename ${OUT} )

singularity exec -B /hpc /hpc/local/Rocky8/pmc_vanboxtel/singularity_cache/depot.galaxyproject.org-singularity-samtools-1.17--h00cdaf9_0.img \
samtools mpileup --fasta-ref ${FASTA} --output ${OUT} ${BAM}



