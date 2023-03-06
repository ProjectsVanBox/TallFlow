#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Telomeres/pt2229

for BAM in /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/BULK_SEQ/BAMS/pt2229/*ALLBULK_*bam; do
 sbatch /hpc/pmc_vanboxtel/projects/TallFlow/2_Code/Telomeres/sbatch_telomerecat.sh ${BAM}
done