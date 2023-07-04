#!/bin/bash

module load bedtools

BED=$1
BED=/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/RSSmotifSearch/test.bed
PT=pt2322
SVDIR=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/*/intermediate/svs/gridss/*/*/

bedtools intersect -a ${BED} -b ${SVDIR}/${PT}*TALL7*unfiltered.vcf.gz -filenames 