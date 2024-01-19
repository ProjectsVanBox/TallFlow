#!/bin/bash

module load bedtools

RATIO=$1
BED=/hpc/pmc_vanboxtel/projects/TallFlow/1_Input/FreecBAF/TCR_subset.sorted.bed
RATIO_BED=${RATIO/.txt/.bed}

paste <(awk '{ print $1"\t"$2-1 }' ${RATIO}) <( cut -f 2- ${RATIO}) > ${RATIO_BED}

RATIO_BED2=${RATIO_BED/.bed/_subset.bed}

bedtools intersect -a <( tail -n +2 ${RATIO_BED} ) -b ${BED} -wa > ${RATIO_BED2}