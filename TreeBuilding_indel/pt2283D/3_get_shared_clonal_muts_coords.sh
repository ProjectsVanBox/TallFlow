#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2283D/"

cat ${mytreedir}/muts_vcfs/*.vcf | \
grep -v "^#" | \
cut -f 1,2 | \
sort | \
uniq > ${mytreedir}/shared_clonal_muts_coords.txt
