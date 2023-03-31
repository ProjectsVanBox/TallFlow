#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283D_new/manual_filter"

cat ${mytreedir}/muts_vcfs/*.vcf | \
grep -v "^#" | \
cut -f 1,2 | \
sort | \
uniq > ${mytreedir}/shared_clonal_muts_coords.txt
