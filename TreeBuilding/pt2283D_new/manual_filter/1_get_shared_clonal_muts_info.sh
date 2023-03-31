#!/bin/bash

myptatodir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/snvs/manual_filter"
mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283D_new/manual_filter"

if [ ! -d ${mytreedir} ]; then
  mkdir -p ${mytreedir}
fi

# Search for shared "clonal" mutations
grep -v "^#" ${myptatodir}/*.vcf | \
 cut -f 1,2,8 | \
 cut -f 2- -d ':' | \
 sort | \
 uniq -c | \
 grep -v ";CLONAL_SAMPLES=1[;,]" | \
 grep -oP ";CLONAL_SAMPLE_NAMES=[\w,-]+" | \
 sort | \
 uniq -c | \
 sort -k1n \
 > ${mytreedir}/shared_clonal_muts_info.txt
