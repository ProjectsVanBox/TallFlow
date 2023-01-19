#!/bin/bash

myptatodir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO/pt2322/indels/pt2322"
mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2322/"

if [ ! -d ${mytreedir} ]; then
  mkdir -p ${mytreedir}
fi

# Search for shared "clonal" mutations
zgrep -v "^#" ${myptatodir}/*/*.vcf.gz | \
 cut -f 1,2,8 | \
 cut -f 2- -d ':' | \
 sort | \
 uniq -c | \
 grep -v ";CLONAL_SAMPLES=1" | \
 grep -oP ";CLONAL_SAMPLE_NAMES=[\w,-]+" | \
 sort | \
 uniq -c | \
 sort -k1n \
 > ${mytreedir}/shared_clonal_muts_info.txt
