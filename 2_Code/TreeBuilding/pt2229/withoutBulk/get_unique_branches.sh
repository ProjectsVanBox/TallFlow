#!/bin/bash

myptatodir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2229/snvs/pt2229"
mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2229_new/withoutBulk/branches_vcfs"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2229/intermediate/short_variants/SMuRF/pt2229/220919_HMFreg1764_pt2229.SMuRF.filtered.vcf"

if [ ! -d ${mytreedir} ]; then
  mkdir -p ${mytreedir}
fi

for VCF in ${myptatodir}/*/*.vcf.gz; do 
    SAMPLE=$( echo ${VCF} | grep -oP "TALL\d+.ptato" | cut -f1 -d'.')
    if [ ! -z ${SAMPLE} ]; then
        zgrep "^#" ${mysmurfvcf} > ${mytreedir}/${SAMPLE}_branch.vcf
        grep -f <( zgrep -P ";CLONAL_SAMPLES=1;" ${VCF} | grep -P ";ABSENT_SAMPLES=(9|10)" | cut -f 1,2 ) ${mysmurfvcf} >> ${mytreedir}/${SAMPLE}_branch.vcf
    fi
    
#    /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO/pt2322/snvs/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-ALLBULK/220919_HMFreg1764_pt2322_pt2322-DX1BM-ALLBULK.ptato.filtered.vcf.gz
    
done

# Search for shared "clonal" mutations
#zgrep -v "^#" ${myptatodir}/*/*.vcf.gz | \
# cut -f 1,2,8 | \
# cut -f 2- -d ':' | \
# sort | \
# uniq -c | \
# grep -P ";CLONAL_SAMPLES=1;" | \
# sort | \
# uniq -c | \
# sort -k1n \
# > ${mytreedir}/unique_clonal_muts_info.txt
