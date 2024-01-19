#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1180/"
myinfofile="${mytreedir}/shared_clonal_muts_info.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt1180/intermediate/short_variants/SMuRF/pt1180/231016_HMFreg2090_pt1180.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf"

grep "^#" ${mysmurfvcf} > ${mytreedir}/header.vcf

while read line; do
  line_grep=$(echo ${line} | awk '{ print $2}' )
  nmuts=$(echo ${line} | awk '{ print $1}' )
  mysamples=$(echo ${line} | cut -f 2 -d'=' | tr ',' '_' | sed 's/pt1180-DX1PB-//g')
  mysharedmutsvcf=${mytreedir}/muts_vcfs/${nmuts}_${mysamples}.vcf
  cat ${mytreedir}/header.vcf > ${mysharedmutsvcf}
  grep ${line_grep} ${mysmurfvcf} | awk '{ if (length($4) == 1 && length($5) == 1) { print $0} }'>> ${mysharedmutsvcf}
done < ${myinfofile}
