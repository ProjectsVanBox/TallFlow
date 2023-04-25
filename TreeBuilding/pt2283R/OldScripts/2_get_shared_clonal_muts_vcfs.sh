#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283R"
myinfofile="${mytreedir}/shared_clonal_muts_info.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/3_Output/BulkSeq_SMuRF/pt2283R/230307_HMF1931_pt2283R.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf"

grep "^#" ${mysmurfvcf} > ${mytreedir}/header.vcf

while read line; do
  line_grep=$(echo ${line} | awk '{ print $2}' )
  nmuts=$(echo ${line} | awk '{ print $1}' )
  mysamples=$(echo ${line} | cut -f 2 -d'=' | tr ',' '_' | sed 's/pt2283R-DX1BM-//g')
  mysharedmutsvcf=${mytreedir}/muts_vcfs/${nmuts}_${mysamples}.vcf
  cat ${mytreedir}/header.vcf > ${mysharedmutsvcf}
  grep ${line_grep} ${mysmurfvcf} | awk '{ if (length($4) == 1 && length($5) == 1) { print $0} }'>> ${mysharedmutsvcf}
done < ${myinfofile}
