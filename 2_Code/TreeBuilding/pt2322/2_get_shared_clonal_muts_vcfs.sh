#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/"
myinfofile="${mytreedir}/shared_clonal_muts_info.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/intermediate/short_variants/SMuRF/pt2322/220919_HMFreg1764_pt2322.SMuRF.filtered.vcf"

grep "^#" ${mysmurfvcf} > ${mytreedir}/header.vcf

while read line; do
  line_grep=$(echo ${line} | awk '{ print $2}' )
  nmuts=$(echo ${line} | awk '{ print $1}' )
  mysamples=$(echo ${line} | cut -f 2 -d'=' | tr ',' '_' | sed 's/pt2322-DX1BM-//g')
  mysharedmutsvcf=${mytreedir}/muts_vcfs/${nmuts}_${mysamples}.vcf
  cat ${mytreedir}/header.vcf > ${mysharedmutsvcf}
  grep ${line_grep} ${mysmurfvcf} | awk '{ if (length($4) == 1 && length($5) == 1) { print $0} }'>> ${mysharedmutsvcf}
done < ${myinfofile}
