#!/bin/bash

module load bedtools bcftools

for VCF in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/*/intermediate/svs/Integration/*/*/pt2322*TALL**filtered.vcf; do
  echo ${VCF}
  PT=$( echo ${VCF} | grep -oP "pt\d+D*-\w+-TALL" | cut -f 1 -d '-' )
  echo ${PT}	
   
  mkdir /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/${PT}

  cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/${PT}

  OUT=$(basename ${VCF})
  OUT=${OUT/%.vcf/_geneAnnotation.vcf}
  bcftools annotate -a /hpc/pmc_vanboxtel/projects/TallFlow/2_Code/RSSmotifSearch/genes.sorted.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') ${VCF} > ${OUT}

  OUT1=${OUT/%.vcf/_intersect.vcf}
  grep "^#" ${VCF} > ${OUT1}
  bedtools intersect -a ${OUT} -b /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/*/intermediate/svs/Integration/*/*/${PT}*TALL**filtered.vcf -wa -c >> ${OUT1}

  OUT2=${OUT1/%.vcf/_uniq.vcf}
  grep -P "^#|\s+1$" ${OUT1} > ${OUT2}
 
#  OUT3=${OUT1/%.vcf/_shared.vcf}
#  grep "^#" ${VCF} > ${OUT3}
#  grep -vP "^#|\s+1$" ${OUT1} >> ${OUT3}
done
