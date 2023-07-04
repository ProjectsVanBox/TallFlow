#!/bin/bash

module load bedtools

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones

for VCF in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/*/intermediate/svs/Integration/*/*/*TALL**filtered.vcf; do
  echo ${VCF}
  OUT=$(basename ${VCF})
  OUT=${OUT/%.vcf/_geneAnnotation.vcf}
  bcftools annotate -a genes.sorted.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') ${VCF} > ${OUT}

  OUT1=${OUT/%.vcf/_intersect.vcf}
  grep "^#" ${VCF} > ${OUT1}
  bedtools intersect -a ${OUT} -b /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/*/intermediate/svs/Integration/*/*/*TALL**filtered.vcf -wa -c >> ${OUT1}

  OUT2=${OUT1/%.vcf/_uniq.vcf}
  grep -P "^#|\s+1$" ${OUT1} > ${OUT2}
done
