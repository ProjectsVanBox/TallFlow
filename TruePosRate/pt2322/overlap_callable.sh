#!/bin/bash

SAMPLE=$1

MISSING_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TruePosRate/pt2322/${SAMPLE}_MISSING.vcf
CL_BED=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/QC/CallableLoci/pt2322/pt2322-DX1BM-${SAMPLE}.callableloci.bed

grep "^#" ${MISSING_VCF} > ${SAMPLE}_MISSING_Callable.vcf
bedtools intersect -a ${MISSING_VCF} -b <(grep "CALLABLE" ${CL_BED}) -wa >> ${SAMPLE}_MISSING_Callable.vcf

grep "^#" ${MISSING_VCF} > ${SAMPLE}_MISSING_NotCallable.vcf
bedtools intersect -a ${MISSING_VCF} -b <(grep -v "CALLABLE" ${CL_BED}) -wa >> ${SAMPLE}_MISSING_NotCallable.vcf

