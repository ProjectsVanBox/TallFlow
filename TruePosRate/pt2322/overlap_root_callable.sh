#!/bin/bash

SAMPLE=$1

ROOT_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/withoutBulk/branches_vcfs/root_branch.vcf
CL_BED=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/QC/CallableLoci/pt2322/pt2322-DX1BM-${SAMPLE}.callableloci.bed

grep "^#" ${ROOT_VCF} > ${SAMPLE}_root_Callable.vcf
bedtools intersect -a ${ROOT_VCF} -b ${CL_BED} -wa -loj >> ${SAMPLE}_root_Callable.vcf

cut -f 1,2,24- ${SAMPLE}_root_Callable.vcf > ${SAMPLE}_root_Callable.txt
