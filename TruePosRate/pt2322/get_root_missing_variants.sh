#!/bin/bash

ROOT_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/withoutBulk/branches_vcfs/root_branch.vcf

SAMPLE=$1
SNV_PTATO_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/snvs/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}.ptato.filtered.vcf.gz
INDEL_PTATO_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/indels/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}.indels.ptato.filtered.vcf.gz

zgrep "^#" ${ROOT_VCF} > ${SAMPLE}_MISSING.vcf

zgrep -v -f <(zgrep -v "^#" ${INDEL_PTATO_VCF} | cut -f1,2) <( zgrep -v -f <(zgrep -v "^#" ${SNV_PTATO_VCF} | cut -f 1,2) <(zgrep -P "^\d+" ${ROOT_VCF} )) >> ${SAMPLE}_MISSING.vcf
