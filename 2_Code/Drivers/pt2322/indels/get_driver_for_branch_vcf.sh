#!/bin/bash

BRANCH_VCF=$1
ANN_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/1_Input/VCF/pt2322/220919_HMFreg1764_pt2322.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf

BRANCH_DRIVERS_VCF=${BRANCH_VCF/.vcf/_MODERATE_HIGH.vcf}
BRANCH_DRIVERS_VCF=$(basename ${BRANCH_DRIVERS_VCF})

grep -f <( grep -v "^#" ${BRANCH_VCF} | cut -f 1,2) <(grep -P "MODERATE|HIGH" ${ANN_VCF}) > ${BRANCH_DRIVERS_VCF}


