#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Drivers/pt2229

RAW_VCF=/hpc/pmc_vanboxtel/processed/220919_HMFreg1764_pt2229/VCFS/VCF/220919_HMFreg1764_pt2229.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf
PTATO_MERGED_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO/pt2229/ptato_vcfs/pt2229/pt2229.ptato.merged.vcf.gz
SMURF_FILTERED_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO/pt2229/intermediate/short_variants/SMuRF/pt2229/220919_HMFreg1764_pt2229.SMuRF.filtered.vcf

MOD_HIGH_VCF=moderate_high.vcf
grep -P "MODERATE|HIGH" ${RAW_VCF} > ${MOD_HIGH_VCF}

PTATO_MOD_HIGH_VCF=ptato_moderate_high.vcf
zgrep -f <( cut -f 1,2 ${MOD_HIGH_VCF} ) ${PTATO_MERGED_VCF} > ${PTATO_MOD_HIGH_VCF}

SMURF_MOD_HIGH_VCF=SMuRF_moderate_high.vcf
grep -f <( cut -f 1,2 ${PTATO_MOD_HIGH_VCF} ) ${SMURF_FILTERED_VCF} > ${SMURF_MOD_HIGH_VCF}


MOD_HIGH_SMURF_VCF=moderate_high_SMuRF.vcf
grep -f <( cut -f 1,2 ${SMURF_MOD_HIGH_VCF} ) ${MOD_HIGH_VCF} > ${MOD_HIGH_SMURF_VCF}


