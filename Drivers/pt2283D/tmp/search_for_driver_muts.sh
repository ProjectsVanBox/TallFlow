#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Drivers/pt2283D

RAW_VCF=/hpc/pmc_vanboxtel/processed/230307_HMF1931_pt2283D/VCFS/VCF/230307_HMF1931_pt2283D.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf
PTATO_MERGED_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2238D/PTATO_vcfs/pt2283D/pt2283D.ptato.merged.vcf.gz
SMURF_FILTERED_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/intermediate/short_variants/SMuRF/pt2283D/230307_HMF1931_pt2283D.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf


MOD_HIGH_VCF=moderate_high.vcf
#grep -P "MODERATE|HIGH" ${RAW_VCF} > ${MOD_HIGH_VCF}

PTATO_MOD_HIGH_VCF=ptato_moderate_high.vcf
zgrep -f <( cut -f 1,2 ${MOD_HIGH_VCF} ) ${PTATO_MERGED_VCF} > ${PTATO_MOD_HIGH_VCF}

SMURF_MOD_HIGH_VCF=SMuRF_moderate_high.vcf
grep -f <( cut -f 1,2 ${PTATO_MOD_HIGH_VCF} ) ${SMURF_FILTERED_VCF} > ${SMURF_MOD_HIGH_VCF}


MOD_HIGH_SMURF_VCF=moderate_high_SMuRF.vcf
grep -f <( cut -f 1,2 ${SMURF_MOD_HIGH_VCF} ) ${MOD_HIGH_VCF} > ${MOD_HIGH_SMURF_VCF}


