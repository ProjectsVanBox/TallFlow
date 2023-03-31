#!/bin/bash

VCF=$1
FINGERPRINT_VCF=$2

#VCF=/hpc/pmc_vanboxtel/projects/TallFlow/1_Input/VCF/pt2283D/230307_HMF1931_pt2283D.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf
VCF=/hpc/pmc_vanboxtel/projects/TallFlow/1_Input/VCF/pt2283R/230307_HMF1931_pt2283R.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf

FINGERPRINT_VCF=/hpc/pmc_vanboxtel/personal/mroosmalen/FingerprintTool/resources/81_snps_mip_design_nijmegen_sort_hg38_liftover.vcf


OUT=${VCF/%.vcf/_fingerprint.txt}

grep -m 1 "^#CHROM" ${VCF} | cut -f 1,2,4,5,10- > ${OUT}

grep -f <(grep -v "^#" ${FINGERPRINT_VCF} | cut -f 1,2 | awk '{ print "^"$1"\t"$2"\t"}') ${VCF} | cut -f 1,2,4,5,10- >> ${OUT}