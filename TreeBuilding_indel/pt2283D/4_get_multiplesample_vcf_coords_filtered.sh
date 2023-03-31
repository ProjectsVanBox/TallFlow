#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2283D/"
mycoordsfile="${mytreedir}/shared_clonal_muts_coords_real_indels.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/intermediate/short_variants/SMuRF/pt2283D/230307_HMF1931_pt2283D.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf"

grep "^#" ${mysmurfvcf} > ${mytreedir}/SMuRF.filtered.coords.vcf

grep -f ${mycoordsfile} ${mysmurfvcf} >> ${mytreedir}/SMuRF.filtered.coords.vcf
