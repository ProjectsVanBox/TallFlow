#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt344/"
mycoordsfile="${mytreedir}/shared_clonal_muts_coords.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt344/intermediate/short_variants/SMuRF/pt344/230612_HMFreg2010_pt344.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf"

grep "^#" ${mysmurfvcf} > ${mytreedir}/SMuRF.filtered.coords.vcf

grep -f ${mycoordsfile} ${mysmurfvcf} >> ${mytreedir}/SMuRF.filtered.coords.vcf
