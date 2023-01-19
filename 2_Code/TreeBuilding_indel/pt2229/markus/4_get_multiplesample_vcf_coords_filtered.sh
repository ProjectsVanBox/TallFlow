#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2229/"
mycoordsfile="${mytreedir}/shared_clonal_muts_coords.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2229/intermediate/short_variants/SMuRF/pt2229/220919_HMFreg1764_pt2229.SMuRF.filtered.vcf"

grep "^#" ${mysmurfvcf} > ${mytreedir}/SMuRF.filtered.coords.vcf

grep -f ${mycoordsfile} ${mysmurfvcf} >> ${mytreedir}/SMuRF.filtered.coords.vcf
