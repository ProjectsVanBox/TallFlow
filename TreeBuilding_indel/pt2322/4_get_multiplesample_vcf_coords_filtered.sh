#!/bin/bash

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2322/"
mycoordsfile="${mytreedir}/shared_clonal_muts_coords.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO/pt2322/intermediate/short_variants/SMuRF/pt2322/220919_HMFreg1764_pt2322.SMuRF.filtered.vcf"

grep "^#" ${mysmurfvcf} > ${mytreedir}/SMuRF.filtered.coords.vcf

grep -f ${mycoordsfile} ${mysmurfvcf} >> ${mytreedir}/SMuRF.filtered.coords.vcf
