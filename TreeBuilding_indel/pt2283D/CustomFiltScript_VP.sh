#!/bin/bash
#SBATCH --job-name=SmurfFilt
#SBATCH -c 8
#SBATCH --mem=75G
#SBATCH --time=4:0:0
#SBATCH -e filt.err
#SBATCH -o filt.log

mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2283D/"
mycoordsfile="${mytreedir}/shared_clonal_muts_coords.txt"
mysmurfvcf="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/intermediate/short_variants/SMuRF/pt2283D/230307_HMF1931_pt2283D.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf"

echo "start" > filt.start

grep "^#" ${mysmurfvcf} > ${mytreedir}/SMuRF.filtered.coords.vcf

grep -f ${mycoordsfile} ${mysmurfvcf} >> ${mytreedir}/SMuRF.filtered.coords.vcf
