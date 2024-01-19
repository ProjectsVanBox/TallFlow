#!/bin/bash

#SBATCH --time=1:0:0

module load bcftools R

CELLPHY=/hpc/pmc_vanboxtel/tools/external/cellphy/cellphy.sh
SUPPORT_MAP=/hpc/pmc_vanboxtel/tools/external/cellphy/script/support-map.R
VCF=$1
VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/HaplotypeCaller/pt344_pt3291/VCFS/VCF/pt344_pt3291.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf
OUTGROUP=$2
OUTGROUP=pt344-DX1PB-MONOBULK
MODEL=GT10+FO+E
SEED=2
THREADS=1
PREFIX=${VCF/%.vcf/}
PROB_MSA=off
BS_TREES=100

bash ${CELLPHY} RAXML --msa ${VCF} --msa-format VCF --model ${MODEL} --seed ${SEED} --threads ${THREADS} --prefix ${PREFIX}.Tree --prob-msa ${PROB_MSA}

bash ${CELLPHY} RAXML --bootstrap --msa ${VCF} --model ${MODEL} --seed ${SEED} --threads ${THREADS} --bs-trees ${BS_TREES} --prefix ${PREFIX}.Boot

bash ${CELLPHY} RAXML --support -tree ${PREFIX}.Tree.raxml.bestTree --bs-trees ${PREFIX}.Boot.raxml.bootstraps --prefix ${PREFIX}.Support --threads ${THREADS} --redo

Rscript ${SUPPORT_MAP} ${PREFIX}.Support.raxml.support ${OUTGROUP}
