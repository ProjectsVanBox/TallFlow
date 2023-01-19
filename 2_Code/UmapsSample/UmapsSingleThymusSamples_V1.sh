#!/bin/bash
#SBATCH --job-name=UmapsSingleThymusSamples_V1.R
#SBATCH --output=/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/UmapsSample/UmapsSingleThymusSamples_V1.R.out
#SBATCH --error=/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/UmapsSample/UmapsSingleThymusSamples_V1.R.err
#SBATCH --partition=cpu
#SBATCH --time=23:59:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 16
#SBATCH --mem=40G
#SBATCH --gres=tmpspace:300G
#SBATCH --nodes=1
#SBATCH --open-mode=append


module load R/4.1.2
Rscript /hpc/pmc_vanboxtel/projects/TallFlow/2_Code/UmapsSample/UmapsSingleThymusSamples_V1.R