#!/bin/bash

#SBATCH --job-name=ProtHeatmaps
#SBATCH --mem=20G
#SBATCH --time=1:0:0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


module load R/4.1.2
Rscript HeatmapsProtein.R

exit