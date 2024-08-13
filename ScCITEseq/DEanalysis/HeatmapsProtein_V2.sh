#!/bin/bash

#SBATCH --job-name=ProtHeatmaps
#SBATCH --mem=35G
#SBATCH --time=0:30:0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


module load R/4.1.2
Rscript HeatmapsProtein_V2.R

exit
