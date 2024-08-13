#!/bin/bash

#SBATCH --job-name=HeatmapFlowMarkers
#SBATCH --mem=60G
#SBATCH --time=0:45:0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


module load R/4.1.2
Rscript HeatmapsProtein_flowMarkers.R

exit
