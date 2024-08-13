#!/bin/bash

#SBATCH --job-name=HeatmapCelltypes
#SBATCH --mem=20G
#SBATCH --time=2:0:0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


Module load R/4.1.2
Rscript heatmapcelltypes.R

exit