#!/bin/bash

#SBATCH --job-name=HeatmapPerPatient
#SBATCH --mem=40G
#SBATCH --time=2:00:0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


module load R/4.1.2
Rscript HeatmapPerPatient_DEextended_V3.R

exit
