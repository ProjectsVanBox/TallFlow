#!/bin/bash

#SBATCH --job-name=prediction_CelltypeCite
#SBATCH --mem=30G
#SBATCH --time=2:00:0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


module load R/4.1.2
Rscript Prediction_celltypesVSGating.R

exit
