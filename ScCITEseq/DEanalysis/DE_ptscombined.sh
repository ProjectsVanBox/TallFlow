#!/bin/bash

#SBATCH --job-name=DE_ptsCombined
#SBATCH --mem=60G
#SBATCH --time=2:00:0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


module load R/4.1.2
Rscript DE_ptscombined_log.R

exit
