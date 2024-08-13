#!/bin/bash

#SBATCH --job-name=ModuleScoringCordes
#SBATCH --mem=15G
#SBATCH --time=0:30:0
#SBATCH --mail-type=START,FAIL,END
#SBATCH --mail-user=v.m.poort@prinsesmaximacentrum.nl


module load R/4.1.2
Rscript ModuleScoringVis.R

exit