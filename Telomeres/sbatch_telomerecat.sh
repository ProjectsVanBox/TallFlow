#!/bin/bash

#SBATCH --job-name=Telomerecat
#SBATCH --time=2:0:0
#SBATCH -c 8

. /hpc/pmc_vanboxtel/tools/external/telomerecat_v4.0.0/venv_3.6/bin/activate

BAM=$1

echo $BAM

telomerecat bam2length -p8 --temp_dir ./ $BAM

echo Done