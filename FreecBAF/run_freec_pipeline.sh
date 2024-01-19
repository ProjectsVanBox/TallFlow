#!/bin/bash

#SBATCH --time=4:0:0
#SBATCH --mem=20G
#SBATCH -c 6

BAM=$1
FASTA=/hpc/pmc_vanboxtel/projects/TallFlow/1_Input/FreecBAF/freec/Homo_sapiens_assembly38.fasta
PILEUP=${BAM/.bam/.mpileup}
PILEUP=$( basename ${PILEUP} )

#bash /hpc/pmc_vanboxtel/projects/TallFlow/2_Code/FreecBAF/run_pileup.sh ${BAM}

CONFIG=${PILEUP/.mpileup/_config.txt}

bash /hpc/pmc_vanboxtel/projects/TallFlow/2_Code/FreecBAF/create_config.sh ${PILEUP} > ${CONFIG}

bash /hpc/pmc_vanboxtel/projects/TallFlow/2_Code/FreecBAF/run_freec.sh ${CONFIG}


