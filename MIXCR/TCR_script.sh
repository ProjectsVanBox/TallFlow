#!/bin/bash

#SBATCH -t 10:0:0
#SBATCH --mem=25G
#SBATCH -c 8
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=m.j.vanroosmalen-3@prinsesmaximacentrum.nl


### Load modules 
module load samtools
module load Java/1.8.0_60


### Variables 
BAM=${1}
OutDir=${2}
TCR="/hpc/pmc_vanboxtel/personal/rhagelaar/TCR_MIXCR/TCR_subset.bed"
PICARD=/hpc/local/CentOS7/pmc_vanboxtel/bin/picard-2.24.1/picard.jar
BGZIP=/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip
PWD=`pwd`/
SampleName=$(basename $BAM | sed 's|_dedup.bam||g')


### Create directories 
mkdir -p $OutDir/input/bam
mkdir -p $OutDir/input/fastq/Merged
mkdir -p $OutDir/Output/$SampleName 


### Check MIXCR for the TCR receptor
/hpc/pmc_vanboxtel/tools/mixcr/mixcr analyze shotgun --species hsa --starting-material dna \
    $OutDir/input/fastq/Merged/${SampleName}_R1.fastq \
    $OutDir/input/fastq/Merged/${SampleName}_R2.fastq \
    $OutDir/Output/$SampleName/${SampleName}_TCR






