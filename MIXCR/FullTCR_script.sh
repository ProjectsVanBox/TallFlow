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


## Create directories 
mkdir -p $OutDir/input/bam
mkdir -p $OutDir/input/fastq/Merged
mkdir -p $OutDir/Output/$SampleName 


## Create a bamslice 
samtools view -L $TCR -bo $OutDir/input/bam/${SampleName}_dedup.bam $BAM


## Create fastq files 
cd $OutDir/input/fastq
java -Xmx60G -jar $PICARD SamToFastq \
    INPUT=$OutDir/input/bam/${SampleName}_dedup.bam \
    OUTPUT_DIR=$PWD \
    RG_TAG=ID \
    OUTPUT_PER_RG=true
## Merge fastq files 
cat $OutDir/input/fastq/${SampleName}*_1.fastq > $OutDir/input/fastq/Merged/${SampleName}_R1.fastq
cat $OutDir/input/fastq/${SampleName}*_2.fastq > $OutDir/input/fastq/Merged/${SampleName}_R2.fastq
cd -


### Check MIXCR for the TCR receptor
#/hpc/pmc_vanboxtel/tools/mixcr/mixcr activate-license

/hpc/pmc_vanboxtel/tools/mixcr/mixcr analyze shotgun --species hsa --starting-material dna \
    $OutDir/input/fastq/Merged/${SampleName}_R1.fastq \
    $OutDir/input/fastq/Merged/${SampleName}_R2.fastq \
    $OutDir/Output/$SampleName/${SampleName}_TCR






