#!/bin/bash

#SBATCH -t 48:0:0
#SBATCH --mem=40G
#SBATCH -c 16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=r.hagelaar@prinsesmaximacentrum.nl
#SBATCH -o /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/pt2322-DX1BM-ALLBULK-SPCD4/logs/pt2322-DX1BM-ALLBULK-SPCD4.log

## Script standardly runs GRIDSS+GRIPSS in tumor-normal and the rest in tumor-only mode

GRIDSS_PURPLE_LINX_SH=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/gridss-purple-linx/gridss-purple-linx_v3.sh
GUIX_PATH=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/.guix-profile

GRIDSS=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/gridss/gridss.jar
AMBER=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/hmftools/amber-3.5.jar
GRIPSS=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/hmftools/gripss-1.9.jar
COBALT=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/hmftools/cobalt-1.11.jar
PURPLE=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/hmftools/purple-2.51.jar
LINX=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/hmftools/sv-linx_v1.12.jar

TUMOR_BAM="/hpc/pmc_vanboxtel/processed/221026_HMFreg1816_pt2322/BAMS/pt2322-DX1BM-ALLBULK-SPCD4_dedup.bam"
OUTPUTDIR="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/pt2322-DX1BM-ALLBULK-SPCD4/"
NORMAL_BAM="/hpc/pmc_vanboxtel/processed/220919_HMFreg1764_pt2322/BAMS/pt2322-DX1BM-MONOBULK_dedup.bam"
SAMPLE_NAME=$(basename $TUMOR_BAM)
SAMPLE_NAME=${SAMPLE_NAME%%_*}
NORMAL_NAME=$(basename $NORMAL_BAM)
NORMAL_NAME=${NORMAL_NAME%%_*}
THREADS=16

## Purple cannot read the SMuRF vcfs at this point, even after removing all other samples and running the AnnotateStrelkaWithAllelicDepth script
#SOMATICSNV=""

echo $SAMPLE_NAME
#mkdir $SAMPLE_NAME
RUN_DIR=${OUTPUTDIR}
echo $RUN_DIR
#cd $RUN_DIR
echo date

guixr load-profile ${GUIX_PATH} -- << EOF
export GRIDSS_JAR=${GRIDSS}
export AMBER_JAR=${AMBER}
export COBALT_JAR=${COBALT}
export GRIPSS_JAR=${GRIPSS}
export PURPLE_JAR=${PURPLE}
export LINX_JAR=${LINX}
bash ${GRIDSS_PURPLE_LINX_SH} \
  -t ${TUMOR_BAM} \
  -n ${NORMAL_BAM} \
  -s ${SAMPLE_NAME} \
  --rundir $RUN_DIR \
  -o $RUN_DIR \
  --nosnvvcf \
  --threads ${THREADS}
EOF
