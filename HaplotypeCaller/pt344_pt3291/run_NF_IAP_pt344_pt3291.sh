#!/bin/bash

#SBATCH --time=10:0:0
#SBATCH --mem=25G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=r.hagelaar@prinsesmaximacentrum.nl

CONFIG=run_NF_IAP_pt344_pt3291.config

module load Java

/hpc/pmc_vanboxtel/tools/external/nextflow_v20.03.0-edge.5288/nextflow run /hpc/pmc_vanboxtel/tools/ToolsVanBox/NF-IAP_v1.3.0/nf-iap.nf \
 -c ${CONFIG} --out_dir /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/HaplotypeCaller/pt344_pt3291/ -profile slurm -resume
