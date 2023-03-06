#!/usr/bin/env bash

#SBATCH --time=48:0:0
#SBATCH --mem=30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=r.hagelaars@prinsesmaximacentrum.nl
#SBATCH -o /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2229/PTATO.out
#SBATCH -e /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2229/PTATO.error


run_config_file=$1
projectDir=$(dirname "$0")

module load Java
module load R/4.1.2
export R_LIBS="/hpc/local/CentOS7/pmc_vanboxtel/R_libs/4.1.2":$R_LIBS


#source ${projectDir}/scripts/bash/params_check.sh

#check_params ${projectDir} ${run_config_file}

/hpc/pmc_vanboxtel/tools/external/nextflow_v21.10.6-all/nextflow run \
/hpc/pmc_vanboxtel/tools/ToolsVanBox_DEV/PTATO/ptato.nf \
-c ${run_config_file} \
--out_dir /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2229 \
-profile slurm -resume
