#!/bin/bash

VCF=$1

# _1  STARTED WITH -o 200 -f 500
# _2  STARTED WITH -o 100 -f 500
# _3  STARTED WITH -o 150 -f 500

bash /hpc/pmc_vanboxtel/tools/ToolsVanBox/PrimerDesign_v.1.0.0/primer_design_pipeline.sh -v ${VCF} -g 38 -o 200 -f 500