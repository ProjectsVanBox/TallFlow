#!/bin/bash

VCF=$1
OUT=$(basename ${VCF})

OUT=${OUT/%.vcf/.drivers.vcf}

grep -P "^#|MODERATE|HIGH" ${VCF} > /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Drivers/pt1180_pt1909/${OUT}
