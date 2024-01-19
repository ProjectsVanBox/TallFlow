#!/bin/bash

VCF=$1
OUT=$(basename ${VCF})

OUT=${OUT/%.vcf/.sv.vcf}

grep -P "^#|MODERATE|HIGH" ${VCF} > /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Drivers/pt344/${OUT}


grep -P "^11" ${VCF} | grep -P "(\[|\})14:" | grep -P "PASS" > /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/pt344/${OUT}