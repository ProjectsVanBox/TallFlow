#!/bin/bash

for VCF in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones/*RSSmotif.vcf; do 
  N=$(grep -v "^#" ${VCF} | wc -l )
  N_MOTIFS=$(grep -P "[HN]M.*?[12]=" ${VCF} | wc -l)
  RANDOM_MOTIFS=$(shuf -n${N} /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones/random_pos_RSSmotif.txt | grep ":" | wc -l)


  PT=$(echo ${VCF} | grep -oP "pt\d+")
  CLONE=$(echo ${VCF} | grep -oP "TALL\d+")
  echo ${PT} ${CLONE} ${N} ${N_MOTIFS} ${RANDOM_MOTIFS}

done