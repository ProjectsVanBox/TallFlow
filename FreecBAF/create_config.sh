#!/bin/bash

PILEUP=$1
CONFIG=/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/FreecBAF/config.txt

while read LINE; do
    if [[ ${LINE} = "mateFile"* ]]; then
      echo "mateFile = ${PILEUP}"
    else
       echo ${LINE} 
    fi

done < ${CONFIG}