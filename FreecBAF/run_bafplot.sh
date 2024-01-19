#!/bin/bash

module load R

BAF=$1
PREFIX=${BAF/_BAF.txt/}

R --slave --file=/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/FreecBAF/plotbaf.R --args ${BAF} 1000 ${PREFIX}