#!/bin/bash

CONF=$1

singularity exec -B /hpc /hpc/local/Rocky8/pmc_vanboxtel/singularity_cache/depot.galaxyproject.org-singularity-control-freec-11.6--h1b792b2_1.img \
freec -conf ${CONF}
