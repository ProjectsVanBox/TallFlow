#!/bin/bash

BED=$1
PT=pt2322

for VCF in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/*/intermediate/svs/gridss/*/*/${PT}*unfiltered.vcf.gz; do
    OUTVCF=$( basename ${VCF} )
    OUTVCF=${OUTVCF/%.vcf.gz/_BULK.vcf}
    zgrep "^#" ${VCF} > ${OUTVCF}
    zgrep -f <(cut -f 1,3 ${BED}) ${VCF} >> ${OUTVCF}
done