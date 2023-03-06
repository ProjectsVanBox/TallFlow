#!/bin/bash

DIR=$1
FLANK=$2

if [ -z ${FLANK} ]; then
    FLANK=250
fi

cat /dev/null > tmp.bed

for VCF in ${DIR}/*/*.filtered.vcf.gz; do 

zgrep -v "^#" ${VCF} | awk -v flank=${FLANK} '{ print $1"\t"$2-flank"\t"$2+flank }' >> tmp.bed

done

bedtools sort -i tmp.bed | bedtools merge > baitset.bed
rm tmp.bed