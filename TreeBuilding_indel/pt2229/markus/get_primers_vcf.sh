#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/primers

TALL7_TALL12_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/branches_vcfs/TALL7_TALL12_branch.vcf
TALL9_TALL11_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/branches_vcfs/TALL9_TALL11_branch.vcf

cat ../header.vcf > TALL7_TALL12_primer.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK" ${TALL7_TALL12_VCF} | \
grep -P ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" ${TALL7_TALL12_VCF} | \
grep -P ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" ${TALL7_TALL12_VCF} | \
grep -P "ABSENT_SAMPLES=10" \
>> TALL7_TALL12_primer.vcf

cat ../header.vcf > TALL9_TALL11_primer.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK" ${TALL9_TALL11_VCF} | \
grep -P ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" ${TALL9_TALL11_VCF} | \
grep -P ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL11" ${TALL9_TALL11_VCF} | \
grep -P "ABSENT_SAMPLES=10" \
>> TALL9_TALL11_primer.vcf
