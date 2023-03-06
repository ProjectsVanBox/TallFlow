#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/SMuRF.filtered.coords.vcf

cat ../header.vcf > ISPCD4_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[45]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[12379]" \
>> ISPCD4_branch.vcf


cat ../header.vcf > ISPCD4_TALL8_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[456]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[12379]" \
>> ISPCD4_TALL8_branch.vcf


cat ../header.vcf > SPCD4_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL10" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(9|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[2345678]" \
>> SPCD4_branch.vcf

cat ../header.vcf > root_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK" ${VCF} | \
grep -P "ABSENT_SAMPLES=[012]" \
>> root_branch.vcf

cat ../header.vcf > TALL9_TALL11_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL11" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[2345678]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL10" \
>> TALL9_TALL11_branch.vcf

cat ../header.vcf > TALL7_TALL12_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[2345689]" \
>> TALL7_TALL12_branch.vcf

