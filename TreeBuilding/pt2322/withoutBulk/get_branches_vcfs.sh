#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/withoutBulk/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/SMuRF.filtered.coords.vcf

cat ../header.vcf > TALL6_TALL4_TALL5_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[45]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[12379]" \
>> TALL6_TALL4_TALL5_branch.vcf

cat ../header.vcf > TALL4_TALL5_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL4" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL5" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[123879]" \
>> TALL4_TALL5_branch.vcf

cat ../header.vcf > TALL6_TALL4_TALL5_TALL8_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[456][\w,-]+TALL[456]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[12379]" \
>> TALL6_TALL4_TALL5_TALL8_branch.vcf

cat ../header.vcf > TALL10_TALL9_TALL12_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL10" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(9|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2345678]|1[;,])" \
>> TALL10_TALL9_TALL12_branch.vcf

cat ../header.vcf > root_branch.vcf
cat ${VCF} | \
grep -P "ABSENT_SAMPLES=([02]|1[;,])" \
>> root_branch.vcf

cat ../header.vcf > TALL9_TALL11_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL11" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2345678]|1[,;])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL10" \
>> TALL9_TALL11_branch.vcf

cat ../header.vcf > TALL7_TALL12_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2345689]|10|11|1[,;])" \
>> TALL7_TALL12_branch.vcf

cat ../header.vcf > TALL1_TALL3_TALL2_TALL7_TALL12_TALL10_TALL9_TALL11_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([237]|12|1[,;])[\w,-]+TALL([237]|12|1[,;])[\w,-]+TALL([237]|12|1[,;])[\w,-]+TALL([237]|12|1[,;])" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(9|10|11)[\w,-]+TALL(9|10|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[4568]" \
>> TALL1_TALL3_TALL2_TALL7_TALL12_TALL10_TALL9_TALL11_branch.vcf

cat ../header.vcf > TALL1_TALL3_TALL2_TALL7_TALL12_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(3|1[,;])" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(2|7|12)[\w,-]+TALL(2|7|12)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([45689]|10|11)" \
>> TALL1_TALL3_TALL2_TALL7_TALL12_branch.vcf

cat ../header.vcf > TALL1_TALL3_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL1[;,]" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2456789]|10|11|12)" \
>> TALL1_TALL3_branch.vcf

cat ../header.vcf > TALL2_TALL7_TALL12_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL2" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(7|12)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([345689]|10|11|1[,;])" \
>> TALL2_TALL7_TALL12_branch.vcf

cat ../header.vcf > TALL7_TALL12_branch.vcf
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2345689]|10|11|1[,;])" \
>> TALL7_TALL12_branch.vcf


