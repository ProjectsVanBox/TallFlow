#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1180/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1180/SMuRF.filtered.coords.vcf

cat ../header.vcf > TALL11_TALL5_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL11" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL5" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2346789]|10|12|1[,;])" \
>> TALL11_TALL5_branch.vcf

cat ../header.vcf > TALL11_TALL5_TALL3_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[11|5]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([246789]|10|12|1[,;])" \
>> TALL11_TALL5_TALL3_branch.vcf

cat ../header.vcf > TALL11_TALL5_TALL3_TALL6_branch.vcf  ## present in TALL6, and in 2 out of the 3 others
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([35]|11)[\w,-]+TALL([35]|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([24789]|10|12|1[,;])" \
>> TALL11_TALL5_TALL3_TALL6_branch.vcf

cat ../header.vcf > TALL4_TALL12_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL4" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([12356789]|10|11|1[,;])" \
>> TALL4_TALL12_branch.vcf

cat ../header.vcf > TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_branch.vcf  ## present in n-1 of each subbranch
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([356]|11)[\w,-]+TALL([356]|11)[\w,-]+TALL([356]|11)" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(4|12)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2789]|10|1[,;])" \
>> TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_branch.vcf

cat ../header.vcf > TALL10_TALL1_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL1[,;]" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL10" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([23456789]|11|12)" \
>> TALL10_TALL1_branch.vcf

cat ../header.vcf > TALL10_TALL1_TALL2_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL2" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(10|1[,;])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3456789]|11|12)" \
>> TALL10_TALL1_TALL2_branch.vcf

cat ../header.vcf > TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_TALL10_TALL1_TALL2_branch.vcf  ## present in n-1 of each subbranch
cat ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3456]|11|12)[\w,-]+TALL([3456]|11|12)[\w,-]+TALL([3456]|11|12)" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(2|10|1[,;])[\w,-]+TALL(2|10|1[,;])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[789]" \
>> TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_TALL10_TALL1_TALL2_branch.vcf

cat ../header.vcf > TALL9_TALL7_branch.vcf ## exclude all TALL starting with 1. 
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[1234568]" \
>> TALL9_TALL7_branch.vcf

cat ../header.vcf > TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_TALL10_TALL1_TALL2_TALL9_TALL7_branch.vcf ## present in n-1 of each subbranch
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([79])" ${VCF}  | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3456]|11|12)[\w,-]+TALL([3456]|11|12)[\w,-]+TALL([3456]|11|12)" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(2|10|1[,;])[\w,-]+TALL(2|10|1[,;])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" \
>> TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_TALL10_TALL1_TALL2_TALL9_TALL7_branch.vcf

cat ../header.vcf > root_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[8]" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679]" \
>> root_branch.vcf

