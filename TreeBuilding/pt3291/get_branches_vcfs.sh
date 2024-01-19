#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt3291/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt3291/SMuRF.filtered.coords.vcf

cat ../header.vcf > TALL9_TALL6_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[1234578]" \
>> TALL9_TALL6_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[69]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[123457]" \
>> TALL9_TALL6_TALL8_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_TALL11_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL11" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[698][\w,-]+TALL[698]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([23457]|10|12|1[,;])" \
>> TALL9_TALL6_TALL8_TALL11_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_TALL11_TALL3_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([698]|11)[\w,-]+TALL([698]|11)[\w,-]+TALL([698]|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2457]|10|12|1[,;])" \
>> TALL9_TALL6_TALL8_TALL11_TALL3_branch.vcf

cat ../header.vcf > TALL12_TALL10_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL10" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([23456789]|11|1[,;])" \
>> TALL12_TALL10_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(10|12)" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([2457]|1[,;])" \
>> TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL5" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(10|12)" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([247]|1[,;])" \
>> TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(5|10|12)[\w,-]+TALL(5|10|12)" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([24]|1[,;])" \
>> TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_TALL2_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL2" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([57]|10|12)[\w,-]+TALL([57]|10|12)[\w,-]+TALL([57]|10|12)" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([4]|1[,;])" \
>> TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_TALL2_branch.vcf

cat ../header.vcf > TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_TALL2_TALL1_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL1[,;]" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([257]|10|12)[\w,-]+TALL([257]|10|12)[\w,-]+TALL([257]|10|12)[\w,-]+TALL([257]|10|12)" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)[\w,-]+TALL([3698]|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL4" \
>> TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_TALL2_TALL1_branch.vcf

cat ../header.vcf > root_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL4" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[12356789][\w,-]+TALL[12356789][\w,-]+TALL[12356789][\w,-]+TALL[12356789][\w,-]+TALL[12356789][\w,-]+TALL[12356789][\w,-]+TALL[12356789][\w,-]+TALL[12356789][\w,-]+TALL[12356789]" \
>> root_branch.vcf


