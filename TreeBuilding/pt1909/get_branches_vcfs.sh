#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1909/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1909/SMuRF.filtered.coords.vcf

cat ../header.vcf > TALL1_TALL5_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL1[,;]" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL5" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[2346789]" \
>> TALL1_TALL5_branch.vcf

cat ../header.vcf > TALL1_TALL5_TALL4_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL4" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(5|1[,;])]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[236789]" \
>> TALL1_TALL5_TALL4_branch.vcf

cat ../header.vcf > TALL1_TALL5_TALL4_TALL9_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([45]|1[,;])[\w,-]+TALL([45]|1[,;])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[23678]" \
>> TALL1_TALL5_TALL4_TALL9_branch.vcf

cat ../header.vcf > TALL1_TALL5_TALL4_TALL9_TALL7_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([459]|1[,;])[\w,-]+TALL([459]|1[,;])[\w,-]+TALL([459]|1[,;])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[2368]" \
>> TALL1_TALL5_TALL4_TALL9_TALL7_branch.vcf

cat ../header.vcf > TALL6_TALL2_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL2" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([345789]|1[,;])" \
>> TALL6_TALL2_branch.vcf

cat ../header.vcf > TALL6_TALL2_TALL3_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[26]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([45789]|1[,;])" \
>> TALL6_TALL2_TALL3_branch.vcf

cat ../header.vcf > TALL1_TALL5_TALL4_TALL9_TALL7_TALL6_TALL2_TALL3_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[236][\w,-]+TALL[236]" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([4579]|1[,;])[\w,-]+TALL([4579]|1[,;])[\w,-]+TALL([4579]|1[,;])[\w,-]+TALL([4579]|1[,;])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[8]" \
>> TALL1_TALL5_TALL4_TALL9_TALL7_TALL6_TALL2_TALL3_branch.vcf

cat ../header.vcf > root_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[8]" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679][\w,-]+TALL[12345679]" \
>> root_branch.vcf