#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt344/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt344/SMuRF.filtered.coords.vcf

cat ../header.vcf > TALL1_TALL5_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL1" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL5" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[2346789]" \
>> TALL1_TALL5_branch.vcf

cat ../header.vcf > TALL1_TALL5_TALL2_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL2" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[15]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[346789]" \
>> TALL1_TALL5_TALL2_branch.vcf

cat ../header.vcf > TALL8_TALL7_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[1234569]" \
>> TALL8_TALL7_branch.vcf

cat ../header.vcf > TALL8_TALL7_TALL9_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[87]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[123456]" \
>> TALL8_TALL7_TALL9_branch.vcf

cat ../header.vcf > TALL1_TALL5_TALL2_TALL8_TALL7_TALL9_TALL3_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[152879][\w,-]+TALL[152879][\w,-]+TALL[152879][\w,-]+TALL[152879][\w,-]+TALL[152879]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[46]" \
>> TALL1_TALL5_TALL2_TALL8_TALL7_TALL9_TALL3_branch.vcf

cat ../header.vcf > TALL1_TALL5_TALL2_TALL8_TALL7_TALL9_TALL3_TALL6_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[1528793][\w,-]+TALL[1528793][\w,-]+TALL[1528793][\w,-]+TALL[1528793][\w,-]+TALL[1528793][\w,-]+TALL[1528793]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[4]" \
>> TALL1_TALL5_TALL2_TALL8_TALL7_TALL9_TALL3_TALL6_branch.vcf



