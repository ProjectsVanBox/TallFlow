#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2283D/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2283D/SMuRF.filtered.coords.vcf

cat ../header.vcf > root_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" | \
grep -P "ABSENT_SAMPLES=([02]|1[;,])" \
>> root_branch.vcf

cat ../header.vcf > TALL1_TALL2_TALL4_TALL5_TALL6_TALL7_TALL8_TALL9_TALL10_TALL11_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)[\w,-]+TALL([245678]|1|10|11)" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL12" \
>> TALL1_TALL2_TALL4_TALL5_TALL6_TALL7_TALL8_TALL9_TALL10_TALL11_branch.vcf

