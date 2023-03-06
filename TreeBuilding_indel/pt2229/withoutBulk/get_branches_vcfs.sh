#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2229/withoutBulk/branches_vcfs

VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2229/SMuRF.filtered.coords.vcf

cat ../header.vcf > TALL6_TALL1_TALL2_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(2|1[;,])" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(10|[895347])" \
>> TALL6_TALL1_TALL2_branch.vcf

cat ../header.vcf > TALL1_TALL2_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL1[;,]" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL2" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL6" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(10|[895347])" \
>> TALL1_TALL2_branch.vcf

cat ../header.vcf > TALL10_TALL8_TALL9_TALL5_TALL3_TALL4_TALL7_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL7" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(10|[89])[\w,-]+TALL(10|[89])" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[543][\w,-]+TALL[543]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([62]|1[;,])" \
>> TALL10_TALL8_TALL9_TALL5_TALL3_TALL4_TALL7_branch.vcf

cat ../header.vcf > TALL10_TALL8_TALL9_TALL5_TALL3_TALL4_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL(10|[89])[\w,-]+TALL(10|[89])" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[543][\w,-]+TALL[543]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([627]|1[;,])" \
>> TALL10_TALL8_TALL9_TALL5_TALL3_TALL4_branch.vcf

cat ../header.vcf > TALL10_TALL8_TALL9_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL10" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[89]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([627534]|1[;,])" \
>> TALL10_TALL8_TALL9_branch.vcf

cat ../header.vcf > TALL8_TALL9_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL8" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL9" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([627534]|1[;,]|10)" \
>> TALL8_TALL9_branch.vcf

cat ../header.vcf > TALL5_TALL3_TALL4_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL5" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL[43]" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([89627]|1[;,]|10)" \
>> TALL5_TALL3_TALL4_branch.vcf

cat ../header.vcf > TALL3_TALL4_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL3" | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL4" | \
grep -vP "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([589627]|1[;,]|10)" \
>> TALL3_TALL4_branch.vcf

cat ../header.vcf > root_branch.vcf
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+" ${VCF} | \
grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+TALL([62954]|10)[\w,-]+TALL([62954]|10)[\w,-]+TALL([62954]|10)[\w,-]+TALL([62954]|10)[\w,-]+TALL([62954]|10)" \
>> root_branch.vcf

