#!/bin/bash

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Drivers/pt2229

SMURF_MOD_HIGH_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Drivers/pt2229/SMuRF_moderate_high.vcf
SMURF_MOD_HIGH_TREE_VCF=SMuRF_moderate_high_tree.vcf

grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK[\w,-]*" ${SMURF_MOD_HIGH_VCF} | \
grep -P ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL[126][\w,-]+" | \
grep -vP ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL[345798][\w,-]*" \
> ${SMURF_MOD_HIGH_TREE_VCF}

grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK[\w,-]*" ${SMURF_MOD_HIGH_VCF} | \
grep -P ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL[345798][\w,-]+" | \
grep -vP ";CLONAL_SAMPLE_NAMES=[\w,-]+TALL[26][\w,-]*" \
>> ${SMURF_MOD_HIGH_TREE_VCF}

grep -P "CLONAL_SAMPLE_NAMES=[\w,-]+BULK[\w,-]*" ${SMURF_MOD_HIGH_VCF} | \
grep -P "ABSENT_SAMPLES=[01234]" \
>> ${SMURF_MOD_HIGH_TREE_VCF}

SMURF_MOD_HIGH_TREE_VCF_SORTED=${SMURF_MOD_HIGH_TREE_VCF/.vcf/_sorted.vcf}
sort -k1n,2n ${SMURF_MOD_HIGH_TREE_VCF} | uniq > ${SMURF_MOD_HIGH_TREE_VCF_SORTED}

MOD_HIGH_SMURF_VCF=moderate_high_SMuRF.vcf
MOD_HIGH_SMURF_TREE_VCF=moderate_high_SMuRF_tree_sorted.vcf
grep -f <( cut -f 1,2 ${SMURF_MOD_HIGH_TREE_VCF} ) ${MOD_HIGH_SMURF_VCF} | sort -k1n,2n > ${MOD_HIGH_SMURF_TREE_VCF}

DRIVERS_TREE_FILE=drivers_tree.txt
paste \
<( cut -f 1-5 ${MOD_HIGH_SMURF_TREE_VCF} ) \
<( grep -oP "ANN=.+?ENSG" ${MOD_HIGH_SMURF_TREE_VCF} | cut -f 4 -d'|' ) \
<( grep -oP ";CLONAL_SAMPLE_NAMES=[\w,-]+" ${SMURF_MOD_HIGH_TREE_VCF_SORTED} ) \
<( grep -oP ";SUBCLONAL_SAMPLE_NAMES=[\w,-]*" ${SMURF_MOD_HIGH_TREE_VCF_SORTED} ) \
<( grep -oP ";ABSENT_SAMPLE_NAMES=[\w,-]*" ${SMURF_MOD_HIGH_TREE_VCF_SORTED} ) \
> ${DRIVERS_TREE_FILE}

