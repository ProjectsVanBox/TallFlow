#!/bin/bash

VCF=$1
CONTROL=$2
SAMPLE=$3

VCF_TREE=${VCF/.vcf/_tree_unique_$SAMPLE.vcf}

zgrep "^#" $VCF > $VCF_TREE

zgrep -P "SUBCLONAL_SAMPLES=[01]" $VCF | \
zgrep -P "SUBCLONAL_SAMPLE_NAMES=($CONTROL)*[;\s]" | \
zgrep -v -P "ABSENT_SAMPLES=0" | \
zgrep -P ";CLONAL_SAMPLES=[12];" | \
zgrep -P ";CLONAL_SAMPLE_NAMES=($CONTROL,)*$SAMPLE(,$CONTROL)*[;\s]" | \
zgrep -P "FAIL_QC_SAMPLES=[01]" | \
zgrep -P "FAIL_QC_SAMPLE_NAMES=($CONTROL)*[;\s]" \
>> $VCF_TREE


VCF_TREE=${VCF/.vcf/Root.vcf}
zgrep "^#" $VCF > $VCF_TREE
zgrep -P "ABSENT_SAMPLES=[01]" $VCF | \
>> $VCF_TREE


#grep -vP ";CLONAL_SAMPLE_NAMES=[\w,]*($CONTROL)[\w,]*" | \

#awk '{ if ( length($4) == 1 && length($5) == 1 ) { print $0 } }' \
