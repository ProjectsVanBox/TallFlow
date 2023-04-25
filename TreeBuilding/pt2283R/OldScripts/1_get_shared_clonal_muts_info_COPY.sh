#!/bin/bash

VCF=$1

grep -v "^#" ${VCF} | \
grep -v ";CLONAL_SAMPLES=1[;,]" | \
grep -oP ";CLONAL_SAMPLE_NAMES=[\w,-]+"
