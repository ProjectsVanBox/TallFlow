#!/bin/bash

BED1=$1
BED2=$2

cat ${BED1} ${BED2} | bedtools sort | bedtools merge > baitsets_merged.bed