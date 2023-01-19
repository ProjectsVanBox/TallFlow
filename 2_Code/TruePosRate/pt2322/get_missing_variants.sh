#!/bin/bash

SAMPLE=$1
TP_BULK_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/intermediate/short_variants/somatic_vcfs/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-ALLBULK.SMuRF.filtered.vcf.gz
TP_SAMPLE_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/intermediate/short_variants/somatic_vcfs/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}.SMuRF.filtered.vcf.gz
SNV_PTATO_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/snvs/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}.ptato.filtered.vcf.gz
INDEL_PTATO_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/indels/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}/220919_HMFreg1764_pt2322_pt2322-DX1BM-${SAMPLE}.indels.ptato.filtered.vcf.gz

OUT_VCF=$(basename ${SNV_PTATO_VCF})
OUT_VCF=${OUT_VCF/.vcf.gz/}

zgrep -f <(zgrep -v "^#" ${TP_SAMPLE_VCF} | cut -f 1,2) <(zgrep -P "^\d+" ${TP_BULK_VCF}) > ${OUT_VCF}_TP.vcf

zgrep -f <(zgrep -v "^#" ${INDEL_PTATO_VCF} | cut -f1,2) <( zgrep -v -f <(zgrep -v "^#" ${SNV_PTATO_VCF} | cut -f 1,2) <(zgrep -P "^\d+" ${OUT_VCF}_TP.vcf )) > ${OUT_VCF}_MISSING.vcf
