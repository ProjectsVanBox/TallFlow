#!/bin/bash

get_callable_size() {
    SAMPLE_FNAME=$1
    CONTROL_FNAME=$2
    
    CMD="bedtools intersect -a <(grep \"CALLABLE\" ${SAMPLE_FNAME} | grep -vP \"^\D\") -b <(grep \"CALLABLE\" ${CONTROL_FNAME} | grep -vP \"^\D+\") | \
    bedtools sort -i stdin | \
    bedtools merge -i stdin | \
    awk 'BEGIN{ SUM=0 }{ SUM+=\$3-\$2 }END{ print SUM }'
    "
    
    SIZE=$(eval ${CMD})
    echo ${SIZE}
}

get_snv_load() {
    VCF=$1
    CMD="zgrep -v \"^#\" ${VCF} | awk '{ if (length(\$4) == 1 && length(\$5) == 1) { print \$0 } }' | wc -l"
    
    SNV_LOAD=$(eval ${CMD})
    echo ${SNV_LOAD}    
}

get_indel_load() {
    VCF=$1
    CMD="zgrep -v \"^#\" ${VCF} | awk '{ if (length(\$4) != 1 || length(\$5) != 1) { print \$0 } }' | wc -l"

    INDEL_LOAD=$(eval ${CMD})
    echo ${INDEL_LOAD}    
}

echo "DONOR SAMPLE CONTROL CALLABLE SNV_LOAD INDEL_LOAD"

DONOR="pt2229"
CONTROL_CL="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CallableLociBulk/pt2229/pt2229-DX1BM-MSCBULK_dedup_CallableLoci.bed"
CONTROL=$( echo ${CONTROL_CL} | grep -oP "[\w,-]+.CallableLoci.bed" | cut -f1 -d'_')
for SAMPLE_CL in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CallableLociBulk/pt2229/*CallableLoci.bed; do
    SAMPLE=$( echo ${SAMPLE_CL} | grep -oP "[\w,-]+.CallableLoci.bed" | cut -f1 -d'_')
    if [[ "${SAMPLE_CL}" == "${CONTROL_CL}" ]]; then # || "${SAMPLE_CL}" == *"BULK"* 
		continue
    fi
    CALLABLE_SIZE=$(get_callable_size ${SAMPLE_CL} ${CONTROL_CL})
    # Something funky with the vcf files 
    for SAMPLE_VCF in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2229inclDP/230124_HMFreg1891_pt2229.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_*${SAMPLE}.SMuRF.filtered.vcf; do
    	#echo $SAMPLE_VCF
		SNV_LOAD=$(get_snv_load ${SAMPLE_VCF})
		INDEL_LOAD=$(get_indel_load ${SAMPLE_VCF})
	
    done    
    
    echo ${DONOR} ${SAMPLE} ${CONTROL} ${CALLABLE_SIZE} ${SNV_LOAD} ${INDEL_LOAD}
done








