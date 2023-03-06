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

for DONOR in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/snvs/*; do
    DONOR=$(basename ${DONOR} )
    for CONTROL_CL in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/QC/CallableLoci/${DONOR}/*MONO*callableloci.bed; do
	CONTROL=$( echo ${CONTROL_CL} | grep -oP "[\w,-]+.callableloci.bed" | cut -f1 -d'.')
	for SAMPLE_CL in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/QC/CallableLoci/${DONOR}/*callableloci.bed; do
	    SAMPLE=$( echo ${SAMPLE_CL} | grep -oP "[\w,-]+.callableloci.bed" | cut -f1 -d'.')
	    if [[ "${SAMPLE_CL}" == "${CONTROL_CL}" || "${SAMPLE_CL}" == *"BULK"* ]]; then
		continue
	    fi
	    CALLABLE_SIZE=$(get_callable_size ${SAMPLE_CL} ${CONTROL_CL})
	    
	    for SAMPLE_VCF in /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/snvs/${DONOR}/220919_HMFreg1764_pt2322_${SAMPLE}/*${SAMPLE}.ptato.filtered.vcf.gz; do
		SNV_LOAD=$(get_snv_load ${SAMPLE_VCF})
		INDEL_LOAD=$(get_indel_load ${SAMPLE_VCF})
		
	    done    
	    
	    echo ${DONOR} ${SAMPLE} ${CONTROL} ${CALLABLE_SIZE} ${SNV_LOAD} ${INDEL_LOAD}
	done
    done
done