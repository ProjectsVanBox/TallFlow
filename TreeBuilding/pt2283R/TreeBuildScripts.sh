#!/bin/bash


#VCF=$1
VCF="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2283R/230307_HMF1931_pt2283R.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf"
VCF_TREE="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283R/oldscript.vcf"
VCF_TREE_simple="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283R/simple.vcf"
VCF_TREE_NoCounts="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283R/NoCounts.vcf"
VCF_TREE_ExcludeInAll="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283R/ExcludeInAll.vcf"
mytreedir="/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283R/"

#grep "^#" $VCF > $VCF_TREE_simple
#grep -v "^#" ${VCF} | \
#grep -v ";CLONAL_SAMPLES=1[;,]" | \
#grep -oP ";CLONAL_SAMPLE_NAMES=[\w,-]+" >> $VCF_TREE_simple

#grep "^#" $VCF > $VCF_TREE_NoCounts
#grep -v "^#" ${VCF} | \
#grep -v ";CLONAL_SAMPLES=1[;,]"  >> $VCF_TREE_NoCounts


grep "^#" $VCF > $VCF_TREE_ExcludeInAll
grep -v "^#" ${VCF} | \
grep -v ";CLONAL_SAMPLES=1[;,]" | \
grep -v ";CLONAL_SAMPLES=5" >> $VCF_TREE_ExcludeInAll

#grep "^#" $VCF > $VCF_TREE
#grep -P "SUBCLONAL_SAMPLES=[01]" $VCF | \
#grep -v -P "ABSENT_SAMPLES=0" | \
#grep -v -P ";CLONAL_SAMPLES=[01];" | \
#grep -P "FAIL_QC_SAMPLES=[01]" | \
#>> $VCF_TREE




#cat ${VCF} | \
#grep -v "^#" | \
#cut -f 1,2 | \
#sort | \
#uniq > ${mytreedir}/shared_clonal_muts_coords.txt


#cat ${VCF} | grep -v "^#" | cut -f 1,2 | sort | uniq > ${mytreedir}/shared_clonal_muts_coords.txt

#grep "^#" $VCF > ${mytreedir}/SMuRF.filtered.coords.vcf
#grep -f ${mytreedir}/SMuRF.filtered.coords.vcf $VCF >> ${mytreedir}/SMuRF.filtered.coords.vcf