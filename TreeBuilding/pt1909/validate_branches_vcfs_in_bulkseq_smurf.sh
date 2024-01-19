#!/bin/bash

SMURF_FILT_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt1909/231017_HMFreg2090_pt1909_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf
SMURF_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt1909/231017_HMFreg2090_pt1909_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.vcf

BRANCHES_VCF=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1909/branches_vcfs/branches_merged.vcf

OUTDIR=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1909/branches_vcfs/overlapBulkSeq

#grep -f <(grep -v "^#" ${BRANCHES_VCF} | cut -f 1,2) ${SMURF_FILT_VCF} > ${OUTDIR}/bulkseqSMuRFfiltVCF_overlap_branches.vcf

#grep -f <(grep -v "^#" ${SMURF_FILT_VCF} | cut -f 1,2) ${BRANCHES_VCF} > ${OUTDIR}/branchesVCF_overlap_bulkSMuRFfilt.vcf

#grep -vf <(grep -v "^#" ${SMURF_FILT_VCF} | cut -f 1,2) ${BRANCHES_VCF} > ${OUTDIR}/branchesVCF_NOTIN_bulkSMuRFfilt.vcf

#grep -vf <(grep -v "^#" ${BRANCHES_VCF} | cut -f 1,2) ${SMURF_FILT_VCF} > ${OUTDIR}/bulkseqSMuRFfiltVCF_NOTIN_branches.vcf

for BRANCH_VCF in ${OUTDIR}/../{root,TALL}*.vcf; do 
    OUT_VCF=${BRANCH_VCF/.vcf/_BulkSeqSMuRF.vcf}
    OUT_VCF=$(basename ${OUT_VCF})
    grep "^#" ${SMURF_VCF} > ${OUTDIR}/${OUT_VCF}
    grep -f <(grep -v "^#" ${BRANCH_VCF} | cut -f 1,2) ${SMURF_VCF} >> ${OUTDIR}/${OUT_VCF}
done
    
for BRANCH_VCF in ${OUTDIR}/../{root,TALL}*.vcf; do 
    BRANCH_BULKSEQ_VCF=${BRANCH_VCF/.vcf/_BulkSeqSMuRF.vcf}
    BRANCH_BULKSEQ_VCF=${BRANCH_BULKSEQ_VCF/../}
    OUT_VCF=${BRANCH_BULKSEQ_VCF/.vcf/_MISSED.vcf}
    OUT_VCF=$(basename ${OUT_VCF})

    echo "grep -vf <(grep -v \"^#\" ${BRANCH_BULKSEQ_VCF} | cut -f 1,2) ${BRANCH_VCF} > ${OUTDIR}/${OUT_VCF}"
    grep -vf <(grep -v "^#" ${BRANCH_BULKSEQ_VCF} | cut -f 1,2) ${BRANCH_VCF} > ${OUTDIR}/${OUT_VCF}
done
