#!/bin/bash
#SBATCH --job-name=SMuRF_pt3291
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --time=12:0:0
#SBATCH -e SMuRF_pt3291.err
#SBATCH -o SMuRF_pt3291.log

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt3291

/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip -c /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/VCFS/VCF/230918_HMFreg2063_pt3291_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf > /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/VCFS/VCF/230918_HMFreg2063_pt3291_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz
/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/tabix /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/VCFS/VCF/230918_HMFreg2063_pt3291_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz

singularity exec -B /hpc /hpc/local/CentOS7/pmc_vanboxtel/singularity_cache/vanboxtelbioinformatics-smurf-3.0.2.img python /smurf/SMuRF.py \
  -i /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/VCFS/VCF/230918_HMFreg2063_pt3291_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz \
  -b /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/BAMS/pt3291-DX1PB-ALLBULK_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/BAMS/pt3291-DX1PB-ALLBULK-DNCD1a-min_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/BAMS/pt3291-DX1PB-ALLBULK-DNCD1a-plus_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/BAMS/pt3291-DX1PB-ALLBULK-DP_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/230918_HMFreg2063_pt3291_Bulk/BAMS/pt3291-DX1PB-ALLBULK-SPCD8_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/230612_HMFreg2010_pt344/BAMS/pt344-DX1PB-MONOBULK_dedup.bam \
  -n pt344-DX1PB-MONOBULK \
  -c /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/config_15cov.ini


 