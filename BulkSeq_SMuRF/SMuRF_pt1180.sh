#!/bin/bash
#SBATCH --job-name=SMuRF_pt1180
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --time=12:0:0
#SBATCH -e SMuRF_pt1180.err
#SBATCH -o SMuRF_pt1180.log

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt1180

/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip -c /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/VCFS/VCF/231017_HMFreg2090_pt1180_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf > /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/VCFS/VCF/231017_HMFreg2090_pt1180_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz
/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/tabix /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/VCFS/VCF/231017_HMFreg2090_pt1180_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz

singularity exec -B /hpc /hpc/local/CentOS7/pmc_vanboxtel/singularity_cache/vanboxtelbioinformatics-smurf-3.0.2.img python /smurf/SMuRF.py \
  -i /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/VCFS/VCF/231017_HMFreg2090_pt1180_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz \
  -b /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/BAMS/pt1180-DX1BM-MSCBULK_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/BAMS/pt1180-DX1PB-ALLBULK_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/BAMS/pt1180-DX1PB-ALLBULK-DNCD1aNeg_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/BAMS/pt1180-DX1PB-ALLBULK-DNCD1aPos_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/BAMS/pt1180-DX1PB-ALLBULK-iSPCD4_dedup.bam \
  -b /hpc/pmc_vanboxtel/processed/231017_HMFreg2090_pt1180_Bulk/BAMS/pt1180-DX1PB-ALLBULK-SPCD4_dedup.bam \
  -n pt1180-DX1BM-MSCBULK \
  -c /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/config_15cov.ini


 