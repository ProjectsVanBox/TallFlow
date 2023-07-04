#!/bin/bash
#SBATCH --job-name=SMuRF
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --time=24:0:0
#SBATCH -e SMuRF.err
#SBATCH -o SMuRF.log

cd /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt344

. /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/venv_3.6/bin/activate
/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip -c /hpc/pmc_vanboxtel/processed/230612_HMFreg2010_pt344_Bulk/VCFS/VCF/230612_HMFreg2010_pt344_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf > /hpc/pmc_vanboxtel/processed/230612_HMFreg2010_pt344_Bulk/VCFS/VCF/230612_HMFreg2010_pt344_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz
/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/tabix /hpc/pmc_vanboxtel/processed/230612_HMFreg2010_pt344_Bulk/VCFS/VCF/230612_HMFreg2010_pt344_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz

python /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/SMuRF.py \
  -i /hpc/pmc_vanboxtel/processed/230612_HMFreg2010_pt344_Bulk/VCFS/VCF/230612_HMFreg2010_pt344_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.vcf.gz \
  -b /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/BULK_SEQ/BAMS/pt344/pt344-DX1BM-ALLBULK_dedup.bam \
  -b /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/BULK_SEQ/BAMS/pt344/pt344-DX1PB-ALLBULK-DNCD1aPos_dedup.bam \
  -b /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/BULK_SEQ/BAMS/pt344/pt344-DX1PB-ALLBULK-DP_dedup.bam \
  -b /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/BULK_SEQ/BAMS/pt344/pt344-DX1PB-ALLBULK-SPCD8_dedup.bam \
  -b /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/BULK_SEQ/BAMS/pt344/pt344-DX1PB-MONOBULK_dedup.bam \
  -n pt344-DX1PB-MONOBULK \
  -c /hpc/pmc_vanboxtel/projects/TallFlow/2_Code/BulkSeq_SMuRF/SMuRF_config15x.ini

#bash /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/scripts/split_mutation_type.sh /hpc/pmc_vanboxtel/projects/Ageline/2_Code/2_scripts/10817/10817.filtered_variants_snpEff_snpSift_Cosmicv89_SMuRF_filtered.vcf

#bash /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/scripts/split_in_single_sample_vcfs.sh /hpc/pmc_vanboxtel/projects/Ageline/2_Code/2_scripts/10817/10817.filtered_variants_snpEff_snpSift_Cosmicv89_SMuRF_filtered.vcf

#/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip -c /hpc/pmc_vanboxtel/projects/Ageline/2_Code/2_scripts/10817/10817.filtered_variants_snpEff_snpSift_Cosmicv89_SMuRF.vcf > /hpc/pmc_vanboxtel/projects/Ageline/2_Code/2_scripts/10817/10817.filtered_variants_snpEff_snpSift_Cosmicv89_SMuRF.vcf.gz
#/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/tabix /hpc/pmc_vanboxtel/projects/Ageline/2_Code/2_scripts/10817/10817.filtered_variants_snpEff_snpSift_Cosmicv89_SMuRF.vcf.gz
#python /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/scripts/driver_mutations_filter.py -i /hpc/pmc_vanboxtel/projects/Ageline/2_Code/2_scripts/10817/10817.filtered_variants_snpEff_snpSift_Cosmicv89_SMuRF.vcf.gz -c /hpc/pmc_vanboxtel/projects/Ageline/2_Code/1_ini_files/SMuRF_config.ini
