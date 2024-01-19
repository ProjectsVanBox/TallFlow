#!/bin/bash
#SBATCH -t 1:0:0
#SBATCH --mem 10G
#SBATCH -c 2

TARGET=/hpc/pmc_vanboxtel/_old_structure/data/ENRICH/MIPs/81_snps_mip_design_nijmegen_sort_hg38_liftover.vcf
GATKDIR=/hpc/local/CentOS7/pmc_vanboxtel/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/
GENOME=/hpc/pmc_vanboxtel/_old_structure/data/homo_sapiens.GRCh38.GATK.illumina/Homo_sapiens_assembly38.fasta
module load Java/1.8.0_60

for BAM in /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/FingerPrint/*.bam; do 
	SAMPLE=$(basename $BAM)
	SAMPLE=${SAMPLE/_dedup.bam/}
	OUT=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FingerPrint/pt3292_pt344/${SAMPLE}.vcf
	TMP=/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FingerPrint/pt3292_pt344/tmp/${SAMPLE}/
	mkdir -p $TMP
	java -Djava.io.tmpdir=$TMP -Xmx10G -jar $GATKDIR/GenomeAnalysisTK.jar \
	-T UnifiedGenotyper \
	-R $GENOME \
	-L $TARGET \
	-D $TARGET \
	-I $BAM \
	-o $OUT \
	--output_mode EMIT_ALL_SITES 
done

exit 0
