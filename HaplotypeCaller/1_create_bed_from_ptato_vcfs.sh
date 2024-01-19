#!/bin/bash

#### pt3291
patientID="pt3291"

cd "/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/HaplotypeCaller/pt3291/"

for VCF in /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/HaplotypeCaller/pt3291/*.vcf; do
	echo ${VCF}
	grep -v "^#" ${VCF} | awk '{print $1"\t"$2-1"\t"$2}' >> /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/HaplotypeCaller/pt3291/${patientID}_interval_list.bed
done

grep "^@" /hpc/pmc_vanboxtel/resources/homo_sapiens.GRCh38.GATK.illumina/genome.interval_list > /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/HaplotypeCaller/pt3291/${patientID}_interval_list_merged.interval_list

bedtools sort -i /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/HaplotypeCaller/pt3291/${patientID}_interval_list.bed | bedtools merge | awk '{print $0"\t+\t."}' >> /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/HaplotypeCaller/pt3291/${patientID}_interval_list_merged.interval_list

