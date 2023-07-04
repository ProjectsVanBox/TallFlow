#!/bin/bash

GUIX_PATH=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/.guix-profile
gripss_somatic_vcf=$1
gripss_somatic_filtered_vcf=${gripss_somatic_vcf/.vcf.gz/.filtered.vcf.gz}

gripss_jar=/hpc/local/CentOS7/pmc_vanboxtel/bin/gridss-purple-linx_v1.3.2/hmftools/gripss-1.9.jar
	
guixr load-profile ${GUIX_PATH} -- << EOF
java -Xmx80G -cp ${gripss_jar} com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
-input_vcf ${gripss_somatic_vcf} \
-output_vcf ${gripss_somatic_filtered_vcf} \
EOF
