#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(MutationalPatterns)
library(BSgenome)
library(ChIPpeakAnno)
library(ggplot2)
library(NMF)
library(RColorBrewer)
library(tibble)
library(reshape2)
library(grid)
mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2 <- brewer.pal(8, "Dark2")


### Read in genomes
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)


### Get signatures 
signatures = get_known_signatures()
pta_sig = read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/signatures/PTA_Artefact_Signature.txt", sep = "\t", header = T)
pta_sig = as.matrix(pta_sig)
PTA <- as.numeric(pta_sig[,"PTA"])
hspc_sig = read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/signatures/sigProfiler_SBS_working_signatures_incl_hspc.txt", sep = "\t", header = T)
hspc_sig = as.matrix(hspc_sig)
HSPC <- as.numeric(hspc_sig[,"HSPC"])
signatures <- cbind(PTA, HSPC, signatures)

colnames(get_known_signatures(source = "COSMIC_v3.2"))
### Set and create directories 
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/")


### Gather the files 
pt2322_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
#pt2322_bulk_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/BulkSeqs/",
#                                    pattern = "*filtered.vcf$", full.names = TRUE, recursive = TRUE )
pt2229_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2229/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
#pt2229_bulk_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2229/BulkSeqs/",
#                                    pattern = "*filtered.vcf$", full.names = TRUE, recursive = TRUE )
#combined.files <- c(pt2322_vcf_files, pt2322_bulk_vcf_files[1:5], pt2229_vcf_files, pt2229_bulk_vcf_files)
pt2283D_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/snvs/pt2283D/",
                                pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt344_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt344/snvs/pt344/",
                                pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt3291_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt3291/snvs/pt3291/",
                              pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt1180_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt1180/snvs/pt1180/",
                              pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt1909_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt1909/snvs/pt1909/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )

### Sample names 
pt2322_sample_names <- c("pt2322_ALLBULK", "pt2322_TALL1", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12", "pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4",
                         "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", "pt2322_TALL9") 
pt2229_sample_names <- c("pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL10", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", 
                  "pt2229_TALL7", "pt2229_TALL8", "pt2229_TALL9")
pt2283D_sample_names <- c("pt2283D-DX1BM-ALLBULK", "pt2283D-DX1BM-TALL1", "pt2283D-DX1BM-TALL10", "pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12", 
                          "pt2283D-DX1BM-TALL2", "pt2283D-DX1BM-TALL3", "pt2283D-DX1BM-TALL4", "pt2283D-DX1BM-TALL5", "pt2283D-DX1BM-TALL6", 
                          "pt2283D-DX1BM-TALL7", "pt2283D-DX1BM-TALL8", "pt2283D-DX1BM-TALL9")
pt344_sample_names <- c("pt344-DX1BM-ALLBULK", "pt344-DX1PB-TALL1", "pt344-DX1PB-TALL2", "pt344-DX1PB-TALL3", "pt344-DX1PB-TALL4", "pt344-DX1PB-TALL5", 
                        "pt344-DX1PB-TALL6", "pt344-DX1PB-TALL7", "pt344-DX1PB-TALL8", "pt344-DX1PB-TALL9")
pt3291_sample_names <- c("pt3291-DX1BM-TALL10", "pt3291-DX1BM-TALL11", "pt3291-DX1BM-TALL12", "pt3291-DX1BM-TALL5", "pt3291-DX1BM-TALL6", 
                         "pt3291-DX1BM-TALL7", "pt3291-DX1BM-TALL8", "pt3291-DX1BM-TALL9", "pt3291-DX1PB-TALL1", "pt3291-DX1PB-TALL2", "pt3291-DX1PB-TALL3", 
                         "pt3291-DX1PB-TALL4")
pt1180_sample_names <- c("pt1180-DX1PB-TALL1", "pt1180-DX1PB-TALL10", "pt1180-DX1PB-TALL11", "pt1180-DX1PB-TALL12", "pt1180-DX1PB-TALL2", "pt1180-DX1PB-TALL3", 
                         "pt1180-DX1PB-TALL4", "pt1180-DX1PB-TALL5", "pt1180-DX1PB-TALL6", "pt1180-DX1PB-TALL7", "pt1180-DX1PB-TALL8", "pt1180-DX1PB-TALL9")
pt1909_sample_names <- c("pt1909-DX1BM-TALL1", "pt1909-DX1BM-TALL2", "pt1909-DX1BM-TALL3", "pt1909-DX1BM-TALL4", "pt1909-DX1BM-TALL5", "pt1909-DX1BM-TALL6", 
                         "pt1909-DX1BM-TALL7", "pt1909-DX1BM-TALL8", "pt1909-DX1BM-TALL9")


### Branched VCF files 
pt2322_vcf_branches <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/Branches_MutPat/pt2322/",
                                  pattern = "*.vcf", full.names = TRUE, recursive = TRUE )
pt2229_vcf_branches <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/Branches_MutPat/pt2229/",
                                  pattern = "*.vcf", full.names = TRUE, recursive = TRUE )
pt2283D_vcf_branches <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283D_new/manual_filter/branches_vcfs/",
                                pattern = "*.vcf$", full.names = TRUE, recursive = FALSE )
pt344_vcf_branches <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt344/branches_vcfs/",
                              pattern = "*.vcf$", full.names = TRUE, recursive = FALSE )
pt3291_vcf_branches <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt3291/branches_vcfs/",
                               pattern = "*.vcf$", full.names = TRUE, recursive = FALSE )
pt1180_vcf_branches <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1180/branches_vcfs/",
                               pattern = "*.vcf$", full.names = TRUE, recursive = FALSE )
pt1909_vcf_branches <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt1909/branches_vcfs/",
                               pattern = "*.vcf$", full.names = TRUE, recursive = FALSE )

# "pt2322_branches", 
pt2322_branch_names <- c("pt2322_root_branch", "pt2322_TALL1_branch", "pt2322_TALL1_TALL3_branch", "pt2322_TALL1_TALL3_TALL2_TALL7_TALL12_branch", 
                         "pt2322_TALL1_TALL3_TALL2_TALL7_TALL12_TALL10_TALL9_TALL11_branch", "pt2322_TALL10_branch", "pt2322_TALL10_TALL9_TALL12_branch", 
                         "pt2322_TALL11_branch", "pt2322_TALL12_branch", "pt2322_TALL2_branch", "pt2322_TALL2_TALL7_TALL12_branch", "pt2322_TALL3_branch", 
                         "pt2322_TALL4_branch", "pt2322_TALL4_TALL5_branch", "pt2322_TALL5_branch", "pt2322_TALL6_branch", "pt2322_TALL6_TALL4_TALL5_branch", 
                         "pt2322_TALL6_TALL4_TALL5_TALL8_branch", "pt2322_TALL7_branch", "pt2322_TALL7_TALL12_branch", "pt2322_TALL8_branch", 
                         "pt2322_TALL9_branch", "pt2322_TALL9_TALL11_branch")
pt2229_branch_names <- c("pt2229_root_branch", "pt2229_TALL1_TALL2_branch", "pt2229_TALL10_TALL8_TALL9_branch", "pt2229_TALL10_TALL8_TALL9_TALL5_TALL3_TALL4_branch", 
                         "pt2229_TALL10_TALL8_TALL9_TALL5_TALL3_TALL4_TALL7_branch", "pt2229_TALL3_TALL4_branch", "pt2229_TALL5_TALL3_TALL4_branch", 
                         "pt2229_TALL6_TALL1_TALL2_branch", "pt2229_TALL8_TALL9_branch")
#sample_names <- c("pt2322_ALLBULK", "pt2322_TALL1", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12", "pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4",
#                  "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", "pt2322_TALL9", "pt2322_ALLBULK-DN", "pt2322_ALLBULK-DP",
#                  "pt2322_ALLBULK-iSPCD4", "pt2322_ALLBULK-SPCD4", "pt2322_ALLBULK_2", 
#                  "pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL10", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", 
#                  "pt2229_TALL7", "pt2229_TALL8", "pt2229_TALL9", "pt2229_BULK_1", "pt2229_BULK_2")
pt2283D_branch_names <- c("pt2283D_root_branch.vcf", "pt2283D_TALL1_branch.vcf", 
                         "pt2283D_TALL1_TALL2_TALL3_TALL4_TALL5_TALL6_TALL7_TALL8_TALL9_TALL10_TALL11_branch.vcf", 
                         "pt2283D_TALL1_TALL2_TALL4_TALL5_TALL6_TALL7_TALL8_TALL9_TALL10_TALL11_branch.vcf", "pt2283D_TALL10_branch.vcf", 
                         "pt2283D_TALL11_branch.vcf", "pt2283D_TALL12_branch.vcf", "pt2283D_TALL2_branch.vcf", "pt2283D_TALL3_branch.vcf", 
                         "pt2283D_TALL4_branch.vcf", "pt2283D_TALL5_branch.vcf", "pt2283D_TALL6_branch.vcf", "pt2283D_TALL7_branch.vcf", 
                         "pt2283D_TALL8_branch.vcf", "pt2283D_TALL9_branch.vcf")
pt344_branch_names <- c("pt344_root.vcf", "pt344_TALL1_branch.vcf", "pt344_TALL1_TALL5_branch.vcf", "pt344_TALL1_TALL5_TALL2_branch.vcf", 
                        "pt344_TALL1_TALL5_TALL2_TALL8_TALL7_TALL9_TALL3_branch.vcf", "pt344_TALL1_TALL5_TALL2_TALL8_TALL7_TALL9_TALL3_TALL6_branch.vcf", 
                        "pt344_TALL2_branch.vcf", "pt344_TALL3_branch.vcf", "pt344_TALL4_branch.vcf", "pt344_TALL5_branch.vcf", "pt344_TALL6_branch.vcf", 
                        "pt344_TALL7_branch.vcf", "pt344_TALL8_branch.vcf", "pt344_TALL8_TALL7_branch.vcf", "pt344_TALL8_TALL7_TALL9_branch.vcf", 
                        "pt344_TALL9_branch.vcf")
pt3291_branch_names <- c("pt3291_root_branch.vcf", "pt3291_TALL1_branch.vcf", "pt3291_TALL10_branch.vcf", "pt3291_TALL11_branch.vcf", 
                         "pt3291_TALL12_branch.vcf", "pt3291_TALL12_TALL10_branch.vcf", "pt3291_TALL2_branch.vcf", "pt3291_TALL3_branch.vcf", 
                         "pt3291_TALL4_branch.vcf", "pt3291_TALL5_branch.vcf", "pt3291_TALL6_branch.vcf", "pt3291_TALL7_branch.vcf", 
                         "pt3291_TALL8_branch.vcf", "pt3291_TALL9_branch.vcf", "pt3291_TALL9_TALL6_branch.vcf", "pt3291_TALL9_TALL6_TALL8_branch.vcf", 
                         "pt3291_TALL9_TALL6_TALL8_TALL11_branch.vcf", "pt3291_TALL9_TALL6_TALL8_TALL11_TALL3_branch.vcf", 
                         "pt3291_TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_branch.vcf", 
                         "pt3291_TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_branch.vcf", 
                         "pt3291_TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_branch.vcf", 
                         "pt3291_TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_TALL2_branch.vcf", 
                         "pt3291_TALL9_TALL6_TALL8_TALL11_TALL3_TALL12_TALL10_TALL5_TALL7_TALL2_TALL1_branch.vcf")
pt1180_branch_names <- c("pt1180_header.vcf", "pt1180_root_branch.vcf", "pt1180_TALL1_branch.vcf", "pt1180_TALL10_branch.vcf", 
                         "pt1180_TALL10_TALL1_branch.vcf", "pt1180_TALL10_TALL1_TALL2_branch.vcf", "pt1180_TALL11_branch.vcf", 
                         "pt1180_TALL11_TALL5_branch.vcf", "pt1180_TALL11_TALL5_TALL3_branch.vcf", "pt1180_TALL11_TALL5_TALL3_TALL6_branch.vcf", 
                         "pt1180_TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_branch.vcf", 
                         "pt1180_TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_TALL10_TALL1_TALL2_branch.vcf", 
                         "pt1180_TALL11_TALL5_TALL3_TALL6_TALL4_TALL12_TALL10_TALL1_TALL2_TALL9_TALL7_branch.vcf", "pt1180_TALL12_branch.vcf", 
                         "pt1180_TALL2_branch.vcf", "pt1180_TALL3_branch.vcf", "pt1180_TALL4_branch.vcf", "pt1180_TALL4_TALL12_branch.vcf", 
                         "pt1180_TALL5_branch.vcf", "pt1180_TALL6_branch.vcf", "pt1180_TALL7_branch.vcf", "pt1180_TALL8_branch.vcf", 
                         "pt1180_TALL9_branch.vcf", "pt1180_TALL9_TALL7_branch.vcf")
pt1909_branch_names <- c("pt1909_root_branch.vcf", "pt1909_TALL1_branch.vcf", "pt1909_TALL1_TALL5_branch.vcf", "pt1909_TALL1_TALL5_TALL4_branch.vcf", 
                         "pt1909_TALL1_TALL5_TALL4_TALL9_branch.vcf", "pt1909_TALL1_TALL5_TALL4_TALL9_TALL7_branch.vcf", 
                         "pt1909_TALL1_TALL5_TALL4_TALL9_TALL7_TALL6_TALL2_TALL3_branch.vcf", "pt1909_TALL2_branch.vcf", "pt1909_TALL3_branch.vcf", 
                         "pt1909_TALL4_branch.vcf", "pt1909_TALL5_branch.vcf", "pt1909_TALL6_branch.vcf", "pt1909_TALL6_TALL2_branch.vcf", 
                         "pt1909_TALL6_TALL2_TALL3_branch.vcf", "pt1909_TALL7_branch.vcf", "pt1909_TALL8_branch.vcf", "pt1909_TALL9_branch.vcf")




pt2322_branch <- c("root_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", 
                   "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "ISP_branch", "ISP_branch", 
                   "ISP_branch", "ISP_branch", "ISP_branch", "ISP_DP_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "ISP_DP_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch")
pt2229_branch <- c("pt2229_root_branch", "pt2229_DN_ISP_branch", "pt2229_Other_branch", "pt2229_Other_branch", "pt2229_Other_branch", "pt2229_Other_branch", 
                   "pt2229_Other_branch", "pt2229_DN_ISP_branch", "pt2229_Other_branch")


### Bulkseq vcfs 
pt2322_vcf_bulks <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2322/",
                                  pattern = ".*_pt2322-DX1BM.*.vcf", full.names = TRUE, recursive = TRUE )
pt2229_vcf_bulks <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2229/",
                               pattern = "*.vcf", full.names = TRUE, recursive = TRUE )
pt2283D_vcf_bulks <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2283D/",
                               pattern = ".*pt2283D-DX1BM.*.vcf", full.names = TRUE, recursive = TRUE )
pt344_vcf_bulks <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt344/",
                               pattern = ".*pt344-DX1.*.vcf", full.names = TRUE, recursive = TRUE )
pt3291_vcf_bulks <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt3291/",
                               pattern = ".*pt3291-DX1PB.*.vcf", full.names = TRUE, recursive = TRUE )
pt1180_vcf_bulks <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt1180/",
                               pattern = ".*pt1180-DX1.*.vcf", full.names = TRUE, recursive = TRUE )
pt1909_vcf_bulks <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt1909/",
                               pattern = ".*pt1909-DX1.*.vcf", full.names = TRUE, recursive = TRUE )
# Create single sample vcfs 
pt2322_bulks_names <- c("pt2322-DX1BM-ALLBULK-DN", "pt2322-DX1BM-ALLBULK-DP", "pt2322-DX1BM-ALLBULK-iSPCD4", "pt2322-DX1BM-ALLBULK-SPCD4", "pt2322-DX1BM-ALLBULK", "pt2322-DX1BM-MONOBULK")
#pt2229_bulks_names <- c()
pt2283D_bulks_names <- c("pt2283D-DX1BM-ALLBULK-DN", "pt2283D-DX1BM-ALLBULK-DP", "pt2283D-DX1BM-ALLBULK-iSPCD4", "pt2283D-DX1BM-ALLBULK-SPCD4", "pt2283D-DX1BM-ALLBULK", "pt2283D-DX1BM-MSCBULK")
pt344_bulks_names <- c("pt344-DX1BM-ALLBULK", "pt344-DX1PB-ALLBULK-DNCD1aPos", "pt344-DX1PB-ALLBULK-DP", "pt344-DX1PB-ALLBULK-SPCD8", "pt344-DX1PB-MONOBULK")
pt3291_bulks_names <- c("pt3291-DX1PB-ALLBULK-DNCD1a-min", "pt3291-DX1PB-ALLBULK-DNCD1a-plus", "pt3291-DX1PB-ALLBULK-DP", "pt3291-DX1PB-ALLBULK-SPCD8", "pt3291-DX1PB-ALLBULK")
pt1180_bulks_names <- c("pt1180-DX1BM-MSCBULK", "pt1180-DX1PB-ALLBULK-DNCD1aNeg", "pt1180-DX1PB-ALLBULK-DNCD1aPos", "pt1180-DX1PB-ALLBULK-iSPCD4", "pt1180-DX1PB-ALLBULK-SPCD4", "pt1180-DX1PB-ALLBULK")
pt1909_bulks_names <- c("pt1909-DX1BM-ALLBULK-DNCD1aNeg", "pt1909-DX1BM-ALLBULK-DNCD1aPos", "pt1909-DX1BM-ALLBULK")


pt2322_bulks_grl <- read_vcfs_as_granges(pt2322_vcf_bulks, pt2322_bulks_names, ref_genome, type = 'all')
pt2322_bulks_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2322_bulks_grl, 'snv'), ref_genome = ref_genome)
#pt2229_bulks_grl <- read_vcfs_as_granges(pt2229_vcf_bulks, pt2229_bulks_names, ref_genome, type = 'all')
#pt2229_bulks_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2229_bulks_grl, 'snv'), ref_genome = ref_genome)
pt2283D_bulks_grl <- read_vcfs_as_granges(pt2283D_vcf_bulks, pt2283D_bulks_names, ref_genome, type = 'all')
pt2283D_bulks_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2283D_bulks_grl, 'snv'), ref_genome = ref_genome)
pt344_bulks_grl <- read_vcfs_as_granges(pt344_vcf_bulks, pt344_bulks_names, ref_genome, type = 'all')
pt344_bulks_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt344_bulks_grl, 'snv'), ref_genome = ref_genome)
pt3291_bulks_grl <- read_vcfs_as_granges(pt3291_vcf_bulks, pt3291_bulks_names, ref_genome, type = 'all')
pt3291_bulks_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt3291_bulks_grl, 'snv'), ref_genome = ref_genome)
pt1180_bulks_grl <- read_vcfs_as_granges(pt1180_vcf_bulks, pt1180_bulks_names, ref_genome, type = 'all')
pt1180_bulks_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1180_bulks_grl, 'snv'), ref_genome = ref_genome)
pt1909_bulks_grl <- read_vcfs_as_granges(pt1909_vcf_bulks, pt1909_bulks_names, ref_genome, type = 'all')
pt1909_bulks_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1909_bulks_grl, 'snv'), ref_genome = ref_genome)

### Create matrices 
#pt2322
pt2322_grl <- read_vcfs_as_granges(pt2322_vcf_files, pt2322_sample_names, ref_genome, type = 'all')
pt2322_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2322_grl, 'snv'), ref_genome = ref_genome)
#pt2229
pt2229_grl <- read_vcfs_as_granges(pt2229_vcf_files, pt2229_sample_names, ref_genome, type = 'all')
pt2229_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2229_grl, 'snv'), ref_genome = ref_genome)
#pt2283D
pt2283D_grl <- read_vcfs_as_granges(pt2283D_vcf_files, pt2283D_sample_names, ref_genome, type = 'all')
pt2283D_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2283D_grl, 'snv'), ref_genome = ref_genome)
#pt344
pt344_grl <- read_vcfs_as_granges(pt344_vcf_files, pt344_sample_names, ref_genome, type = 'all')
pt344_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt344_grl, 'snv'), ref_genome = ref_genome)
#pt3291
pt3291_grl <- read_vcfs_as_granges(pt3291_vcf_files, pt3291_sample_names, ref_genome, type = 'all')
pt3291_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt3291_grl, 'snv'), ref_genome = ref_genome)
#pt1180
pt1180_grl <- read_vcfs_as_granges(pt1180_vcf_files, pt1180_sample_names, ref_genome, type = 'all')
pt1180_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1180_grl, 'snv'), ref_genome = ref_genome)
#pt1909
pt1909_grl <- read_vcfs_as_granges(pt1909_vcf_files, pt1909_sample_names, ref_genome, type = 'all')
pt1909_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1909_grl, 'snv'), ref_genome = ref_genome)


### Create mutational matrix  
nmf_ref_mat = read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/signatures/NMF_reference_matrix.txt", sep = "\t", header = T)
mut_mat_own <- cbind(pt2322_mut_mat, pt2229_mut_mat, pt2283D_mut_mat, pt344_mut_mat, pt3291_mut_mat, pt1180_mut_mat, pt1909_mut_mat)
mut_mat <- cbind(mut_mat_own, nmf_ref_mat)
### De novo mutational signature extraction using NMF
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank = 2:10, method = "brunet", 
                nrun = 50, seed = 123456, .opt = "v-p")
pdf("./3_Output/FullDenovoMutPatterns/estimates.pdf")
plot(estimate)
dev.off()


### Extract signatures using the additional NMF samples 
# 3 ranks
nmf_res3 <- extract_signatures(mut_mat, rank = 3, nrun = 50, single_core = TRUE)
colnames(nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
# 4 ranks
nmf_res4 <- extract_signatures(mut_mat, rank = 4, nrun = 50, single_core = TRUE)
colnames(nmf_res4$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res4$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
# 5 ranks
nmf_res5 <- extract_signatures(mut_mat, rank = 5, nrun = 50, single_core = TRUE)
colnames(nmf_res5$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(nmf_res5$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
# 6 ranks
nmf_res6 <- extract_signatures(mut_mat, rank = 6, nrun = 50, single_core = TRUE)
colnames(nmf_res6$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(nmf_res6$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
# 7 ranks
nmf_res7 <- extract_signatures(mut_mat, rank = 7, nrun = 50, single_core = TRUE)
colnames(nmf_res7$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
rownames(nmf_res7$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
# 8 ranks
nmf_res8 <- extract_signatures(mut_mat, rank = 8, nrun = 50, single_core = TRUE)
colnames(nmf_res8$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H")
rownames(nmf_res8$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H")
# 9 ranks
nmf_res9 <- extract_signatures(mut_mat, rank = 9, nrun = 50, single_core = TRUE)
colnames(nmf_res9$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H", "Signature I")
rownames(nmf_res9$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H", "Signature I")
# 10 ranks
nmf_res10 <- extract_signatures(mut_mat, rank = 10, nrun = 50, single_core = TRUE)
colnames(nmf_res10$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                    "Signature H", "Signature I", "Signature J")
rownames(nmf_res10$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                      "Signature H", "Signature I", "Signature J")

### Check with known signatures 
# 3 ranks
nmf_res3 <- rename_nmf_signatures(nmf_res3, signatures, cutoff = 0.85)
colnames(nmf_res3$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank3_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank3_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res3$contribution, nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank3_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res3$contribution[,1:length(colnames(mut_mat_own))], nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 4 ranks
nmf_res4 <- rename_nmf_signatures(nmf_res4, signatures, cutoff = 0.85)
colnames(nmf_res4$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank4_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res4$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank4_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res4$contribution, nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank4_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res4$contribution[,1:length(colnames(mut_mat_own))], nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 5 ranks
nmf_res5 <- rename_nmf_signatures(nmf_res5, signatures, cutoff = 0.85)
colnames(nmf_res5$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank5_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res5$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank5_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res5$contribution, nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank5_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res5$contribution[,1:length(colnames(mut_mat_own))], nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 6 ranks
nmf_res6 <- rename_nmf_signatures(nmf_res6, signatures, cutoff = 0.85)
colnames(nmf_res6$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank6_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res6$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank6_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res6$contribution, nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank6_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res6$contribution[,1:length(colnames(mut_mat_own))], nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 7 ranks
nmf_res7 <- rename_nmf_signatures(nmf_res7, signatures, cutoff = 0.85)
colnames(nmf_res7$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank7_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res7$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank7_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res7$contribution, nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank7_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res7$contribution[,1:length(colnames(mut_mat_own))], nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 8 ranks
nmf_res8 <- rename_nmf_signatures(nmf_res8, signatures, cutoff = 0.85)
colnames(nmf_res8$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank8_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res8$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank8_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res8$contribution, nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank8_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res8$contribution[,1:length(colnames(mut_mat_own))], nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 9 ranks
nmf_res9 <- rename_nmf_signatures(nmf_res9, signatures, cutoff = 0.91)
colnames(nmf_res9$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank9_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res9$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank9_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res9$contribution, nmf_res9$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank9_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res9$contribution[,1:length(colnames(mut_mat_own))], nmf_res9$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 10 ranks
nmf_res10 <- rename_nmf_signatures(nmf_res10, signatures, cutoff = 0.85)
colnames(nmf_res10$signatures)
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank10_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res10$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank10_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res10$contribution, nmf_res10$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/RankSignatures/rank10_signature_contribution_short.pdf", width = 12)
plot_contribution(nmf_res10$contribution[,1:length(colnames(mut_mat_own))], nmf_res10$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




plot_compare_profiles(mut_mat[, 1],
                      nmf_res3$reconstructed[, 1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)
plot_compare_profiles(mut_mat[, 1],
                      nmf_res5$reconstructed[, 1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)
plot_compare_profiles(mut_mat[, 1],
                      nmf_res7$reconstructed[, 1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)
plot_compare_profiles(mut_mat[, 1],
                      nmf_res8$reconstructed[, 1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

### Match with COSMIC signatures
SBS_cosim <- cos_sim_matrix(mut_mat_own, signatures)
plot_cosine_heatmap(SBS_cosim, 
                    cluster_rows = TRUE, cluster_cols = TRUE)

cos_list <- c()
for (c in colnames(SBS_cosim)){
  if (any(SBS_cosim[,c] > 0.75 )){
    print(c)
    cos_list <- append(cos_list, c)
  }
}
pdf("./3_Output/FullDenovoMutPatterns/SBS_Cosine_0.75.pdf", width = 14, height = 9)
plot_cosine_heatmap(SBS_cosim[,colnames(SBS_cosim) %in% cos_list ], 
                    cluster_rows = FALSE, cluster_cols = TRUE, plot_values = TRUE)
dev.off()

plot_96_profile(mut_mat_own[,c("pt3291-DX1BM-TALL5", "pt3291-DX1BM-TALL10", "pt1909-DX1BM-TALL3", "pt1909-DX1BM-TALL7")])

# Bootstrapped refit 
contri_boots <- fit_to_signatures_bootstrapped(mut_mat_own,
                                                  signatures,
                                                  n_boots = 50,
                                                  method = "strict"
)


plot_bootstrapped_contribution(contri_boots)
for (sample_name in colnames(mut_mat_own)){
  print(sample_name)
  bs_plot <- plot_bootstrapped_contribution(contri_boots[rownames(contri_boots) %in% c(paste0(sample_name, "_", 1:50)),])
  pdf(paste0("./3_Output/FullDenovoMutPatterns/Bootstrapped/", sample_name, "_bootstrapped.pdf"), width = 14, height = 9)
  print(bs_plot)
  dev.off()
}



### De novo mutational signature extraction using NMF
mut_mat_own <- mut_mat_own + 0.0001
estimate_own <- nmf(mut_mat_own, rank = 2:10, method = "brunet", 
                nrun = 50, seed = 123456, .opt = "v-p")
pdf("./3_Output/FullDenovoMutPatterns/estimates_own.pdf")
plot(estimate_own)
dev.off()


### Extract signatures using our own samples 
# 3 ranks
own_nmf_res3 <- extract_signatures(mut_mat_own, rank = 3, nrun = 50, single_core = TRUE)
colnames(own_nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(own_nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
# 4 ranks
own_nmf_res4 <- extract_signatures(mut_mat_own, rank = 4, nrun = 50, single_core = TRUE)
colnames(own_nmf_res4$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(own_nmf_res4$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
# 5 ranks
own_nmf_res5 <- extract_signatures(mut_mat_own, rank = 5, nrun = 50, single_core = TRUE)
colnames(own_nmf_res5$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(own_nmf_res5$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
# 6 ranks
own_nmf_res6 <- extract_signatures(mut_mat_own, rank = 6, nrun = 50, single_core = TRUE)
colnames(own_nmf_res6$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(own_nmf_res6$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
# 7 ranks
own_nmf_res7 <- extract_signatures(mut_mat_own, rank = 7, nrun = 50, single_core = TRUE)
colnames(own_nmf_res7$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
rownames(own_nmf_res7$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
# 8 ranks
own_nmf_res8 <- extract_signatures(mut_mat_own, rank = 8, nrun = 50, single_core = TRUE)
colnames(own_nmf_res8$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                       "Signature H")
rownames(own_nmf_res8$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                         "Signature H")
# 9 ranks
own_nmf_res9 <- extract_signatures(mut_mat_own, rank = 9, nrun = 50, single_core = TRUE)
colnames(own_nmf_res9$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                       "Signature H", "Signature I")
rownames(own_nmf_res9$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                         "Signature H", "Signature I")
# 10 ranks
own_nmf_res10 <- extract_signatures(mut_mat_own, rank = 10, nrun = 50, single_core = TRUE)
colnames(own_nmf_res10$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                        "Signature H", "Signature I", "Signature J")
rownames(own_nmf_res10$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                          "Signature H", "Signature I", "Signature J")


### Check with known signatures 
# 3 ranks
own_nmf_res3 <- rename_nmf_signatures(own_nmf_res3, signatures, cutoff = 0.85)
colnames(own_nmf_res3$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank3_tri_nuc_profiles.pdf")
plot_96_profile(own_nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank3_signature_contribution.pdf", width = 12)
plot_contribution(own_nmf_res3$contribution, own_nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank3_signature_contribution_short.pdf", width = 12)
plot_contribution(own_nmf_res3$contribution[,1:length(colnames(mut_mat_own))], own_nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 4 ranks
own_nmf_res4 <- rename_nmf_signatures(own_nmf_res4, signatures, cutoff = 0.85)
colnames(own_nmf_res4$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank4_tri_nuc_profiles.pdf")
plot_96_profile(own_nmf_res4$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank4_signature_contribution.pdf", width = 12)
plot_contribution(own_nmf_res4$contribution, own_nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank4_signature_contribution_short.pdf", width = 12)
plot_contribution(own_nmf_res4$contribution[,1:length(colnames(mut_mat_own))], own_nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 5 ranks
own_nmf_res5 <- rename_nmf_signatures(own_nmf_res5, signatures, cutoff = 0.85)
colnames(own_nmf_res5$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank5_tri_nuc_profiles.pdf")
plot_96_profile(own_nmf_res5$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank5_signature_contribution.pdf", width = 12)
plot_contribution(own_nmf_res5$contribution, own_nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank5_signature_contribution_short.pdf", width = 12)
plot_contribution(own_nmf_res5$contribution[,1:length(colnames(mut_mat_own))], own_nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 6 ranks
own_nmf_res6 <- rename_nmf_signatures(own_nmf_res6, signatures, cutoff = 0.85)
colnames(own_nmf_res6$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank6_tri_nuc_profiles.pdf")
plot_96_profile(own_nmf_res6$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank6_signature_contribution.pdf", width = 12)
plot_contribution(own_nmf_res6$contribution, own_nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank6_signature_contribution_short.pdf", width = 12)
plot_contribution(own_nmf_res6$contribution[,1:length(colnames(mut_mat_own))], own_nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 7 ranks
own_nmf_res7 <- rename_nmf_signatures(own_nmf_res7, signatures, cutoff = 0.85)
colnames(own_nmf_res7$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank7_tri_nuc_profiles.pdf")
plot_96_profile(own_nmf_res7$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank7_signature_contribution.pdf", width = 12)
plot_contribution(own_nmf_res7$contribution, own_nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank7_signature_contribution_short.pdf", width = 12)
plot_contribution(own_nmf_res7$contribution[,1:length(colnames(mut_mat_own))], own_nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 8 ranks
own_nmf_res8 <- rename_nmf_signatures(own_nmf_res8, signatures, cutoff = 0.85)
colnames(own_nmf_res8$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank8_tri_nuc_profiles.pdf")
plot_96_profile(own_nmf_res8$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank8_signature_contribution.pdf", width = 12)
plot_contribution(own_nmf_res8$contribution, own_nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Own_RankSignatures/rank8_signature_contribution_short.pdf", width = 12)
plot_contribution(own_nmf_res8$contribution[,1:length(colnames(mut_mat_own))], own_nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



###################################################### Branch part ######################################################
#pt2322 --> why is there a merged_branch in here? 
pt2322_br_grl <- read_vcfs_as_granges(pt2322_vcf_branches[2:length(pt2322_vcf_branches)], pt2322_branch_names, ref_genome, type = 'all')
pt2322_br_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2322_br_grl, 'snv'), ref_genome = ref_genome)
#pt2229 --> Empty vcfs? 
pt2229_br_grl <- read_vcfs_as_granges(pt2229_vcf_branches, pt2229_branch_names, ref_genome, type = 'all')
pt2229_br_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2229_br_grl, 'snv'), ref_genome = ref_genome)
#pt2283D
pt2283D_br_grl <- read_vcfs_as_granges(pt2283D_vcf_branches, pt2283D_branch_names, ref_genome, type = 'all')
pt2283D_br_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt2283D_br_grl, 'snv'), ref_genome = ref_genome)
#pt344
pt344_br_grl <- read_vcfs_as_granges(pt344_vcf_branches, pt344_branch_names, ref_genome, type = 'all')
pt344_br_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt344_br_grl, 'snv'), ref_genome = ref_genome)
#pt3291 --> empty vcf? 
pt3291_br_grl <- read_vcfs_as_granges(pt3291_vcf_branches, pt3291_branch_names, ref_genome, type = 'all')
pt3291_br_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt3291_br_grl, 'snv'), ref_genome = ref_genome)
#pt1180 --> empty vcf?
pt1180_br_grl <- read_vcfs_as_granges(pt1180_vcf_branches, pt1180_branch_names, ref_genome, type = 'all')
pt1180_br_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1180_br_grl, 'snv'), ref_genome = ref_genome)
#pt1909 --> empty vcf?
pt1909_br_grl <- read_vcfs_as_granges(pt1909_vcf_branches, pt1909_branch_names, ref_genome, type = 'all')
pt1909_br_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1909_br_grl, 'snv'), ref_genome = ref_genome)
### Combine mutational matrices
mut_mat_br_own <- cbind(pt2322_br_mut_mat, pt2229_br_mut_mat, pt2283D_br_mut_mat, pt344_br_mut_mat, pt3291_br_mut_mat, pt1180_br_mut_mat, pt1909_br_mut_mat)
denovo_mut_mat_br <- mut_mat_br_own[,(colSums(mut_mat_br_own) > 100) > 0] 


### De novo mutational signature extraction using NMF
denovo_mut_mat_br <- denovo_mut_mat_br + 0.0001
br_estimate <- nmf(denovo_mut_mat_br, rank = 2:10, method = "brunet", 
                nrun = 50, seed = 123456, .opt = "v-p")
pdf("./3_Output/FullDenovoMutPatterns/br_estimates.pdf")
plot(br_estimate)
dev.off()


### Extract signatures for the branches of each samples  
# 3 ranks
br_nmf_res3 <- extract_signatures(denovo_mut_mat_br, rank = 3, nrun = 50, single_core = TRUE)
colnames(br_nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(br_nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
# 4 ranks
br_nmf_res4 <- extract_signatures(denovo_mut_mat_br, rank = 4, nrun = 50, single_core = TRUE)
colnames(br_nmf_res4$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(br_nmf_res4$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
# 5 ranks
br_nmf_res5 <- extract_signatures(denovo_mut_mat_br, rank = 5, nrun = 50, single_core = TRUE)
colnames(br_nmf_res5$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(br_nmf_res5$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
# 6 ranks
br_nmf_res6 <- extract_signatures(denovo_mut_mat_br, rank = 6, nrun = 50, single_core = TRUE)
colnames(br_nmf_res6$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(br_nmf_res6$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
# 7 ranks
br_nmf_res7 <- extract_signatures(denovo_mut_mat_br, rank = 7, nrun = 50, single_core = TRUE)
colnames(br_nmf_res7$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
rownames(br_nmf_res7$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
# 8 ranks
br_nmf_res8 <- extract_signatures(denovo_mut_mat_br, rank = 8, nrun = 50, single_core = TRUE)
colnames(br_nmf_res8$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                      "Signature H")
rownames(br_nmf_res8$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                        "Signature H")
# 9 ranks
br_nmf_res9 <- extract_signatures(denovo_mut_mat_br, rank = 9, nrun = 50, single_core = TRUE)
colnames(br_nmf_res9$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                      "Signature H", "Signature I")
rownames(br_nmf_res9$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                        "Signature H", "Signature I")
# 10 ranks
br_nmf_res10 <- extract_signatures(denovo_mut_mat_br, rank = 10, nrun = 50, single_core = TRUE)
colnames(br_nmf_res10$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                       "Signature H", "Signature I", "Signature J")
rownames(br_nmf_res10$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                         "Signature H", "Signature I", "Signature J")


### Check with known signatures 
# 3 ranks
br_nmf_res3 <- rename_nmf_signatures(br_nmf_res3, signatures, cutoff = 0.85)
colnames(br_nmf_res3$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank3_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank3_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res3$contribution, br_nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank3_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res3$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 4 ranks
br_nmf_res4 <- rename_nmf_signatures(br_nmf_res4, signatures, cutoff = 0.85)
colnames(br_nmf_res4$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank4_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res4$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank4_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res4$contribution, br_nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank4_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res4$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 5 ranks
br_nmf_res5 <- rename_nmf_signatures(br_nmf_res5, signatures, cutoff = 0.85)
colnames(br_nmf_res5$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank5_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res5$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank5_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res5$contribution, br_nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank5_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res5$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 6 ranks
br_nmf_res6 <- rename_nmf_signatures(br_nmf_res6, signatures, cutoff = 0.85)
colnames(br_nmf_res6$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank6_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res6$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank6_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res6$contribution, br_nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank6_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res6$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 7 ranks
br_nmf_res7 <- rename_nmf_signatures(br_nmf_res7, signatures, cutoff = 0.85)
colnames(br_nmf_res7$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank7_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res7$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank7_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res7$contribution, br_nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank7_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res7$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 8 ranks
br_nmf_res8 <- rename_nmf_signatures(br_nmf_res8, signatures, cutoff = 0.85)
colnames(br_nmf_res8$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank8_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res8$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank8_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res8$contribution, br_nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank8_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res8$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 9 ranks
br_nmf_res9 <- rename_nmf_signatures(br_nmf_res9, signatures, cutoff = 0.85)
colnames(br_nmf_res9$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank9_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res9$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank9_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res9$contribution, br_nmf_res9$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank9_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res9$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res9$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# 10 ranks
br_nmf_res10 <- rename_nmf_signatures(br_nmf_res10, signatures, cutoff = 0.85)
colnames(br_nmf_res10$signatures)
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank10_tri_nuc_profiles.pdf")
plot_96_profile(br_nmf_res10$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank10_signature_contribution.pdf", width = 12)
plot_contribution(br_nmf_res10$contribution, br_nmf_res10$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Br_RankSignatures/rank10_signature_contribution_short.pdf", width = 12)
plot_contribution(br_nmf_res10$contribution[,1:length(colnames(denovo_mut_mat_br))], br_nmf_res10$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


### Check strict refit with COSMIC signatures 
strict_refit <- fit_to_signatures_strict(denovo_mut_mat_br, signatures, max_delta = 0.05)
fig_list <- strict_refit$sim_decay_fig
for (i in 1:length(fig_list)){
  pdf(paste0("./3_Output/FullDenovoMutPatterns/br_StrictRefit/FigList/fig_list", i, ".pdf"))
  print(fig_list[[i]])
  dev.off()
}
fit_res_strict <- strict_refit$fit_res
pdf("./3_Output/FullDenovoMutPatterns/br_StrictRefit/br_StrictRefit_d0.05.pdf", height = 14, width = 14)
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


### Match with COSMIC signatures
SBS_cosim_br <- cos_sim_matrix(mut_mat_br_own, signatures)
cos_br_list <- c()
for (c in colnames(SBS_cosim_br)){
  if (any(SBS_cosim_br[,c] > 0.75, na.rm = TRUE)){
    print(c)
    cos_br_list <- append(cos_br_list, c)
  }
}
pdf("./3_Output/FullDenovoMutPatterns/SBS_branched_Cosine_0.75.pdf", width = 14, height = 9)
plot_cosine_heatmap(SBS_cosim_br[,colnames(SBS_cosim_br) %in% cos_br_list ], 
                    cluster_rows = FALSE, cluster_cols = FALSE, plot_values = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/ProfilePerRoot.pdf", width = 14, height = 9)
plot_96_profile(mut_mat_br_own[,c("pt1180_root_branch.vcf", "pt1909_root_branch.vcf", "pt344_root.vcf", "pt3291_root_branch.vcf", 
                                  "pt2283D_root_branch.vcf", "pt2229_root_branch", "pt2322_root_branch")])
dev.off()

# Bootstrapped refit 
contri_boots_br <- fit_to_signatures_bootstrapped(mut_mat_br_own,
                                               signatures,
                                               n_boots = 50,
                                               method = "strict"
)




plot_compare_profiles(mut_mat_br_own[, c("pt3291_root_branch.vcf")],
                      signatures[,"SBS87"],
                      profile_names = c("Root", "SBS87"),
                      condensed = TRUE
)

# Add bootstrapping, both for branches as single sample


# Remove low mut samples for the contri_boots function
NoLowMuts_br <- mut_mat_br_own[,(colSums(mut_mat_br_own) > 10) > 0] 
# Bootstrapped refit 
contri_boots_br <- fit_to_signatures_bootstrapped(NoLowMuts_br,
                                               signatures,
                                               n_boots = 50,
                                               method = "strict"
)


plot_bootstrapped_contribution(contri_boots_br)
for (sample_name in colnames(NoLowMuts_br)){
  print(sample_name)
  bs_plot <- plot_bootstrapped_contribution(contri_boots_br[rownames(contri_boots_br) %in% c(paste0(sample_name, "_", 1:50)),])
  pdf(paste0("./3_Output/FullDenovoMutPatterns/Bootstrapped_branched/", sample_name, "_bootstrapped.pdf"), width = 14, height = 9)
  print(bs_plot)
  dev.off()
}

br_nmf_res7$signatures = t(t(br_nmf_res7$signatures)/colSums(br_nmf_res7$signatures))
cust_strict_refit <- fit_to_signatures_strict(denovo_mut_mat_br, br_nmf_res7$signatures, max_delta = 0.1)
cust_fit_res_strict <- cust_strict_refit$fit_res
plot_contribution(cust_fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


####### Test refitting 
refit_cosine = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, br_nmf_res3$reconstructed))
data.frame(
  cosine = refit_cosine,
  nmuts = colSums(denovo_mut_mat_br)
) %>%
  ggplot(aes(x = nmuts, y  = cosine)) +
  geom_point() +
  lims(y = c(0.5, 1)) +
  geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
refit_cosine = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, br_nmf_res4$reconstructed))
data.frame(
  cosine = refit_cosine,
  nmuts = colSums(denovo_mut_mat_br)
) %>%
  ggplot(aes(x = nmuts, y  = cosine)) +
  geom_point() +
  lims(y = c(0.5, 1)) +
  geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
refit_cosine = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, br_nmf_res5$reconstructed))
data.frame(
  cosine = refit_cosine,
  nmuts = colSums(denovo_mut_mat_br)
) %>%
  ggplot(aes(x = nmuts, y  = cosine)) +
  geom_point() +
  lims(y = c(0.5, 1)) +
  geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
refit_cosine = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, br_nmf_res6$reconstructed))
data.frame(
  cosine = refit_cosine,
  nmuts = colSums(denovo_mut_mat_br)
) %>%
  ggplot(aes(x = nmuts, y  = cosine)) +
  geom_point() +
  lims(y = c(0.5, 1)) +
  geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
refit_cosine = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, br_nmf_res7$reconstructed))
data.frame(
  cosine = refit_cosine,
  nmuts = colSums(denovo_mut_mat_br)
) %>%
  ggplot(aes(x = nmuts, y  = cosine)) +
  geom_point() +
  lims(y = c(0.5, 1)) +
  geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)

plot_cosine_heatmap(cos_sim_matrix(br_nmf_res6$signatures, br_nmf_res6$signatures))
plot_cosine_heatmap(cos_sim_matrix(br_nmf_res7$signatures, br_nmf_res7$signatures))
plot_cosine_heatmap(cos_sim_matrix(br_nmf_res3$signatures, br_nmf_res3$signatures))
pdf("./3_Output/FullDenovoMutPatterns/Br_CosineRefit/CosineSim_COSMIC.pdf")
plot_cosine_heatmap(cos_sim_matrix(signatures, br_nmf_res6$signatures), plot_values = T)
dev.off()

#colnames(br_nmf_res6$signatures)
#[1] "SBS87-like" "HSPC-like"  "SBS1-like"  "PTA-like"   "SBS29-like" "SBSA" 

# Original
refit_cosine = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, br_nmf_res6$reconstructed))
pdf("./3_Output/FullDenovoMutPatterns/Br_CosineRefit/br_nmf_res6.pdf")
data.frame(cosine = refit_cosine, nmuts = colSums(denovo_mut_mat_br)) %>%
  ggplot(aes(x = nmuts, y  = cosine)) + geom_point() + lims(y = c(0.5, 1)) + geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
dev.off()
# Expected muts  
replaced_sigs_1 <- cbind(signatures[ ,c('SBS1', 'SBS87', 'HSPC')], br_nmf_res6$signatures[ ,c(4,5,6)])
repl_contri_1 = fit_to_signatures_strict(mut_matrix = denovo_mut_mat_br, signatures = replaced_sigs_1, max_delta = 0.01)$fit_res$reconstructed
refit_cosine_1 = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, repl_contri_1))
pdf("./3_Output/FullDenovoMutPatterns/Br_CosineRefit/br_nmf_res6_expected.pdf")
data.frame(cosine = refit_cosine_1, nmuts = colSums(denovo_mut_mat_br)) %>%
  ggplot(aes(x = nmuts, y  = cosine)) + geom_point() + lims(y = c(0.5, 1)) + geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
dev.off()
# SBS 18 
replaced_sigs_2 <- cbind(signatures[ ,c('SBS1', 'SBS87', 'HSPC', "SBS18")], br_nmf_res6$signatures[ ,c(4,6)])
repl_contri_2 = fit_to_signatures_strict(mut_matrix = denovo_mut_mat_br, signatures = replaced_sigs_2, max_delta = 0.01)$fit_res$reconstructed
refit_cosine_2 = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, repl_contri_2))
pdf("./3_Output/FullDenovoMutPatterns/Br_CosineRefit/br_nmf_res6_expected_SBS18.pdf")
data.frame(cosine = refit_cosine_2, nmuts = colSums(denovo_mut_mat_br)) %>%
  ggplot(aes(x = nmuts, y  = cosine)) + geom_point() + lims(y = c(0.5, 1)) + geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
dev.off()
# SBS 29 
replaced_sigs_3 <- cbind(signatures[ ,c('SBS1', 'SBS87', 'HSPC', 'SBS29')], br_nmf_res6$signatures[ ,c(4,6)])
repl_contri_3 = fit_to_signatures_strict(mut_matrix = denovo_mut_mat_br, signatures = replaced_sigs_3, max_delta = 0.01)$fit_res$reconstructed
refit_cosine_3 = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, repl_contri_3))
pdf("./3_Output/FullDenovoMutPatterns/Br_CosineRefit/br_nmf_res6_expected_SBS29.pdf")
data.frame(cosine = refit_cosine_3, nmuts = colSums(denovo_mut_mat_br)) %>%
  ggplot(aes(x = nmuts, y  = cosine)) + geom_point() + lims(y = c(0.5, 1)) + geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)
dev.off()


replaced_sigs  = cbind(signatures[ ,c('SBS1', 'SBS87', 'HSPC')], br_nmf_res6$signatures[ ,c(4,5,6)])
repl_contri = fit_to_signatures_strict(mut_matrix = denovo_mut_mat_br, signatures = replaced_sigs, max_delta = 0.01)$fit_res$contribution
#repl_contri = fit_to_signatures_strict(mut_matrix = denovo_mut_mat_br, signatures = replaced_sigs, max_delta = 0.01)$fit_res$reconstructed
reconstructed_replace = replaced_sigs %*% repl_contri
refit_cosine = diag(cos_sim_matrix(mut_matrix1 = denovo_mut_mat_br, reconstructed_replace))
data.frame(
  cosine = refit_cosine,
  nmuts = colSums(denovo_mut_mat_br)
) %>%
  ggplot(aes(x = nmuts, y  = cosine)) +
  geom_point() +
  lims(y = c(0.5, 1)) +
  geom_hline(yintercept = 0.9) + geom_vline(xintercept = 100)

#colSums(1-(replaced_sigs_1/colSums(replaced_sigs_1)))

# Bootstrapped refit 
contri_boots_br_orig <- fit_to_signatures_bootstrapped(NoLowMuts_br, br_nmf_res6$signatures, n_boots = 50, method = "strict")
contri_boots_br_1 <- fit_to_signatures_bootstrapped(NoLowMuts_br, replaced_sigs_1, n_boots = 50, method = "strict")
contri_boots_br_2 <- fit_to_signatures_bootstrapped(NoLowMuts_br, replaced_sigs_2, n_boots = 50, method = "strict")
contri_boots_br_3 <- fit_to_signatures_bootstrapped(NoLowMuts_br, replaced_sigs_3, n_boots = 50, method = "strict")

plot_bootstrapped_contribution(contri_boots_br_orig)
plot_bootstrapped_contribution(contri_boots_br_1)
plot_bootstrapped_contribution(contri_boots_br_2)
plot_bootstrapped_contribution(contri_boots_br_3)

colSums(replaced_sigs_1)
plot_bootstrapped_contribution(contri_boots_br_orig[rownames(contri_boots_br_orig) %in% c(paste0(sample_name, "_", 1:50)),])
plot_bootstrapped_contribution(contri_boots_br_1[rownames(contri_boots_br_1) %in% c(paste0(sample_name, "_", 1:50)),])
plot_bootstrapped_contribution(contri_boots_br_2[rownames(contri_boots_br_2) %in% c(paste0(sample_name, "_", 1:50)),])
plot_bootstrapped_contribution(contri_boots_br_3[rownames(contri_boots_br_3) %in% c(paste0(sample_name, "_", 1:50)),])
####### Test refitting 


cust_signatures <- signatures[, colnames(signatures) %in% c("PTA", "SBS1", "SBS5", "SBS18", "SBS29", "SBS87", "HSPC") ]
cust_strict_refit <- fit_to_signatures_strict(denovo_mut_mat_br, cust_signatures, max_delta = 0.1)
cust_fit_res_strict <- cust_strict_refit$fit_res
#pdf("./3_Output/FullDenovoMutPatterns/custom_branchedRefit_MD0.1.pdf", width = 14, height = 14)
plot_contribution(cust_fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#dev.off()








# transcription strand bias ------------------ 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
seqlevels(genes, pruning.mode = 'tidy') = paste0('chr', c(1:22, "X", "Y"))
genes = sortSeqlevels(genes)
genes = sort(genes)
# do strand bias
br_grl <- c(pt2322_br_grl, pt2229_br_grl, pt2283D_br_grl, pt344_br_grl, pt3291_br_grl, pt1180_br_grl, pt1909_br_grl)
br_grl_snv <- get_mut_type(br_grl, "snv")
br_mut_mat_s <- mut_matrix_stranded(vcf_list = br_grl_snv, ref_genome, genes)
# plot profile
plot_192_profile(br_mut_mat_s)
# do enrichment analysis
br_strand_counts <- strand_occurrences(br_mut_mat_s, by = names(br_grl_snv))
br_strand_bias <- strand_bias_test(br_strand_counts)
br_strand_bias %>% arrange(fdr)
cowplot::plot_grid(
  plot_strand(br_strand_counts),
  plot_strand_bias(br_strand_bias),
  ncol = 1
)
colnames(br_strand_bias)
unique(br_strand_bias$group)
plot_strand_bias(br_strand_bias[br_strand_bias$group %in% c("pt1909_root_branch.vcf", "pt1909_TALL1_TALL5_TALL4_TALL9_TALL7_TALL6_TALL2_TALL3_branch.vcf"),])
plot_strand(br_strand_counts[br_strand_counts$group %in% c("pt1909_root_branch.vcf", "pt1909_TALL1_TALL5_TALL4_TALL9_TALL7_TALL6_TALL2_TALL3_branch.vcf"),])
###################################################### Branch part ######################################################



###################################################### Custom refit part ######################################################
# Create the correct order 
pt2322_names <- c("pt2322_ALLBULK", "pt2322_TALL1","pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4", "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", 
                  "pt2322_TALL9", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12")#, "root_branch", "ISP_DP_branch", "ISP_branch", "DN_DP_SP4_branch")
pt2229_names <- c("pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", "pt2229_TALL7", 
                  "pt2229_TALL8", "pt2229_TALL9", "pt2229_TALL10")#, "pt2229_root_branch", "pt2229_DN_ISP_branch", "pt2229_Other_branch")
pt2283D_names <- c("pt2283D-DX1BM-ALLBULK", "pt2283D-DX1BM-TALL1", "pt2283D-DX1BM-TALL2", "pt2283D-DX1BM-TALL3", "pt2283D-DX1BM-TALL4", 
                          "pt2283D-DX1BM-TALL5", "pt2283D-DX1BM-TALL6", "pt2283D-DX1BM-TALL7", "pt2283D-DX1BM-TALL8", "pt2283D-DX1BM-TALL9", 
                          "pt2283D-DX1BM-TALL10", "pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")
pt344_names <- c("pt344-DX1BM-ALLBULK", "pt344-DX1PB-TALL1", "pt344-DX1PB-TALL2", "pt344-DX1PB-TALL3", "pt344-DX1PB-TALL4", "pt344-DX1PB-TALL5", 
                        "pt344-DX1PB-TALL6", "pt344-DX1PB-TALL7", "pt344-DX1PB-TALL8", "pt344-DX1PB-TALL9")
pt3291_names <- c("pt3291-DX1PB-TALL1", "pt3291-DX1PB-TALL2", "pt3291-DX1PB-TALL3", "pt3291-DX1PB-TALL4", "pt3291-DX1BM-TALL5", "pt3291-DX1BM-TALL6", 
                  "pt3291-DX1BM-TALL7", "pt3291-DX1BM-TALL8", "pt3291-DX1BM-TALL9", "pt3291-DX1BM-TALL10", "pt3291-DX1BM-TALL11", "pt3291-DX1BM-TALL12")
pt1180_names <- c("pt1180-DX1PB-TALL1", "pt1180-DX1PB-TALL2", "pt1180-DX1PB-TALL3", "pt1180-DX1PB-TALL4", "pt1180-DX1PB-TALL5", "pt1180-DX1PB-TALL6", 
                  "pt1180-DX1PB-TALL7", "pt1180-DX1PB-TALL8", "pt1180-DX1PB-TALL9", "pt1180-DX1PB-TALL10", "pt1180-DX1PB-TALL11", "pt1180-DX1PB-TALL12")
pt1909_names <- c("pt1909-DX1BM-TALL1", "pt1909-DX1BM-TALL2", "pt1909-DX1BM-TALL3", "pt1909-DX1BM-TALL4", "pt1909-DX1BM-TALL5", "pt1909-DX1BM-TALL6", 
                  "pt1909-DX1BM-TALL7", "pt1909-DX1BM-TALL8", "pt1909-DX1BM-TALL9")

### Type of occurrences 
COLORS6 <- c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
COLORS7 <- c("#2EBAED", "#000000", "#E98C7B", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
type_occurrences <- mut_type_occurrences(c(pt2322_grl, pt2229_grl, pt2283D_grl, pt344_grl, pt3291_grl, pt1180_grl, pt1909_grl), ref_genome)

### Refit, and visualization 
# pt2322
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18")]
fit_res <- fit_to_signatures(mut_mat_own[, pt2322_names], sub.sig)
strict_refit <- fit_to_signatures_strict(mut_mat_own[, pt2322_names], sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_spectrum.pdf")
plot_spectrum_stacked(type_occurrences[pt2322_names,] , max_value_to_show = 10)
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt2229
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18")]
fit_res <- fit_to_signatures(mut_mat_own[, pt2229_names], sub.sig)
strict_refit <- fit_to_signatures_strict(mut_mat_own[, pt2229_names], sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_spectrum.pdf")
plot_spectrum_stacked(type_occurrences[pt2229_names,] , max_value_to_show = 10)
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt2283D
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18")]
fit_res <- fit_to_signatures(mut_mat_own[, pt2283D_names], sub.sig)
strict_refit <- fit_to_signatures_strict(mut_mat_own[, pt2283D_names], sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2283D/pt2283D_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2283D/pt2283D_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2283D/pt2283D_spectrum.pdf")
plot_spectrum_stacked(type_occurrences[pt2283D_names,] , max_value_to_show = 10)
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2283D/pt2283D_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2283D/pt2283D_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt344
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
fit_res <- fit_to_signatures(mut_mat_own[, pt344_names], sub.sig)
strict_refit <- fit_to_signatures_strict(mut_mat_own[, pt344_names], sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt344/pt344_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt344/pt344_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt344/pt344_spectrum.pdf")
plot_spectrum_stacked(type_occurrences[pt344_names,] , max_value_to_show = 10)
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt344/pt344_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt344/pt344_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt3291
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
custom_mat <- mut_mat_own[, pt3291_names]
custom_mat <- cbind(custom_mat, as.matrix(pt3291_bulks_mut_mat[,"pt3291-DX1PB-ALLBULK"]))
colnames(custom_mat) <- c(pt3291_names, "pt3291-DX1PB-ALLBULK")
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt3291/pt3291_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt3291/pt3291_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt3291/pt3291_spectrum.pdf")
plot_spectrum_stacked(type_occurrences[pt3291_names,] , max_value_to_show = 10)
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt3291/pt3291_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt3291/pt3291_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt1180
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
custom_mat <- mut_mat_own[, pt1180_names]
custom_mat <- cbind(custom_mat, as.matrix(pt1180_bulks_mut_mat[,"pt1180-DX1PB-ALLBULK"]))
colnames(custom_mat) <- c(pt1180_names, "pt1180-DX1PB-ALLBULK")
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1180/pt1180_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1180/pt1180_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1180/pt1180_spectrum.pdf")
plot_spectrum_stacked(type_occurrences[pt1180_names,] , max_value_to_show = 10)
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1180/pt1180_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1180/pt1180_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt1909
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
custom_mat <- mut_mat_own[, pt1909_names]
custom_mat <- cbind(custom_mat, as.matrix(pt1909_bulks_mut_mat[,"pt1909-DX1BM-ALLBULK"]))
colnames(custom_mat) <- c(pt1909_names, "pt1909-DX1BM-ALLBULK")
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1909/pt1909_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1909/pt1909_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1909/pt1909_spectrum.pdf")
plot_spectrum_stacked(type_occurrences[pt1909_names,] , max_value_to_show = 10)
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1909/pt1909_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt1909/pt1909_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

################### Including SBS29
# pt344
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS29", "SBS87")]
fit_res <- fit_to_signatures(mut_mat_own[, pt344_names], sub.sig)
strict_refit <- fit_to_signatures_strict(mut_mat_own[, pt344_names], sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt344/pt344_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt344/pt344_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt344/pt344_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt344/pt344_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt3291
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS29", "SBS87")]
custom_mat <- mut_mat_own[, pt3291_names]
custom_mat <- cbind(custom_mat, as.matrix(pt3291_bulks_mut_mat[,"pt3291-DX1PB-ALLBULK"]))
colnames(custom_mat) <- c(pt3291_names, "pt3291-DX1PB-ALLBULK")
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt3291/pt3291_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt3291/pt3291_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt3291/pt3291_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt3291/pt3291_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt1180
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS29", "SBS87")]
custom_mat <- mut_mat_own[, pt1180_names]
custom_mat <- cbind(custom_mat, as.matrix(pt1180_bulks_mut_mat[,"pt1180-DX1PB-ALLBULK"]))
colnames(custom_mat) <- c(pt1180_names, "pt1180-DX1PB-ALLBULK")
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1180/pt1180_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1180/pt1180_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1180/pt1180_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1180/pt1180_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# pt1909
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS29", "SBS87")]
custom_mat <- mut_mat_own[, pt1909_names]
custom_mat <- cbind(custom_mat, as.matrix(pt1909_bulks_mut_mat[,"pt1909-DX1BM-ALLBULK"]))
colnames(custom_mat) <- c(pt1909_names, "pt1909-DX1BM-ALLBULK")
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1909/pt1909_refit_absolute.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1909/pt1909_refit_relative.pdf")
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1909/pt1909_strict_refit_absolute.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/InclSBS29/pt1909/pt1909_strict_refit_relative.pdf")
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

###################################################### Custom refit part ######################################################


###################################################### Custom branch refit part ######################################################
#mut_mat_br_own
#denovo_mut_mat_br

### Type of occurrences 
COLORS6 <- c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
COLORS7 <- c("#2EBAED", "#000000", "#E98C7B", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")


# pt1909
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
custom_mat <- mut_mat_br_own[, pt1909_branch_names]
custom_mat <- custom_mat[, !colnames(custom_mat) %in% c("pt1909_TALL1_TALL5_TALL4_branch.vcf") ]
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit//pt1909/pt1909_refit_absolute.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt1909/pt1909_refit_relative.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt1909/pt1909_strict_refit_absolute.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt1909/pt1909_strict_refit_relative.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# pt1180
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
custom_mat <- mut_mat_br_own[, pt1180_branch_names]
custom_mat <- custom_mat[, !colnames(custom_mat) %in% c("pt1180_TALL4_TALL12_branch.vcf", "pt1180_header.vcf") ]
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit//pt1180/pt1180_refit_absolute.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt1180/pt1180_refit_relative.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt1180/pt1180_strict_refit_absolute.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt1180/pt1180_strict_refit_relative.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# pt344
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
custom_mat <- mut_mat_br_own[, pt344_branch_names]
custom_mat <- custom_mat[, !colnames(custom_mat) %in% c("pt344_TALL1_TALL5_branch.vcf") ]
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit//pt344/pt344_refit_absolute.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt344/pt344_refit_relative.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt344/pt344_strict_refit_absolute.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt344/pt344_strict_refit_relative.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# pt3291
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS87")]
custom_mat <- mut_mat_br_own[, pt3291_branch_names]
custom_mat <- custom_mat[, !colnames(custom_mat) %in% c("pt3291_TALL9_TALL6_TALL8_TALL11_TALL3_branch.vcf") ]
fit_res <- fit_to_signatures(custom_mat, sub.sig)
strict_refit <- fit_to_signatures_strict(custom_mat, sub.sig, max_delta = 0.001)
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit//pt3291/pt3291_refit_absolute.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt3291/pt3291_refit_relative.pdf", height = 10)
plot_contribution(fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt3291/pt3291_strict_refit_absolute.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/BranchRefit/pt3291/pt3291_strict_refit_relative.pdf", height = 10)
plot_contribution(strict_refit$fit_res$contribution, coord_flip = FALSE,mode = "relative"
) + scale_fill_manual(values=mycols_paired) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
###################################################### Custom branch refit part ######################################################








indel_signatures = get_known_signatures(muttype = "indel")

indel_cont = get_indel_context(grl_indel, ref_genome)
indel_count = count_indel_contexts(indel_cont)







###################################################### Old part ######################################################

#colSums(pooled_pt2322_branch_mat)
#pt2322_branch <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
#pooled_mut_mat <- pool_mut_mat(mut_mat, grouping = tissue)


### Branched matrices 
#pt2322
pt2322_branch_grl <- read_vcfs_as_granges(pt2322_vcf_branches[2:24], pt2322_branch_names, ref_genome)
pt2322_branch_mat <- mut_matrix(vcf_list = pt2322_branch_grl, ref_genome = ref_genome)
pooled_pt2322_branch_mat <- pool_mut_mat(pt2322_branch_mat, grouping = pt2322_branch)
#pt2229
pt2229_branch_grl <- read_vcfs_as_granges(pt2229_vcf_branches, pt2229_branch_names, ref_genome)
pt2229_branch_mat <- mut_matrix(vcf_list = pt2229_branch_grl, ref_genome = ref_genome)
pooled_pt2229_branch_mat <- pool_mut_mat(pt2229_branch_mat, grouping = pt2229_branch)






####################################################### subset signatures 
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18")]
fit_res <- fit_to_signatures(mut_mat, sub.sig)
pdf("./3_Output/FullDenovoMutPatterns/Refit/Default_Refit.pdf")
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### I added this 
contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
                                               signatures,
                                               n_boots = 100,
                                               method = "strict"
)
contri_boots_sub <- fit_to_signatures_bootstrapped(mut_mat,
                                               sub.sig,
                                               n_boots = 100,
                                               method = "strict"
)
pdf("./3_Output/FullDenovoMutPatterns/Refit/BootstrappedContributions_n100_Sub.pdf", height = 25)
plot_bootstrapped_contribution(contri_boots_sub)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Refit/BootstrappedContributionsRel_n100_sub.pdf", width = 15)
plot_bootstrapped_contribution(contri_boots_sub, 
                               mode = "relative", 
                               plot_type = "dotplot")
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Refit/BootstrappedContributions_n100.pdf", height = 25)
plot_bootstrapped_contribution(contri_boots)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/Refit/BootstrappedContributionsRel_n100.pdf", width = 15)
plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "dotplot")
dev.off()
### End here 
pdf("./3_Output/FullDenovoMutPatterns/Refit/Default_Refit_relative.pdf")
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

write.table(fit_res$contribution, file = "./3_Output/FullDenovoMutPatterns/Refit/Refit_Contributions.txt", sep = "\t", col.names = NA)

# Not usefull
contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
                                               sub.sig,
                                               n_boots = 50,
                                               method = "strict"
)
plot_bootstrapped_contribution(contri_boots)
plot_correlation_bootstrap(contri_boots)




cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, sub.sig)
plot_cosine_heatmap(cos_sim_samples_signatures, 
                    cluster_rows = TRUE, cluster_cols = TRUE)
cos_sim_samples <- cos_sim_matrix(mut_mat, mut_mat)
jpeg("./3_Output/FullDenovoMutPatterns/Refit/cosine_sim_Samples.jpeg")
plot_cosine_heatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()

p1 <- plot_cosine_heatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE)
ggsave("./3_Output/FullDenovoMutPatterns/Refit/cosine_sim_Samples.pdf", p1)




strict_refit <- fit_to_signatures_strict(mut_mat, sub.sig, max_delta = 0.004)
fit_res_strict <- strict_refit$fit_res
pdf("./3_Output/FullDenovoMutPatterns/Refit/Strict_Refit.pdf")
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("./3_Output/FullDenovoMutPatterns/Refit/Strict_Refit_relative.pdf")
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
write.table(fit_res_strict$contribution, file = "./3_Output/FullDenovoMutPatterns/Refit/Strict_Refit_Contributions.txt", sep = "\t", col.names = NA)




#plot_spectrum_stacked(fit_res$contribution)
#plot_contribution2()
plot_contribution2(fit_res$contribution,
                   coord_flip = FALSE,
                   mode = "absolute"
) #+ 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_contribution(fit_res$contribution[,1:17],
                   coord_flip = FALSE,
                   mode = "absolute"
) + scale_fill_manual(values=mycols_dark2) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_contribution(fit_res$contribution[,1:17],
                  coord_flip = FALSE,
                  mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


plot_contribution(fit_res$contribution[,1:17],
                  coord_flip = FALSE,
                  mode = "relative"
) + scale_fill_manual(values=palette("Dark2")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


palette("Dark2")


round(colSums(fit_res$contribution[,1:17]))
colSums(mut_mat[,1:17])


test_fit <- fit_res$contribution[,1:17]
colnames(test_fit)


testnames <- c("pt2322_ALLBULK", "pt2322_TALL1","pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4", "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", 
               "pt2322_TALL9", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12", 
               "root_branch", "ISP_DP_branch", "ISP_branch", "DN_DP_SP4_branch")
test_fit[,testnames]

pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_default_refit_colors.pdf")
plot_contribution(test_fit[,testnames],
                  coord_flip = FALSE,
                  mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_default_refit_colors_relative.pdf")
plot_contribution(test_fit[,testnames],
                  coord_flip = FALSE,
                  mode = "relative"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


grl_list <- c(pt2322_grl, pt2322_branch_grl)#, pt2229_mut_mat, pooled_pt2229_branch_mat)
type_occurrences <- mut_type_occurrences(grl_list, ref_genome)
type_occurrences
p1 <- plot_spectrum(type_occurrences)
p1
p2 <- plot_spectrum(type_occurrences, by = rownames(type_occurrences), legend = TRUE)


COLORS6 <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
)
COLORS7 <- c(
  "#2EBAED", "#000000", 
  "#E98C7B", "#DE1C14","#D4D2D2", "#ADCC54",
  "#F0D0CE"
)
plot_spectrum_stacked(type_occurrences)
type_occurrences[testnames,]
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_StackedSpectrum.pdf")
plot_spectrum_stacked(type_occurrences[c("pt2322_ALLBULK", "pt2322_TALL1","pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4", "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", 
                                         "pt2322_TALL9", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12"),] , max_value_to_show = 10)
dev.off()


colnames(fit_res$contribution)
pt2229_fit <- fit_res$contribution[,c("pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", "pt2229_TALL7", 
                                     "pt2229_TALL8", "pt2229_TALL9", "pt2229_TALL10", "pt2229_root_branch", "pt2229_DN_ISP_branch", "pt2229_Other_branch")]
### Now for pt2229
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_default_refit_colors.pdf")
plot_contribution(pt2229_fit,
                  coord_flip = FALSE,
                  mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_default_refit_colors_relative.pdf")
plot_contribution(pt2229_fit,
                  coord_flip = FALSE,
                  mode = "relative"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
grl_list <- c(pt2229_grl, pt2229_branch_grl)#, pt2229_mut_mat, pooled_pt2229_branch_mat)
type_occurrences <- mut_type_occurrences(grl_list, ref_genome)
type_occurrences
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_StackedSpectrum.pdf")
plot_spectrum_stacked(type_occurrences[c("pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", "pt2229_TALL7", 
                                         "pt2229_TALL8", "pt2229_TALL9", "pt2229_TALL10"),] , max_value_to_show = 10)
dev.off()




##mut_matrix(signatures[,"PTA"])

#plot_96_profile(as.matrix(signatures[,"PTA"]))
#plot_spectrum(as.matrix(signatures[,"PTA"]))

pta_vis <- as.matrix(signatures[, "PTA"])
colnames(pta_vis) <- "PTA"
rownames(pta_vis) <- rownames(pt2322_branch_mat)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Signatures/PTA_signature.pdf")
plot_96_profile(pta_vis)
dev.off()

pta_vis <- as.matrix(signatures[, "HSPC"])
colnames(pta_vis) <- "HSPC"
rownames(pta_vis) <- rownames(pt2322_branch_mat)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Signatures/HSPC_signature.pdf")
plot_96_profile(pta_vis)
dev.off()

pta_vis <- as.matrix(signatures[, "SBS1"])
colnames(pta_vis) <- "SBS1"
rownames(pta_vis) <- rownames(pt2322_branch_mat)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Signatures/SBS1_signature.pdf")
plot_96_profile(pta_vis)
dev.off()

pta_vis <- as.matrix(signatures[, "SBS5"])
colnames(pta_vis) <- "SBS5"
rownames(pta_vis) <- rownames(pt2322_branch_mat)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Signatures/SBS5_signature.pdf")
plot_96_profile(pta_vis)
dev.off()

pta_vis <- as.matrix(signatures[, "SBS18"])
colnames(pta_vis) <- "SBS18"
rownames(pta_vis) <- rownames(pt2322_branch_mat)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Signatures/SBS18_signature.pdf")
plot_96_profile(pta_vis)
dev.off()



