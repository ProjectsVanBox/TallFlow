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

################################################################ Walker end branch ################################################################
### Gather the files 
High_vcf_files <- list.files(path = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/1_Input/pt3291/High/",
                               pattern = "*High.vcf$", full.names = TRUE, recursive = TRUE )
Intermediate_vcf_files <- list.files(path = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/1_Input/pt3291/Intermediate/",
                             pattern = "*Intermediate.vcf$", full.names = TRUE, recursive = TRUE )
Low_vcf_files <- list.files(path = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/1_Input/pt3291/Low/",
                             pattern = "*Low.vcf$", full.names = TRUE, recursive = TRUE )

### Sample names 
High_names <- gsub(pattern = ".walker.bed_High.vcf", replacement = "", 
     x = gsub(pattern = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/1_Input/pt3291/High//231016_HMFreg2090_pt3291.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_", 
     replacement = "", High_vcf_files))
Intermediate_names <- gsub(pattern = ".walker.bed_Intermediate.vcf", replacement = "", 
                   x = gsub(pattern = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/1_Input/pt3291/Intermediate//231016_HMFreg2090_pt3291.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_", 
                            replacement = "", Intermediate_vcf_files))
Low_names <- gsub(pattern = ".walker.bed_Low.vcf", replacement = "", 
                   x = gsub(pattern = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/1_Input/pt3291/Low//231016_HMFreg2090_pt3291.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_", 
                            replacement = "", Low_vcf_files))

# High
High_grl <- read_vcfs_as_granges(High_vcf_files, High_names, ref_genome, type = 'all')
High_mut_mat <- mut_matrix(vcf_list = get_mut_type(High_grl, 'snv'), ref_genome = ref_genome)
# Intermediate
Intermediate_grl <- read_vcfs_as_granges(Intermediate_vcf_files, Intermediate_names, ref_genome, type = 'all')
Intermediate_mut_mat <- mut_matrix(vcf_list = get_mut_type(Intermediate_grl, 'snv'), ref_genome = ref_genome)
# Low
Low_grl <- read_vcfs_as_granges(Low_vcf_files, Low_names, ref_genome, type = 'all')
Low_mut_mat <- mut_matrix(vcf_list = get_mut_type(Low_grl, 'snv'), ref_genome = ref_genome)
# Combined
Comb_mut_mat <- cbind(High_mut_mat, Intermediate_mut_mat, Low_mut_mat)
Comb_Group <- c(rep.int(x = "High", times = length(colnames(High_mut_mat))), 
                rep.int(x = "Intermediate", times = length(colnames(Intermediate_mut_mat))),
                rep.int(x = "Low", times = length(colnames(Low_mut_mat))) )


### Check signatures
pooled_Comb_mat <- pool_mut_mat(Comb_mut_mat, grouping = Comb_Group)
Comb_strict_refit <- fit_to_signatures_strict(pooled_Comb_mat, signatures, max_delta = 0.004)
sub_sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS29", "SBS87")]
Comb_strict_refit_sub <- fit_to_signatures_strict(pooled_Comb_mat, sub_sig, max_delta = 0.004)

### Plot
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/WalkerScoreEndBranch/96_profile.pdf")
plot_96_profile(pooled_Comb_mat)
dev.off()
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/WalkerScoreEndBranch/absolute_contibution_allSig.pdf")
plot_contribution(Comb_strict_refit$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/WalkerScoreEndBranch/absolute_contibution_subSig.pdf")
plot_contribution(Comb_strict_refit_sub$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
################################################################ Walker end branch ################################################################






################################################################ Per cell analysis ################################################################
### Gather the files 

#Vera:
pt3291_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt3291/snvs/pt3291/",
                             pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt3291_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = pt3291_files))
pt1909_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt1909/snvs/pt1909/",
                           pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt1909_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = pt1909_files))
PB30602_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Tsites/3_Output/PTATO/PB30602/snvs/PB30602/",
                           pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
PB30602_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = PB30602_files))
# V1 
pt344_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt344/snvs/pt344/",
                              pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt344_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = pt344_files))
pt1180_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt1180/snvs/pt1180/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt1180_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = pt1180_files))
# Anais:
PB30885_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/SteMDS/3_Output/PTATO/PB30885/all_ptato_filtered/",
                            pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
PB30885_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = PB30885_files))
# V1 
PB32292_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/SteMDS/3_Output/PTATO/PB32292/snvs/PB32292/all.ptato.filtered/",
                            pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
PB32292_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = PB32292_files))
# Alexander:
PB14458_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458/snvs/PB14458/",
                           pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
PB14458_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = PB14458_files))
PB11197_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB11197/snvs/PB11197/",
                            pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
PB11197_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                             x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = PB11197_files))
#/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB08410_new V1
#/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB11197 V2 (the new PB11197 samples that we are currently processing are actually V1 but they should not be in the ptato_vcf file yet, so for now all the PB11197 are V2)
#/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458 V2 

# Lucca/Jurrian:
CB29_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/IBFM_tAML/3_Output/PTATO_133/STEMCELLCB29/snvs/STEMCELLCB29/",
                            pattern = ".*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
CB29_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                            x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = CB29_files))

# Joske/Liza
Mmorganni_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/STE0072_Mmorganni/3_Output/PTATO/STE0072_Mmorganni/snvs/STE0072_Mmorganni/",
                          pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
Mmorganni_sample_names <- gsub(pattern = ".snvs.ptato.filtered.vcf.gz", replacement = "", 
                           x = gsub(pattern = ".*hg38.sorted_", replacement = "", x = Mmorganni_files))

### Create the matrices 
# V2
pt3291_grl <- read_vcfs_as_granges(pt3291_files, pt3291_sample_names, ref_genome, type = 'all')
pt3291_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt3291_grl, 'snv'), ref_genome = ref_genome)
pt1909_grl <- read_vcfs_as_granges(pt1909_files, pt1909_sample_names, ref_genome, type = 'all')
pt1909_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1909_grl, 'snv'), ref_genome = ref_genome)
PB30602_grl <- read_vcfs_as_granges(PB30602_files, PB30602_sample_names, ref_genome, type = 'all')
PB30602_mut_mat <- mut_matrix(vcf_list = get_mut_type(PB30602_grl, 'snv'), ref_genome = ref_genome)
PB30885_grl <- read_vcfs_as_granges(PB30885_files, PB30885_sample_names, ref_genome, type = 'all')
PB30885_mut_mat <- mut_matrix(vcf_list = get_mut_type(PB30885_grl, 'snv'), ref_genome = ref_genome)
PB32292_grl <- read_vcfs_as_granges(PB32292_files, PB32292_sample_names, ref_genome, type = 'all')
PB32292_mut_mat <- mut_matrix(vcf_list = get_mut_type(PB32292_grl, 'snv'), ref_genome = ref_genome)
PB14458_grl <- read_vcfs_as_granges(PB14458_files, PB14458_sample_names, ref_genome, type = 'all')
PB14458_mut_mat <- mut_matrix(vcf_list = get_mut_type(PB14458_grl, 'snv'), ref_genome = ref_genome)
PB11197_grl <- read_vcfs_as_granges(PB11197_files, PB11197_sample_names, ref_genome, type = 'all')
PB11197_mut_mat <- mut_matrix(vcf_list = get_mut_type(PB11197_grl, 'snv'), ref_genome = ref_genome)
CB29_grl <- read_vcfs_as_granges(CB29_files, CB29_sample_names, ref_genome, type = 'all')
CB29_mut_mat <- mut_matrix(vcf_list = get_mut_type(CB29_grl, 'snv'), ref_genome = ref_genome)
Mmorganni_grl <- read_vcfs_as_granges(Mmorganni_files, Mmorganni_sample_names, ref_genome, type = 'all')
Mmorganni_mut_mat <- mut_matrix(vcf_list = get_mut_type(Mmorganni_grl, 'snv'), ref_genome = ref_genome)
# V1 
pt344_grl <- read_vcfs_as_granges(pt344_files, pt344_sample_names, ref_genome, type = 'all')
pt344_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt344_grl, 'snv'), ref_genome = ref_genome)
pt1180_grl <- read_vcfs_as_granges(pt1180_files, pt1180_sample_names, ref_genome, type = 'all')
pt1180_mut_mat <- mut_matrix(vcf_list = get_mut_type(pt1180_grl, 'snv'), ref_genome = ref_genome)
# Combine all
cells_mut_mat <- cbind(pt3291_mut_mat, pt1909_mut_mat, PB30602_mut_mat, PB30885_mut_mat, PB14458_mut_mat, PB11197_mut_mat, CB29_mut_mat, Mmorganni_mut_mat, 
                       pt344_mut_mat, pt1180_mut_mat, PB32292_mut_mat)


#plot_96_profile(cells_mut_mat) # Too many samples 
### Check signatures 
sub_sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18", "SBS29", "SBS87")]
Cells_strict_refit_sub <- fit_to_signatures_strict(cells_mut_mat, sub_sig, max_delta = 0.004)
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/AllCells_strictrefit_absolute.pdf", width = 20)
plot_contribution(Cells_strict_refit_sub$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/AllCells_strictrefit_relative.pdf", width = 20)
plot_contribution(Cells_strict_refit_sub$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/DxRel_couples.pdf", width = 20)
plot_contribution(Cells_strict_refit_sub$fit_res$contribution[, c(pt344_sample_names, pt3291_sample_names, pt1180_sample_names, pt1909_sample_names)],
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


### Pool samples 
pooled_samples <- c(rep.int(x = "V2_pt3291", length(colnames(pt3291_mut_mat))),
                    rep.int(x = "V2_pt1909", length(colnames(pt1909_mut_mat))),
                    rep.int(x = "V2_PB30602", length(colnames(PB30602_mut_mat))),
                    rep.int(x = "V2_PB30885", length(colnames(PB30885_mut_mat))),
                    rep.int(x = "V2_PB14458", length(colnames(PB14458_mut_mat))),
                    rep.int(x = "V2_PB11197", length(colnames(PB11197_mut_mat))),
                    rep.int(x = "V2_CB29", length(colnames(CB29_mut_mat))),
                    rep.int(x = "V2_Mmorganni", length(colnames(Mmorganni_mut_mat))),
                    rep.int(x = "V1_pt344", length(colnames(pt344_mut_mat))),
                    rep.int(x = "V1_pt1180", length(colnames(pt1180_mut_mat))),
                    rep.int(x = "V1_PB32292", length(colnames(PB32292_mut_mat))))
pooled_cells_mut_mat <- pool_mut_mat(cells_mut_mat, grouping = pooled_samples)
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/Pooled_96_profile.pdf", height = 10)
plot_96_profile(pooled_cells_mut_mat, ymax = 0.15)
dev.off()
pooled_cells_strict_refit_sub <- fit_to_signatures_strict(pooled_cells_mut_mat, sub_sig, max_delta = 0.004)
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/Pooled_strictrefit_absolute.pdf")
plot_contribution(pooled_cells_strict_refit_sub$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/Pooled_strictrefit_relative.pdf")
plot_contribution(pooled_cells_strict_refit_sub$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


### Check filtered out PTATO stuff
V1_filt_grl <-  read_vcfs_as_granges("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/vcf/isec_V1/0000.vcf", "V1_Filtered", ref_genome, type = 'all')
V1_filt_mut_mat <- mut_matrix(vcf_list = get_mut_type(V1_filt_grl, 'snv'), ref_genome = ref_genome)
V2_filt_grl <-  read_vcfs_as_granges("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/CheckArtefactsPTA/vcf/isec_V2/0001.vcf", "V2_Filtered", ref_genome, type = 'all')
V2_filt_mut_mat <- mut_matrix(vcf_list = get_mut_type(V2_filt_grl, 'snv'), ref_genome = ref_genome)



filt_strict_refit_sub <- fit_to_signatures_strict(cbind(V2_filt_mut_mat, V1_filt_mut_mat), sub_sig, max_delta = 0.004)
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/PTA_Filtered_1cell_absolute.pdf")
plot_contribution(filt_strict_refit_sub$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/PTA_Filtered_1cell_relative.pdf")
plot_contribution(filt_strict_refit_sub$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()






plot_contribution(Cells_strict_refit_sub$fit_res$contribution[, c(CB29_sample_names)],
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_contribution(Cells_strict_refit_sub$fit_res$contribution[, c(Mmorganni_sample_names)],
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






# Combine all
cells_mut_mat <- cbind(pt3291_mut_mat, pt1909_mut_mat, PB30602_mut_mat, PB30885_mut_mat, PB14458_mut_mat, PB11197_mut_mat, CB29_mut_mat, Mmorganni_mut_mat, 
                       pt344_mut_mat, pt1180_mut_mat, PB32292_mut_mat)
cells_mut_mat_renamed <- cells_mut_mat
renamed_names <- c(rep.int(x = "V2_Kit1_S1", length(colnames(pt3291_mut_mat))),
                   rep.int(x = "V2_Kit1_S2", length(colnames(pt1909_mut_mat))),
                   rep.int(x = "V2_Kit1_S3", length(colnames(PB30602_mut_mat))),
                   rep.int(x = "V2_Kit1_S4", length(colnames(PB30885_mut_mat))),
                   rep.int(x = "V2_Kit1+2", length(colnames(PB14458_mut_mat))),
                   rep.int(x = "V2_NA", length(colnames(PB11197_mut_mat))),
                   rep.int(x = "V2_Kit4_S1", length(colnames(CB29_mut_mat))),
                   rep.int(x = "V2_Kit1_S5", length(colnames(Mmorganni_mut_mat))),
                   rep.int(x = "V1_1", length(colnames(pt344_mut_mat))),
                   rep.int(x = "V1_2", length(colnames(pt1180_mut_mat))),
                   rep.int(x = "V1_3", length(colnames(PB32292_mut_mat))))
colnames(cells_mut_mat_renamed) <- renamed_names
#

Cells_renamed_strict_refit_sub <- fit_to_signatures_strict(cells_mut_mat_renamed, sub_sig, max_delta = 0.004)
bioskryp_contribution <- Cells_renamed_strict_refit_sub$fit_res$contribution
colnames(bioskryp_contribution) <- renamed_names
pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/bioskryp_contribution.pdf", width = 20)
plot_contribution(bioskryp_contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
renamed_pooled_cells_mut_mat <- pool_mut_mat(cells_mut_mat, grouping = renamed_names)

pdf("~/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/CheckArtefactsPTA/bioskryp_96_profile.pdf", height = 10)
plot_96_profile(renamed_pooled_cells_mut_mat, ymax = 0.15)
dev.off()
################################################################ Per cell analysis ################################################################

### Order presentation:
#Vera DxRel
# small TallFlow introduction
# small PTATO introduction
#--> What does this mean?
#--> V2
#talk:
#--> Other samples
#What to do now?

#- PTA Kit?
#- changes in work-up?
#- Lab differences?












