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


### Set and create directories 
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/")


### Gather the files 
pt2322_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt2322_vcf_branches <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/Branches_MutPat/pt2322/",
                                  pattern = "*.vcf", full.names = TRUE, recursive = TRUE )

#pt2322_bulk_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/BulkSeqs/",
#                                    pattern = "*filtered.vcf$", full.names = TRUE, recursive = TRUE )
pt2229_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2229/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt2229_vcf_branches <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/Branches_MutPat/pt2229/",
                                  pattern = "*.vcf", full.names = TRUE, recursive = TRUE )
#pt2229_bulk_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2229/BulkSeqs/",
#                                    pattern = "*filtered.vcf$", full.names = TRUE, recursive = TRUE )
#combined.files <- c(pt2322_vcf_files, pt2322_bulk_vcf_files[1:5], pt2229_vcf_files, pt2229_bulk_vcf_files)


### Sample names 
pt2322_sample_names <- c("pt2322_ALLBULK", "pt2322_TALL1", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12", "pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4",
                         "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", "pt2322_TALL9") 
pt2229_sample_names <- c("pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL10", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", 
                  "pt2229_TALL7", "pt2229_TALL8", "pt2229_TALL9")
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

pt2322_branch <- c("root_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", 
                   "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "ISP_branch", "ISP_branch", 
                   "ISP_branch", "ISP_branch", "ISP_branch", "ISP_DP_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch", "ISP_DP_branch", "DN_DP_SP4_branch", "DN_DP_SP4_branch")
pt2229_branch <- c("pt2229_root_branch", "pt2229_DN_ISP_branch", "pt2229_Other_branch", "pt2229_Other_branch", "pt2229_Other_branch", "pt2229_Other_branch", 
                   "pt2229_Other_branch", "pt2229_DN_ISP_branch", "pt2229_Other_branch")


### Create matrices 
#pt2322
pt2322_grl <- read_vcfs_as_granges(pt2322_vcf_files, pt2322_sample_names, ref_genome)
pt2322_mut_mat <- mut_matrix(vcf_list = pt2322_grl, ref_genome = ref_genome)
pt2322_branch_grl <- read_vcfs_as_granges(pt2322_vcf_branches[2:24], pt2322_branch_names, ref_genome)
pt2322_branch_mat <- mut_matrix(vcf_list = pt2322_branch_grl, ref_genome = ref_genome)
pooled_pt2322_branch_mat <- pool_mut_mat(pt2322_branch_mat, grouping = pt2322_branch)
#pt2229
pt2229_grl <- read_vcfs_as_granges(pt2229_vcf_files, pt2229_sample_names, ref_genome)
pt2229_mut_mat <- mut_matrix(vcf_list = pt2229_grl, ref_genome = ref_genome)
pt2229_branch_grl <- read_vcfs_as_granges(pt2229_vcf_branches, pt2229_branch_names, ref_genome)
pt2229_branch_mat <- mut_matrix(vcf_list = pt2229_branch_grl, ref_genome = ref_genome)
pooled_pt2229_branch_mat <- pool_mut_mat(pt2229_branch_mat, grouping = pt2229_branch)
#colSums(pooled_pt2322_branch_mat)
#pt2322_branch <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
#pooled_mut_mat <- pool_mut_mat(mut_mat, grouping = tissue)




mut_mat <- cbind(pt2322_mut_mat, pooled_pt2322_branch_mat, pt2229_mut_mat, pooled_pt2229_branch_mat)
### De novo mutational signature extraction using NMF
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", 
                nrun = 50, seed = 123456, .opt = "v-p")
pdf("./3_Output/FullDenovoMutPatterns/OnlyPatient/OnlyPatient_estimates.pdf")
plot(estimate)
dev.off()

### Extract signatures 
nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 50, single_core = TRUE)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C")


### Check with known signatures 
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
colnames(nmf_res$signatures)
pdf("./3_Output/FullDenovoMutPatterns/OnlyPatient/OnlyPatient_tri_nuc_profiles.pdf")
plot_96_profile(nmf_res$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/OnlyPatient/OnlyPatient_signature_contribution.pdf", width = 12)
plot_contribution(nmf_res$contribution, nmf_res$signature,
                  mode = "relative"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/OnlyPatient/OnlyPatient_original_vs_reconstructed.pdf", width = 12)
plot_original_vs_reconstructed(mut_mat, nmf_res$reconstructed, 
                               y_intercept = 0.95)
dev.off()











### Create mutational matrix  
grl <- read_vcfs_as_granges(combined.files, sample_names, ref_genome)
nmf_ref_mat = read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/signatures/NMF_reference_matrix.txt", sep = "\t", header = T)
mut_mat_own <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat <- cbind(mut_mat_own, nmf_ref_mat)


### De novo mutational signature extraction using NMF
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", 
                nrun = 50, seed = 123456, .opt = "v-p")
pdf("./3_Output/FullDenovoMutPatterns/estimates.pdf")
plot(estimate)
dev.off()

### Extract signatures 
nmf_res <- extract_signatures(mut_mat, rank = 4, nrun = 50, single_core = TRUE)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")


### Check with known signatures 
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
colnames(nmf_res$signatures)
pdf("./3_Output/FullDenovoMutPatterns/tri_nuc_profiles.pdf")
plot_96_profile(nmf_res$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/FullDenovoMutPatterns/signature_contribution.pdf", width = 12)
plot_contribution(nmf_res$contribution, nmf_res$signature,
                  mode = "relative"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



pdf("./3_Output/FullDenovoMutPatterns/original_vs_reconstructed.pdf", width = 12)
plot_original_vs_reconstructed(mut_mat, nmf_res$reconstructed, 
                               y_intercept = 0.95)
dev.off()







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



