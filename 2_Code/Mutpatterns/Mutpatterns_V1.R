#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(MutationalPatterns)
library(BSgenome)
library(ChIPpeakAnno)
library(ggplot2)
library("NMF")


### Read in genomes
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)


### Get signatures 
signatures = get_known_signatures()
pta_sig = read.table("/Users/ricohagelaar/Documents/Thymus_Project/PTATO/PTA_Artefact_Signature.txt", sep = "\t", header = T)
pta_sig = as.matrix(pta_sig)
PTA <- as.numeric(pta_sig[,"PTA"])
signatures <- cbind(PTA, signatures)


### Set and create directories 
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/")
dir.create("./3_Output/Mutpatterns")
dir.create("./3_Output/Mutpatterns/RefittingSignatures")
dir.create("./3_Output/Mutpatterns/RefittingSignatures/pt2322")
dir.create("./3_Output/Mutpatterns/RefittingSignatures/pt2229")
dir.create("./3_Output/Mutpatterns/96trinuc_profile")
dir.create("./3_Output/Mutpatterns/96trinuc_profile/pt2322")
dir.create("./3_Output/Mutpatterns/96trinuc_profile/pt2229")


### Read in data and create mutational matrix 
pt2322_sample_names <- c("ALLBULK", "TALL1", "TALL10", "TALL11", "TALL12", "TALL2", "TALL3", "TALL4", "TALL5", "TALL6", "TALL7", "TALL8", "TALL9", 
                         "ALLBULK-DN", "ALLBULK-DP", "ALLBULK-iSPCD4", "ALLBULK-SPCD4", "ALLBULK_2")#, "MONOBULK")
pt2322_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/",
                            pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE
)
pt2322_vcf_files <- c(pt2322_vcf_files[1:13], list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/BulkSeqs/",
                                                   pattern = "*filtered.vcf$", full.names = TRUE, recursive = TRUE
))
pt2322_grl <- read_vcfs_as_granges(pt2322_vcf_files[1:18], pt2322_sample_names, ref_genome)
pt2322_mut_mat <- mut_matrix(vcf_list = pt2322_grl, ref_genome = ref_genome)


### Print the tri-nucleotide plot per sample 
for (s in pt2322_sample_names){
  print(s)
  plotmat <- as.matrix(pt2322_mut_mat[, s])
  colnames(plotmat) <- s
  pdf(paste0("./3_Output/Mutpatterns/96trinuc_profile/pt2322/", s, "_trinuc.pdf"))
  print(plot_96_profile(plotmat))
  dev.off()
}


### De novo signature 
pt2322_mut_mat <- pt2322_mut_mat + 0.0001
estimate <- nmf(pt2322_mut_mat, rank = 1:8, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")
plot(estimate)


estimate2 <- nmf(pt2322_mut_mat, rank = 1:8, method = "brunet", 
                nrun = 50, seed = 123456, .opt = "v-p")
pdf("./3_Output/Mutpatterns/DeNovo/estimates.pdf")
plot(estimate2)
dev.off()



nmf_res <- extract_signatures(pt2322_mut_mat, rank = 4, nrun = 50, single_core = TRUE)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
# Check if signatures are similar, spoiler they are not
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
colnames(nmf_res$signatures)
pdf("./3_Output/Mutpatterns/DeNovo/tri_nuc_profiles.pdf")
plot_96_profile(nmf_res$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo/signature_contribution.pdf")
plot_contribution(nmf_res$contribution, nmf_res$signature,
                  mode = "relative"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo/plot_contribution_heatmap.pdf")
plot_contribution_heatmap(nmf_res$contribution, 
                          cluster_samples = TRUE, 
                          cluster_sigs = TRUE)
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo/plot_compare_profiles.pdf")
plot_compare_profiles(pt2322_mut_mat[, 1],
                      nmf_res$reconstructed[, 1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo/plot_original_vs_reconstructed.pdf")
plot_original_vs_reconstructed(pt2322_mut_mat, nmf_res$reconstructed, 
                               y_intercept = 0.95)
dev.off()



#### Now with 3 ranks 
nmf_res3 <- extract_signatures(pt2322_mut_mat, rank = 3, nrun = 50, single_core = TRUE)
colnames(nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
# Check if signatures are similar, spoiler they are not
nmf_res3 <- rename_nmf_signatures(nmf_res3, signatures, cutoff = 0.85)
colnames(nmf_res3$signatures)
pdf("./3_Output/Mutpatterns/DeNovo3/tri_nuc_profiles.pdf")
plot_96_profile(nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo3/signature_contribution.pdf")
plot_contribution(nmf_res3$contribution, nmf_res3$signature,
                  mode = "relative"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo3/plot_contribution_heatmap.pdf")
plot_contribution_heatmap(nmf_res3$contribution, 
                          cluster_samples = TRUE, 
                          cluster_sigs = TRUE)
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo3/plot_compare_profiles.pdf")
plot_compare_profiles(pt2322_mut_mat[, 1],
                      nmf_res3$reconstructed[, 1],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)
dev.off()
pdf("./3_Output/Mutpatterns/DeNovo3/plot_original_vs_reconstructed.pdf")
plot_original_vs_reconstructed(pt2322_mut_mat, nmf_res3$reconstructed, 
                               y_intercept = 0.95)
dev.off()





sub.sig <- signatures[, c("PTA", "SBS1", "SBS5", "SBS18")]
### extra2 strict refitting
pt2322_extra2strictfit  <- fit_to_signatures_strict(pt2322_mut_mat, sub.sig, max_delta = 0.1)
pt2322_extra2strictfit_res <- pt2322_extra2strictfit$fit_res
pdf("./3_Output/Mutpatterns/pt2322_SubSetSignatures.pdf")
plot_contribution(pt2322_extra2strictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




pt2322_fit  <- fit_to_signatures(pt2322_mut_mat, sub.sig)
pdf("./3_Output/Mutpatterns/pt2322_SubSetSignatures_Default.pdf")
plot_contribution(pt2322_fit$contribution, coord_flip = FALSE, mode = "absolute") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()






library("ccfindR")
sc <- scNMFSet(count = pt2322_mut_mat)
set.seed(4)
estimate_bayes <- vb_factorize(sc, ranks = 1:8, nrun = 50, 
                               progress.bar = FALSE, verbose = 0)
pdf("./3_Output/Mutpatterns/DeNovo/estimates_Bayes.pdf")
plot(estimate_bayes)
dev.off()






### Refit signatures 
pt2322_fit  <- fit_to_signatures(pt2322_mut_mat, signatures)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_default_refit.pdf")
plot_contribution(pt2322_fit$contribution, coord_flip = FALSE, mode = "absolute") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### Change strictness of refitting
pt2322_strictfit  <- fit_to_signatures_strict(pt2322_mut_mat, signatures, max_delta = 0.004)
pt2322_strictfit_res <- pt2322_strictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_delta0.004_refit.pdf")
plot_contribution(pt2322_strictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### Less strict refitting
pt2322_Lessstrictfit  <- fit_to_signatures_strict(pt2322_mut_mat, signatures, max_delta = 0.001)
pt2322_Lessstrictfit_res <- pt2322_Lessstrictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_delta0.001_refit.pdf")
plot_contribution(pt2322_Lessstrictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### Extra strict refitting
pt2322_extrastrictfit  <- fit_to_signatures_strict(pt2322_mut_mat, signatures, max_delta = 0.01)
pt2322_extrastrictfit_res <- pt2322_extrastrictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_delta0.01_refit.pdf")
plot_contribution(pt2322_extrastrictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### extra2 strict refitting
pt2322_extra2strictfit  <- fit_to_signatures_strict(pt2322_mut_mat, signatures, max_delta = 0.1)
pt2322_extra2strictfit_res <- pt2322_extra2strictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/pt2322_delta0.1_refit.pdf")
plot_contribution(pt2322_extra2strictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()












################################################################# now for pt2229
### Read in data and create mutational matrix 
pt2229_sample_names <- c("ALLBULK", "TALL1", "TALL10", "TALL2", "TALL3", "TALL4", "TALL5", "TALL6", "TALL7", "TALL8", "TALL9")
pt2229_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2229/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE
)
pt2229_grl <- read_vcfs_as_granges(pt2229_vcf_files, pt2229_sample_names, ref_genome)
pt2229_mut_mat <- mut_matrix(vcf_list = pt2229_grl, ref_genome = ref_genome)


### Print the tri-nucleotide plot per sample 
for (s in pt2229_sample_names){
  print(s)
  plotmat <- as.matrix(pt2229_mut_mat[, s])
  colnames(plotmat) <- s
  pdf(paste0("./3_Output/Mutpatterns/96trinuc_profile/pt2229/", s, "_trinuc.pdf"))
  print(plot_96_profile(plotmat))
  dev.off()
}


### Refit signatures 
pt2229_fit  <- fit_to_signatures(pt2229_mut_mat, signatures)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_default_refit.pdf")
plot_contribution(pt2229_fit$contribution, coord_flip = FALSE, mode = "absolute") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### Change strictness of refitting
pt2229_strictfit  <- fit_to_signatures_strict(pt2229_mut_mat, signatures, max_delta = 0.004)
pt2229_strictfit_res <- pt2229_strictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_delta0.004_refit.pdf")
plot_contribution(pt2229_strictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


### Less strict refitting
pt2229_Lessstrictfit  <- fit_to_signatures_strict(pt2229_mut_mat, signatures, max_delta = 0.001)
pt2229_Lessstrictfit_res <- pt2229_Lessstrictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_delta0.001_refit.pdf")
plot_contribution(pt2229_Lessstrictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### Extra strict refitting
pt2229_extrastrictfit  <- fit_to_signatures_strict(pt2229_mut_mat, signatures, max_delta = 0.01)
pt2229_extrastrictfit_res <- pt2229_extrastrictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_delta0.01_refit.pdf")
plot_contribution(pt2229_extrastrictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### extra2 strict refitting
pt2229_extra2strictfit  <- fit_to_signatures_strict(pt2229_mut_mat, signatures, max_delta = 0.1)
pt2229_extra2strictfit_res <- pt2229_extra2strictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2229/pt2229_delta0.1_refit.pdf")
plot_contribution(pt2229_extra2strictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()








# De vcf van de verschillende branches staan hier:
#/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/branches_vcfs
#Dit is het script hoe ze gemaakt zijn:
#/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/TreeBuilding/pt2322/get_branches_vcfs.sh


### Read in data and create mutational matrix 
pt2322branches_sample_names <- c("ISPCD4_branch", "ISPCD4_TALL8_branch", "root_branch", "SPCD4_branch", "TALL7_TALL12_branch","TALL9_TALL11_branch")
pt2322branches_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/Branches_MutPat/pt2322",
                                       pattern = "*.vcf$", full.names = TRUE, recursive = TRUE
)
pt2322branches_grl <- read_vcfs_as_granges(pt2322branches_vcf_files, pt2322branches_sample_names, ref_genome)
pt2322branches_mut_mat <- mut_matrix(vcf_list = pt2322branches_grl, ref_genome = ref_genome)


### Print the tri-nucleotide plot per sample 
for (s in pt2322branches_sample_names){
  print(s)
  plotmat <- as.matrix(pt2322branches_mut_mat[, s])
  colnames(plotmat) <- s
  pdf(paste0("./3_Output/Mutpatterns/96trinuc_profile/pt2322/branches/", s, "_trinuc.pdf"))
  print(plot_96_profile(plotmat))
  dev.off()
}

### Refit signatures 
pt2322branches_fit  <- fit_to_signatures(pt2322branches_mut_mat, signatures)
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/branches/pt2322branches_default_refit.pdf")
plot_contribution(pt2322branches_fit$contribution, coord_flip = FALSE, mode = "absolute") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### Change strictness of refitting
pt2322branches_strictfit  <- fit_to_signatures_strict(pt2322branches_mut_mat, signatures, max_delta = 0.004)
pt2322branches_strictfit_res <- pt2322branches_strictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/branches/pt2322branches_delta0.004_refit.pdf")
plot_contribution(pt2322branches_strictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


### Change strictness of refitting
pt2322branches_ultrastrictfit  <- fit_to_signatures_strict(pt2322branches_mut_mat, signatures, max_delta = 0.1)
pt2322branches_ultrastrictfit_res <- pt2322branches_ultrastrictfit$fit_res
pdf("./3_Output/Mutpatterns/RefittingSignatures/pt2322/branches/pt2322branches_delta0.1_refit.pdf")
plot_contribution(pt2322branches_ultrastrictfit_res$contribution, coord_flip = FALSE, mode = "absolute" ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



cos_sim_pt2322branches <- cos_sim_matrix(pt2322branches_mut_mat, pt2322branches_mut_mat)
plot_cosine_heatmap(cos_sim_pt2322branches, cluster_rows = TRUE, cluster_cols = TRUE)




