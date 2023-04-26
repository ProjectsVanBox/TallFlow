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


## Plot mutation spectra as stacked bar
plot_spectrum_stacked <- function(type_occurrences, mode = "absolute", max_value_to_show = 30, facet = "", total_value = TRUE, round_number = 0, title = ""){
  type_occurrences_stack <- type_occurrences[,c(1,2, 7, 8,4,5,6)]
  ylabel <- "Absolute contribution"
  y_extension <- 0.1
  if(mode == "relative"){
    type_occurrences_stack <- type_occurrences_stack / rowSums(type_occurrences_stack)
    round_number <- 2
    ylabel <- "Relative contribution"
    y_extension <- 0
  }
  
  type_occurrences_stack$Sample <- row.names(type_occurrences_stack)
  type_occurrences_stack$Sample <- factor(type_occurrences_stack$Sample, levels = row.names(type_occurrences_stack))
  if(length(facet > 1)){
    type_occurrences_stack$Facet <- facet
    type_occurrences_stack_m <- melt(type_occurrences_stack, id.vars = c("Sample", "Facet"))
    
  } else {
    type_occurrences_stack_m <- melt(type_occurrences_stack, id.vars = "Sample")
    
  }
  type_occurrences_stack_m$Label <- type_occurrences_stack_m$value
  # The values below max_value_to_show are not plotted
  type_occurrences_stack_m$Label[type_occurrences_stack_m$Label < max_value_to_show] <- ""
  # type_occurrences_total <- aggregate(type_occurrences_stack_m$value ~ type_occurrences_stack_m$Sample, FUN = "sum")
  # colnames(type_occurrences_total) <- c("Sample", "Count")
  if(length(facet) < 2){
    
    plot <- ggplot(data = type_occurrences_stack_m, aes(x = Sample, 
                                                        y = value, 
                                                        label = as.character(round(as.numeric(Label),round_number )))) + 
      geom_bar(stat = "identity", aes(fill = variable), width = 0.8) +
      geom_text(aes(fill = variable), size = 2, position = position_stack(vjust = 0.5)) +
      #facet_grid(.~Type, space = "free",scales = "free", switch = "both") +
      scale_fill_manual(values = COLORS7) +
      scale_y_continuous(expand = expansion(mult = c(0,y_extension))) +
      theme_classic() +
      ggtitle(title) + 
      labs(x = "Sample", y = ylabel, fill = "Point mutation type")  +
      theme(axis.text.x= element_text(angle = 45, hjust = 1), 
            strip.placement = "outside", 
            strip.background = element_rect(fill = "lightgray"), 
            plot.title = element_text(hjust = 0.5))
    if(total_value == TRUE){
      plot +
        geom_text(aes(label = round(stat(y), round_number), group = Sample), 
                  stat = 'summary', fun = sum, vjust = -1, size = 3) 
    } else {
      plot
    }
  } else {
    ggplot(data = type_occurrences_stack_m, aes(x = Sample, y = value, label = as.character(round(as.numeric(Label),round_number )))) + 
      geom_bar(stat = "identity", aes(fill = variable), width = 0.8) +
      geom_text(aes(fill = variable), size = 2, position = position_stack(vjust = 0.5)) +
      geom_text(aes(label = round(stat(y), 0), group = Sample), 
                stat = 'summary', fun = sum, vjust = -1, size = 3) + 
      facet_grid(.~Facet, scales = "free", space = "free") +
      #facet_grid(.~Type, space = "free",scales = "free", switch = "both") +
      scale_fill_manual(values = COLORS7) +
      scale_y_continuous(expand = expansion(mult = c(0,y_extension))) +
      theme_classic() +
      labs(x = "Sample", y = ylabel, fill = "Point mutation type") +
      ggtitle(title) + 
      theme(axis.text.x= element_text(angle = 45, hjust = 1), 
            strip.placement = "outside", 
            strip.background = element_rect(fill = "lightgray"), 
            plot.title = element_text(hjust = 0.5))
  }
}


### Get signatures 
signatures = get_known_signatures()
pta_sig = read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/signatures/PTA_Artefact_Signature.txt", sep = "\t", header = T)
pta_sig = as.matrix(pta_sig)
PTA <- as.numeric(pta_sig[,"PTA"])
hspc_sig = read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/signatures/sigProfiler_SBS_working_signatures_incl_hspc.txt", sep = "\t", header = T)
hspc_sig = as.matrix(hspc_sig)
HSPC <- as.numeric(hspc_sig[,"HSPC"])
signatures <- cbind(PTA, HSPC, signatures)

COLORS6 <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
)
COLORS7 <- c(
  "#2EBAED", "#000000", 
  "#E98C7B", "#DE1C14","#D4D2D2", "#ADCC54",
  "#F0D0CE"
)



########################################################### DeNovo analysis ###########################################################
### Read in data 
vcf_files <- list.files(path = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/snvs/pt2283D/",
                        pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
sample_names <- c("pt2283D-DX1BM-ALLBULK", "pt2283D-DX1BM-TALL1", "pt2283D-DX1BM-TALL10", "pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12", 
                  "pt2283D-DX1BM-TALL2", "pt2283D-DX1BM-TALL3", "pt2283D-DX1BM-TALL4", "pt2283D-DX1BM-TALL5", "pt2283D-DX1BM-TALL6", 
                  "pt2283D-DX1BM-TALL7", "pt2283D-DX1BM-TALL8", "pt2283D-DX1BM-TALL9")
mut_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat_own <- mut_matrix(vcf_list = mut_grl, ref_genome = ref_genome)
nmf_ref_mat = read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/signatures/NMF_reference_matrix.txt", sep = "\t", header = T)
mut_mat <- cbind(mut_mat_own, nmf_ref_mat)

mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", 
                nrun = 50, seed = 123456, .opt = "v-p")
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/DeNovo/pt2283D_Estimates.pdf")
plot(estimate)
dev.off()

### Extract signatures 
nmf_res <- extract_signatures(mut_mat, rank = 4, nrun = 50, single_core = TRUE)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")

### Check with known signatures 
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
colnames(nmf_res$signatures)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/DeNovo/pt2283D_NMF_96profiles.pdf")
plot_96_profile(nmf_res$signatures, condensed = TRUE)
dev.off()


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/DeNovo/pt2283D_NMF_RelContribution.pdf")
plot_contribution(nmf_res$contribution, nmf_res$signature,
                  mode = "relative"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/DeNovo/pt2283D_NMF_OrigRecon.pdf")
plot_original_vs_reconstructed(mut_mat, nmf_res$reconstructed, 
                               y_intercept = 0.95)
dev.off()

nmf_res$contribution[,1:13]
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/DeNovo/pt2283D_NMF_OrigRecon_ourSample.pdf")
plot_contribution(nmf_res$contribution[,1:13], nmf_res$signature,
                  mode = "relative"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

########################################################### DeNovo analysis ###########################################################



plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





########################################################### Standard fit ###########################################################
### Read in data 
vcf_files <- list.files(path = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/snvs/pt2283D/",
                        pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
sample_names <- c("pt2283D-DX1BM-ALLBULK", "pt2283D-DX1BM-TALL1", "pt2283D-DX1BM-TALL10", "pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12", 
                  "pt2283D-DX1BM-TALL2", "pt2283D-DX1BM-TALL3", "pt2283D-DX1BM-TALL4", "pt2283D-DX1BM-TALL5", "pt2283D-DX1BM-TALL6", 
                  "pt2283D-DX1BM-TALL7", "pt2283D-DX1BM-TALL8", "pt2283D-DX1BM-TALL9")
mut_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat_own <- mut_matrix(vcf_list = mut_grl, ref_genome = ref_genome)







type_occurrences <- mut_type_occurrences(mut_grl, ref_genome)
type_occurrences
plot_spectrum(type_occurrences, CT = TRUE)
plot_spectrum(type_occurrences, indv_points = TRUE, CT = TRUE, legend = TRUE)
#plot_spectrum(type_occurrences, by = rownames(type_occurrences), legend = TRUE)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/MutSpectrumStacked.pdf")
plot_spectrum_stacked(type_occurrences)
dev.off()
write.table(type_occurrences, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/MutSpectrum.txt", sep = "\t", col.names = NA)

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/96_profile_Merged.pdf")
plot_96_profile(mut_mat_own)
dev.off()

for (i in colnames(mut_mat_own)){
  print(i)
  sub.mat <- as.matrix(mut_mat_own[,i])
  colnames(sub.mat) <- i
  pdf(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/96_profiles_SingleSample/", i, "_96_profiles.pdf"))
  print(plot_96_profile(sub.mat))
  dev.off()
}


fit_res <- fit_to_signatures(mut_mat_own, signatures)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/ReFit_Cosmic/pt2283D_DefaultRefit.pdf")
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



strict_refit <- fit_to_signatures_strict(mut_mat_own, signatures, max_delta = 0.004)
fig_list <- strict_refit$sim_decay_fig
fig_list[[1]]

strict_refit <- fit_to_signatures_strict(mut_mat_own, signatures, max_delta = 0.1)
fit_res_strict <- strict_refit$fit_res
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/ReFit_Cosmic/pt2283D_StrictRefit_abs.pdf")
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/ReFit_Cosmic/pt2283D_StrictRefit_rel.pdf")
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



contri_boots <- fit_to_signatures_bootstrapped(mut_mat_own,
                                               signatures,
                                               n_boots = 50,
                                               method = "strict",
                                               max_delta = 0.1
)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/ReFit_Cosmic/pt2283D_StrictRefit_bootstrap.pdf")
plot_bootstrapped_contribution(contri_boots)
dev.off()

########################################################### Standard fit ###########################################################




########################################################### subselect ###########################################################
sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18")]
fit_res <- fit_to_signatures(mut_mat_own, sub.sig)

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/SubSelect/pt2283D_abs_contribution.pdf")
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/SubSelect/pt2283D_rel_contribution.pdf")
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

########################################################### subselect ###########################################################





########################################################### Branches ###########################################################
### Read in data 
vcf_files <- list.files(path = "/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283D_new/manual_filter/branches_vcfs/",
                        pattern = "*.vcf$", full.names = TRUE, recursive = FALSE )
sample_names <- c("root_branch", "TALL1_branch", "TALL1_TALL2_TALL3_TALL4_TALL5_TALL6_TALL7_TALL8_TALL9_TALL10_TALL11_branch", 
                  "TALL1_TALL2_TALL4_TALL5_TALL6_TALL7_TALL8_TALL9_TALL10_TALL11_branch", "TALL10_branch", "TALL11_branch", 
                  "TALL12_branch", "TALL2_branch", "TALL3_branch", "TALL4_branch", "TALL5_branch", "TALL6_branch", "TALL7_branch", 
                  "TALL8_branch", "TALL9_branch")
mut_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat_own <- mut_matrix(vcf_list = mut_grl, ref_genome = ref_genome)


sub.sig <- signatures[, c("PTA", "HSPC", "SBS1", "SBS5", "SBS18")]
fit_res <- fit_to_signatures(mut_mat_own, sub.sig)

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/Branched/pt2283D_Branched_abs_contribution.pdf")
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/pt2283D/Branched/pt2283D_Branched_rel_contribution.pdf")
plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + scale_fill_manual(values=mycols_paired) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#/Users/ricohagelaar/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2283D_new/manual_filter/branches_vcfs/



