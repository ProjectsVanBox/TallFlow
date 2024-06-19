library(tidyverse)
library(ggtree)
library(treeio)
library(stringi)
library(stringdist)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ggnewscale)
ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38'
library(MutationalPatterns)
library(VariantAnnotation)
#theme_set(theme_tree()) --> does not work with 96 plot
theme_set(theme_classic())

library(cellPhyWrapperPlotting) #devtools::install_local("~/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/Scripts/cellPhyWrapperPlotting/",force = TRUE)
source("/Users/v.m.poort/surfdrive - Vera Poort@surfdrive.surf.nl/Shared/pmc_vanboxtel/general/2_Bioinformatics/colors/Jurrians_colors.R")

setwd("/Users/v.m.poort/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/CellPhyWrapper/DxRelsamples/pt344_pt3291/")

tree = readRDS("TreeObject0.3.RDS")
vcf = VariantAnnotation::readVcf("../../../../../1_Input/TreeBuilding/pt344_pt3291/pt344_pt3291_filtered_snv_NoBulk.vcf")

# prepare tree
tree = prepare_tree(tree)

# plot bare tree
plot_gg_tree_base(tree)
plot_gg_tree(tree)
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = TRUE)

# rename the branches
old = tree@phylo$tip.label # get old labels
new = c("pt344-TALL16", "pt344-MONOBULK", "pt344-TALL6", "pt344-TALL7", "pt344-TALL8", "pt344-TALL9", "pt344-TALL5", 
        "pt344-TALL1", "pt344-TALL2", "pt344-TALL3", "pt344-TALL4", "pt344-TALL13", "pt344-TALL10", "pt344-TALL11", 
        "pt344-TALL19", "pt344-TALL18", "pt344-TALL15", "pt344-TALL20", "pt344-TALL17", "pt344-TALL21", "pt344-TALL12", 
        "pt344-TALL14")
tip_names_replace = data.frame(old, new)
tree_rename = change_tip_labels(tree = tree, id_change_tb = tip_names_replace)
plot_gg_tree_base(tree_rename, common_name = 'PB09204 | Male | 7yo | Lymph Node & Bone Marrow')
ggsave("../../pt344_CellphyPloting_branchNames.pdf")
plot_gg_tree(tree_rename, add_branch_length = F, add_bootstrap = F, branch_text_size = 7)
ggsave("../../pt344_CellphyPloting.pdf")
#tree_rename2 = change_tip_labels(tree = tree, id_change_tb = tip_names_replace[1:4, ]) # you can also replace only a few names
#plot_gg_tree_base(tree_rename2)

# VCF per branch
branch_vcf = extract_vcf_per_branch(tree = tree_rename, vcf = vcf, ref_genome = ref_genome)
branch_grl = convert_vcf_to_granges(branch_vcf_list = branch_vcf, ref_genome = ref_genome)

## save branches VCFs
dir.create("../../pt344_branches_vcf")
for ( i in names(branch_vcf)){
  print(i)
  writeVcf(branch_vcf[[i]], filename = paste0("../../pt344_branches_vcf/", i, "_branch.vcf"))
}

# profile per branch --> can be used for any MutationalPatterns analysis separate from the tree
branch_mm = mut_matrix(branch_grl, ref_genome)
#plot_96_profile2(branch_mm, ncol = 5)
#plot_96_profile2(branch_mm[ ,c("P", "B"), drop = FALSE], ncol = 1)

# contributions per branch | plot separate from tree
#cosmic = get_known_signatures()
#pta = readRDS("~/hpc/pmc_vanboxtel/resources/pta_vanboxtel_TP_PTATO_predicted.RDS")[[1]][ ,'False', drop = F] %>% `colnames<-`("PTA") 
#sigs_slct = cosmic[ ,c("SBS1", "SBS5", "SBS18", "SBS29", "SBS40")] # better to first extract signatures in a bigger dataset
#sigs_slct = cbind(sigs_slct, pta)
#contribution = fit_to_signatures_strict_tree(mut_matrix = branch_mm, signatures = sigs_slct, max_delta = 0.01, remove_min = 20)
#plot_contribution(contribution, palette = dist_cols50)

### Get signatures 
signatures = get_known_signatures()
#pta_sig = read.table("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_lymphoma/1_Input/WGS/PTA_Artefact_Signature.txt", sep = "\t", header = T)
#pta_sig = as.matrix(pta_sig)
#PTA <- as.numeric(pta_sig[,"PTA"])
#PTA <- PTA[!is.na(PTA)]
#hspc_sig = read.table("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_lymphoma/1_Input/WGS/sigProfiler_SBS_working_signatures_incl_hspc.txt", sep = "\t", header = T)
hspc_sig = read.table("~/hpc/projects/Tsites/1_Input/sigProfiler_SBS_working_signatures_incl_hspc.txt", sep = "\t", header = T)
hspc_sig = as.matrix(hspc_sig)
HSPC <- as.numeric(hspc_sig[,"HSPC"])
signatures <- cbind(HSPC, signatures)

### Refitting part 
sub_sig <- signatures[,c("HSPC", "SBS1", "SBS5", "SBS40", "SBS18", "SBS29")]
contribution <- fit_to_signatures_strict_tree(mut_matrix = branch_mm, signatures = sub_sig, max_delta = 0.01, remove_min = 40)
plot_contribution(contribution, palette = dist_cols50)

# check if all branches are explained well with these signatures: cosine < 0.85 with > 200 mutations would suggest you miss a mutation
p1a <- check_reconstructed_cosine(contribution, branch_mm, sub_sig, tree_rename)
ggsave( "cosine.pdf")

# add signature contributions to your tree
tree_rename = add_contribution(tree_rename, contribution = contribution) # if you already did signature_fitting
tree = add_contribution(tree, signatures = sigs_slct, mut_matrix = branch_mm) # if you did not yet fit signatures

# contribution per branch | plot in tree
#plot_tree_contribution(tree = tree_rename, signature = 'SBS9', pie_size = 1, common_name = 'SBS9 signature contribution')
#plot_tree_contribution(tree = tree_rename, signature = 'SBS9', type = 'color', legend_pos = 'bottom')
#plot_tree_contribution(tree = tree_rename, signature = 'SBS9', type = 'color_dot', legend_pos = 'bottom')
#plot_tree_contribution(tree = tree_rename, signature = 'PTA', pie_size = 1,common_name = 'PTA signature contribution')
#plot_tree_contribution(tree = tree_rename, signature = 'HSPC', pie_size = 1, common_name = 'HSPC signature contribution')


# plot bars with all signatures in the branches of the tree
p2 <- plot_tree_contribution_bars(tree = tree_rename, signatures = sub_sig, mut_matrix = branch_mm, scaling = 0.3) # give signatures and branch_mut_matrix
ggsave( "PB09204_all_contributionTree.pdf")
#plot_tree_contribution_bars(tree = tree_rename, contribution = contribution) # or give contribution if you have already calculated that
#plot_tree_contribution_bars(tree = tree_rename, contribution = contribution, remove_min = 20) # removing branches with less than 50/100 mutations might be good, as those refits are unreliable

# plot drivers per branch
#tree_rename = add_tree_drivers(tree = tree_rename, branch_vcf = branch_vcf, add_mod_high = T, add_cosmic = T)
#plot_gg_tree(tree_rename, branch_text_param = 'drivers', branch_text_size = 2)

# add drivers manually --> e.g. CNVs: check the "letters" in the default 'plot_gg_tree(tree)'
#manual_drivers = setNames(c('FOXO1_miss', 'IG_SHM', 'IG_SHM', 'IG_SHM', 'MYC_IGH_trans', "SMARCA4_miss/FOXO1_miss"),
#                          c('U', 'K', 'H', 'E', 'M', "O"))
#tree = add_info_to_tree(tree, manual_drivers = manual_drivers)
#plot_gg_tree(tree, branch_text_param = 'manual_drivers', branch_text_size = 3)

# find genes
#find_gene(branch_vcf, 'FOXO1')
#vafs = vcf@assays@data$AD[names(vcf) == '13:40666151_C/G', ] %>% sapply(function(x) x[2]/sum(x)) %>% setNames(samples(header(vcf))) %>% sort

# plot only shared branches
#plot_gg_tree(tree, only_shared_branches = TRUE)
#plot_tree_contribution(tree = tree, signature = 'SBS9', type = 'color', only_shared_branches = TRUE)

# color branches
#plot_gg_tree(tree, branch_color_param = 'manual_drivers') # discreet
#plot_gg_tree(tree, branch_color_param = 'branch_length') # continuous

# group branches -------------------------------
# make labels: each number corresponds to the levels, e.g., 2 = 'naive_B'
# labels should be ordered by the letters in the "plot_gg_tree"
#branch_cat = factor(c(9,9,2,3,2,  2,3,1,3,3,  1,3,4,2,5,  8,8,8,10,8, 7,7,6,9,9,  10,9,9,9,10, 10,10,10,9,10, 10), 
#                    labels = c('GC-B', 'naive_B', 'other', 'MYC_branch', 'tumor_clonal', 
#                               'tumor_subclonal_LN', 'tumor_subclonal_BM', 'tumor_BM_private', 'tumor_LN_private', 'other'))
#branch_cat = setNames(branch_cat, pull(tree@data, branch_id))
#tree = add_info_to_tree(tree, category = branch_cat)
#plot_gg_tree(tree, branch_color_param = 'category')

# summarize mutational_matrix per category and fit signatures
#grouped_branch_mm = group_tree_mut_matrix(tree, mut_matrix = branch_mm, group = 'category')
#plot_96_profile2(grouped_branch_mm)

#grouped_contribution = fit_to_signatures_strict_tree(mut_matrix = grouped_branch_mm, signatures = sigs_slct, max_delta = 0.01, remove_min = 20)
#plot_contribution(grouped_contribution, palette = dist_cols50) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# alternatively, directly summarize contributions per category
#grouped_branch_contri = group_tree_contributions(tree, contribution = contribution, group = 'category')
#plot_contribution(grouped_branch_contri, palette = dist_cols50) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#plot_tree_contribution(tree, branch_color_param = 'manual_drivers', signature = 'SBS9', legend_pos = 'bottom') # discreet
#ggsave("~/Documents/tree_example.png", width = 9, height = 6)



group_tree_contributions <- function(tree, group, contribution = NULL, signatures = NULL, mut_matrix = NULL, ...) 
{
  contri = as.data.frame(contribution)
  branches_missing = setdiff(tree@data$branch_id, colnames(contribution))
  contri = cbind(contri, as.data.frame(matrix(0, ncol = length(branches_missing), nrow = nrow(contri),
                                              dimnames = list(rownames(contri), branches_missing))))
  group = pull(tree@data, group)
  if (!any(class(group) == 'factor')) { group = factor(group) }
  grouped_contri = lapply(levels(group), function(gr) { rowSums(contri[ ,tree@data$branch_id[group == gr],drop=FALSE]) }) %>%
    do.call(cbind, .) %>% `colnames<-`(levels(group))
}
return(grouped_contri)
}

plot_tree_contribution_bars <- function(tree, signatures = NULL, contribution = NULL, 
                                        mut_matrix = NULL, signature_colors = NULL, scaling = 1, remove_min = 20, ...) {
  # colors
  if (is.null(signature_colors)) {
    n_sig = nrow(contribution)
    if (n_sig <= 22) { cols = dist_cols } else 
      if (n_sig <= 31) { cols = dist_cols31 } else 
        if (n_sig <= 50) { cols = dist_cols31 } else 
        { cols = scales::hue_pal()(n_sig) }
    signature_colors = setNames(cols, rownames(contribution))
    signature_colors = signature_colors[!is.na(names(signature_colors))]
  }
  if (is.null(names(signature_colors))) {
    names(signature_colors) = rownames(contribution)
    signature_colors = signature_colors[!is.na(names(signature_colors))]
  }
  # calculate the fractions
  fraction = contribution %>%
    t() %>% as.data.frame %>%
    rownames_to_column('node') %>%
    pivot_longer(cols = -starts_with('node'), names_to = 'sig', values_to = 'contribution')
  fraction$node = tree@data$node[match(fraction$node, tree@data$branch_id)] # use node_id, this is what the "inset" function need later
  # this part of the code has been adapted from ggtree::nodepie
  # all credits go to Guangchuang Yu:https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html
  fractions_split <- split(fraction, fraction$node)
  bars = lapply(fractions_split, function(df) {
    ggplot(df, aes(x = '', y = contribution, fill = sig)) + 
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = signature_colors) +
      theme_void() +
      theme(legend.position = 'none')
  })
  nodes_slct = fractions_split %>% sapply(function(x) pull(x, node)[1]) %>% unname
  # determine the bar lengths
  bar_lengths = tree@data$branch_length/max(tree@data$branch_length)
  bar_lengths = bar_lengths[match(nodes_slct, tree@data$node)] * scaling
  plot = plot_gg_tree(tree, ...)
  tree_plot = ggtree::inset(tree_view = plot, insets = bars, width = bar_lengths, height = 0.1, x = 'branch')
  # add legend
  fig_for_leg = ggplot(fraction, aes(x = '', y = contribution, fill = sig)) +
    geom_col() +
    scale_fill_manual(values = signature_colors)
  legend = cowplot::get_legend(fig_for_leg)
  cowplot::plot_grid(tree_plot, legend, nrow = 1, rel_widths = c(1, 0.2))
