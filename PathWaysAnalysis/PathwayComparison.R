#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(dplyr)
library(fgsea)
set.seed(54321)

# Global variables 
GO_file <- "/Users/ricohagelaar/Downloads/c5.all.v2023.2.Hs.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)


ThyDEgenes <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/DEgenes_thymi.csv", sep = ",", header = T, row.names = 1)
pt10138 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt10138.csv", sep = ",", header = T, row.names = 1)
pt1179 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt1179.csv", sep = ",", header = T, row.names = 1)
pt11801 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt11801.csv", sep = ",", header = T, row.names = 1)
pt2283 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt2283.csv", sep = ",", header = T, row.names = 1)
pt2337 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt2337.csv", sep = ",", header = T, row.names = 1)
pt3045 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt3045.csv", sep = ",", header = T, row.names = 1)
pt315 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt315.csv", sep = ",", header = T, row.names = 1)
pt3291 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt3291.csv", sep = ",", header = T, row.names = 1)
pt335 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt335.csv", sep = ",", header = T, row.names = 1)
pt344 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt344.csv", sep = ",", header = T, row.names = 1)
pt4068 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt4068.csv", sep = ",", header = T, row.names = 1)
pt5242 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt5242.csv", sep = ",", header = T, row.names = 1)
pt5438 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt5438.csv", sep = ",", header = T, row.names = 1)
pt5676 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt5676.csv", sep = ",", header = T, row.names = 1)
pt9160 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt9160.csv", sep = ",", header = T, row.names = 1)
pt9175 <- read.table("~/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/PerPatient/DEtestingNew/DEgenesWilC_pt9175.csv", sep = ",", header = T, row.names = 1)




GSEA = function(gene_list, GO_file, pval) {
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        maxSize=400, ## maximum gene set size
                        nperm=10000, 
                        scoreType = "pos") %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  #message("Collapsing Pathways -----")
  #concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
  #                                    pathways = myGO,
  #                                    stats = gene_list)
  #fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  #message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) + 
    theme_bw()
  #g1
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}








# Thymus DE genes 
DEgenes_selected <- ThyDEgenes[ThyDEgenes$cluster == "SPCD8", ]
Thy_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- ThyDEgenes[ThyDEgenes$cluster == "DNearly", ]
Thy_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- ThyDEgenes[ThyDEgenes$cluster == "SPCD4", ]
Thy_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- ThyDEgenes[ThyDEgenes$cluster == "DN3", ]
Thy_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- ThyDEgenes[ThyDEgenes$cluster == "DP", ]
Thy_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- ThyDEgenes[ThyDEgenes$cluster == "NK", ]
Thy_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- ThyDEgenes[ThyDEgenes$cluster == "gdTcell", ]
Thy_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt10138 DE genes 
DEgenes_selected <- pt10138[pt10138$cluster == "SPCD8", ]
pt10138_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt10138[pt10138$cluster == "DNearly", ]
pt10138_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt10138[pt10138$cluster == "SPCD4", ]
pt10138_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt10138[pt10138$cluster == "DN3", ]
pt10138_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt10138[pt10138$cluster == "DP", ]
pt10138_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt10138[pt10138$cluster == "NK", ]
pt10138_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt10138[pt10138$cluster == "gdTcell", ]
pt10138_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt1179 DE genes 
DEgenes_selected <- pt1179[pt1179$cluster == "SPCD8", ]
pt1179_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt1179[pt1179$cluster == "DNearly", ]
pt1179_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt1179[pt1179$cluster == "SPCD4", ]
pt1179_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt1179[pt1179$cluster == "DN3", ]
pt1179_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt1179[pt1179$cluster == "DP", ]
pt1179_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt1179[pt1179$cluster == "NK", ]
pt1179_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt1179[pt1179$cluster == "gdTcell", ]
pt1179_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt11801 DE genes 
DEgenes_selected <- pt11801[pt11801$cluster == "SPCD8", ]
pt11801_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt11801[pt11801$cluster == "DNearly", ]
pt11801_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt11801[pt11801$cluster == "SPCD4", ]
pt11801_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt11801[pt11801$cluster == "DN3", ]
pt11801_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt11801[pt11801$cluster == "DP", ]
pt11801_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt11801[pt11801$cluster == "NK", ]
pt11801_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt11801[pt11801$cluster == "gdTcell", ]
pt11801_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt2283 DE genes 
DEgenes_selected <- pt2283[pt2283$cluster == "SPCD8", ]
pt2283_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2283[pt2283$cluster == "DNearly", ]
pt2283_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2283[pt2283$cluster == "SPCD4", ]
pt2283_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2283[pt2283$cluster == "DN3", ]
pt2283_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2283[pt2283$cluster == "DP", ]
pt2283_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2283[pt2283$cluster == "NK", ]
pt2283_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2283[pt2283$cluster == "gdTcell", ]
pt2283_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt2337 DE genes 
DEgenes_selected <- pt2337[pt2337$cluster == "SPCD8", ]
pt2337_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2337[pt2337$cluster == "DNearly", ]
pt2337_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2337[pt2337$cluster == "SPCD4", ]
pt2337_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2337[pt2337$cluster == "DN3", ]
pt2337_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2337[pt2337$cluster == "DP", ]
pt2337_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2337[pt2337$cluster == "NK", ]
pt2337_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt2337[pt2337$cluster == "gdTcell", ]
pt2337_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt3045 DE genes 
DEgenes_selected <- pt3045[pt3045$cluster == "SPCD8", ]
pt3045_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3045[pt3045$cluster == "DNearly", ]
pt3045_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3045[pt3045$cluster == "SPCD4", ]
pt3045_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3045[pt3045$cluster == "DN3", ]
pt3045_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3045[pt3045$cluster == "DP", ]
pt3045_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3045[pt3045$cluster == "NK", ]
pt3045_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3045[pt3045$cluster == "gdTcell", ]
pt3045_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt315 DE genes 
DEgenes_selected <- pt315[pt315$cluster == "SPCD8", ]
pt315_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt315[pt315$cluster == "DNearly", ]
pt315_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt315[pt315$cluster == "SPCD4", ]
pt315_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt315[pt315$cluster == "DN3", ]
pt315_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt315[pt315$cluster == "DP", ]
pt315_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt315[pt315$cluster == "NK", ]
pt315_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt315[pt315$cluster == "gdTcell", ]
pt315_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt3291 DE genes 
DEgenes_selected <- pt3291[pt3291$cluster == "SPCD8", ]
pt3291_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3291[pt3291$cluster == "DNearly", ]
pt3291_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3291[pt3291$cluster == "SPCD4", ]
pt3291_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3291[pt3291$cluster == "DN3", ]
pt3291_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3291[pt3291$cluster == "DP", ]
pt3291_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3291[pt3291$cluster == "NK", ]
pt3291_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt3291[pt3291$cluster == "gdTcell", ]
pt3291_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt335 DE genes 
DEgenes_selected <- pt335[pt335$cluster == "SPCD8", ]
pt335_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt335[pt335$cluster == "DNearly", ]
pt335_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt335[pt335$cluster == "SPCD4", ]
pt335_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt335[pt335$cluster == "DN3", ]
pt335_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt335[pt335$cluster == "DP", ]
pt335_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt335[pt335$cluster == "NK", ]
pt335_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt335[pt335$cluster == "gdTcell", ]
pt335_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt344 DE genes 
DEgenes_selected <- pt344[pt344$cluster == "SPCD8", ]
pt344_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt344[pt344$cluster == "DNearly", ]
pt344_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt344[pt344$cluster == "SPCD4", ]
pt344_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt344[pt344$cluster == "DN3", ]
pt344_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt344[pt344$cluster == "DP", ]
pt344_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt344[pt344$cluster == "NK", ]
pt344_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt344[pt344$cluster == "gdTcell", ]
pt344_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt4068 DE genes 
DEgenes_selected <- pt4068[pt4068$cluster == "SPCD8", ]
pt4068_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt4068[pt4068$cluster == "DNearly", ]
pt4068_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt4068[pt4068$cluster == "SPCD4", ]
pt4068_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt4068[pt4068$cluster == "DN3", ]
pt4068_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt4068[pt4068$cluster == "DP", ]
pt4068_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt4068[pt4068$cluster == "NK", ]
pt4068_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt4068[pt4068$cluster == "gdTcell", ]
pt4068_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt5242 DE genes 
DEgenes_selected <- pt5242[pt5242$cluster == "SPCD8", ]
pt5242_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5242[pt5242$cluster == "DNearly", ]
pt5242_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5242[pt5242$cluster == "SPCD4", ]
pt5242_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5242[pt5242$cluster == "DN3", ]
pt5242_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5242[pt5242$cluster == "DP", ]
pt5242_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5242[pt5242$cluster == "NK", ]
pt5242_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5242[pt5242$cluster == "gdTcell", ]
pt5242_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt5438 DE genes 
DEgenes_selected <- pt5438[pt5438$cluster == "SPCD8", ]
pt5438_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5438[pt5438$cluster == "DNearly", ]
pt5438_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5438[pt5438$cluster == "SPCD4", ]
pt5438_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5438[pt5438$cluster == "DN3", ]
pt5438_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5438[pt5438$cluster == "DP", ]
pt5438_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5438[pt5438$cluster == "NK", ]
pt5438_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5438[pt5438$cluster == "gdTcell", ]
pt5438_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt5676 DE genes 
DEgenes_selected <- pt5676[pt5676$cluster == "SPCD8", ]
pt5676_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5676[pt5676$cluster == "DNearly", ]
pt5676_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5676[pt5676$cluster == "SPCD4", ]
pt5676_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5676[pt5676$cluster == "DN3", ]
pt5676_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5676[pt5676$cluster == "DP", ]
pt5676_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5676[pt5676$cluster == "NK", ]
pt5676_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt5676[pt5676$cluster == "gdTcell", ]
pt5676_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt9160 DE genes 
DEgenes_selected <- pt9160[pt9160$cluster == "SPCD8", ]
pt9160_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9160[pt9160$cluster == "DNearly", ]
pt9160_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9160[pt9160$cluster == "SPCD4", ]
pt9160_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9160[pt9160$cluster == "DN3", ]
pt9160_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9160[pt9160$cluster == "DP", ]
pt9160_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9160[pt9160$cluster == "NK", ]
pt9160_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9160[pt9160$cluster == "gdTcell", ]
pt9160_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)

# pt9175 DE genes 
DEgenes_selected <- pt9175[pt9175$cluster == "SPCD8", ]
pt9175_SPCD8 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9175[pt9175$cluster == "DNearly", ]
pt9175_DNearly <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9175[pt9175$cluster == "SPCD4", ]
pt9175_SPCD4 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9175[pt9175$cluster == "DN3", ]
pt9175_DN3 <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9175[pt9175$cluster == "DP", ]
pt9175_DP <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9175[pt9175$cluster == "NK", ]
pt9175_NK <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)
DEgenes_selected <- pt9175[pt9175$cluster == "gdTcell", ]
pt9175_gdTcell <- GSEA(setNames(c(DEgenes_selected$avg_log2FC), c(DEgenes_selected$gene)), GO_file, 0.3)



# SPCD8
if (nrow(Thy_SPCD8$Results) > 0){
  write.table(x = Thy_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/Thy_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/Thy_SPCD8.pdf", width = 12)
  print(Thy_SPCD8$Plot)
  dev.off()
}
if (nrow(pt10138_SPCD8$Results) > 0){
  write.table(x = pt10138_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt10138_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt10138_SPCD8.pdf", width = 12)
  print(pt10138_SPCD8$Plot)
  dev.off()
}
if (nrow(pt1179_SPCD8$Results) > 0){
  write.table(x = pt1179_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt1179_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt1179_SPCD8.pdf", width = 12)
  print(pt1179_SPCD8$Plot)
  dev.off()
}
if (nrow(pt11801_SPCD8$Results) > 0){
  write.table(x = pt11801_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt11801_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt11801_SPCD8.pdf", width = 12)
  print(pt11801_SPCD8$Plot)
  dev.off()
}
if (nrow(pt2283_SPCD8$Results) > 0){
  write.table(x = pt2283_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt2283_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt2283_SPCD8.pdf", width = 12)
  print(pt2283_SPCD8$Plot)
  dev.off()
}
if (nrow(pt2337_SPCD8$Results) > 0){
  write.table(x = pt2337_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt2337_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt2337_SPCD8.pdf", width = 12)
  print(pt2337_SPCD8$Plot)
  dev.off()
}
if (nrow(pt3045_SPCD8$Results) > 0){
  write.table(x = pt3045_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt3045_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt3045_SPCD8.pdf", width = 12)
  print(pt3045_SPCD8$Plot)
  dev.off()
}
if (nrow(pt315_SPCD8$Results) > 0){
  write.table(x = pt315_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt315_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt315_SPCD8.pdf", width = 12)
  print(pt315_SPCD8$Plot)
  dev.off()
}
if (nrow(pt3291_SPCD8$Results) > 0){
  write.table(x = pt3291_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt3291_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt3291_SPCD8.pdf", width = 12)
  print(pt3291_SPCD8$Plot)
  dev.off()
}
if (nrow(pt335_SPCD8$Results) > 0){
  write.table(x = pt335_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt335_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt335_SPCD8.pdf", width = 12)
  print(pt335_SPCD8$Plot)
  dev.off()
}
if (nrow(pt344_SPCD8$Results) > 0){
  write.table(x = pt344_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt344_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt344_SPCD8.pdf", width = 12)
  print(pt344_SPCD8$Plot)
  dev.off()
}
if (nrow(pt4068_SPCD8$Results) > 0){
  write.table(x = pt4068_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt4068_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt4068_SPCD8.pdf", width = 12)
  print(pt4068_SPCD8$Plot)
  dev.off()
}
if (nrow(pt5242_SPCD8$Results) > 0){
  write.table(x = pt5242_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt5242_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt5242_SPCD8.pdf", width = 12)
  print(pt5242_SPCD8$Plot)
  dev.off()
}
if (nrow(pt5438_SPCD8$Results) > 0){
  write.table(x = pt5438_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt5438_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt5438_SPCD8.pdf", width = 12)
  print(pt5438_SPCD8$Plot)
  dev.off()
}
if (nrow(pt5676_SPCD8$Results) > 0){
  write.table(x = pt5676_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt5676_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt5676_SPCD8.pdf", width = 12)
  print(pt5676_SPCD8$Plot)
  dev.off()
}
if (nrow(pt9160_SPCD8$Results) > 0){
  write.table(x = pt9160_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt9160_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt9160_SPCD8.pdf", width = 12)
  print(pt9160_SPCD8$Plot)
  dev.off()
}
if (nrow(pt9175_SPCD8$Results) > 0){
  write.table(x = pt9175_SPCD8$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt9175_SPCD8.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD8/pt9175_SPCD8.pdf", width = 12)
  print(pt9175_SPCD8$Plot)
  dev.off()
}


# SPCD4
if (nrow(Thy_SPCD4$Results) > 0){
  write.table(x = Thy_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/Thy_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/Thy_SPCD4.pdf", width = 12)
  print(Thy_SPCD4$Plot)
  dev.off()
}
if (nrow(pt10138_SPCD4$Results) > 0){
  write.table(x = pt10138_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt10138_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt10138_SPCD4.pdf", width = 12)
  print(pt10138_SPCD4$Plot)
  dev.off()
}
if (nrow(pt1179_SPCD4$Results) > 0){
  write.table(x = pt1179_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt1179_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt1179_SPCD4.pdf", width = 12)
  print(pt1179_SPCD4$Plot)
  dev.off()
}
if (nrow(pt11801_SPCD4$Results) > 0){
  write.table(x = pt11801_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt11801_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt11801_SPCD4.pdf", width = 12)
  print(pt11801_SPCD4$Plot)
  dev.off()
}
if (nrow(pt2283_SPCD4$Results) > 0){
  write.table(x = pt2283_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt2283_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt2283_SPCD4.pdf", width = 12)
  print(pt2283_SPCD4$Plot)
  dev.off()
}
if (nrow(pt2337_SPCD4$Results) > 0){
  write.table(x = pt2337_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt2337_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt2337_SPCD4.pdf", width = 12)
  print(pt2337_SPCD4$Plot)
  dev.off()
}
if (nrow(pt3045_SPCD4$Results) > 0){
  write.table(x = pt3045_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt3045_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt3045_SPCD4.pdf", width = 12)
  print(pt3045_SPCD4$Plot)
  dev.off()
}
if (nrow(pt315_SPCD4$Results) > 0){
  write.table(x = pt315_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt315_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt315_SPCD4.pdf", width = 12)
  print(pt315_SPCD4$Plot)
  dev.off()
}
if (nrow(pt3291_SPCD4$Results) > 0){
  write.table(x = pt3291_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt3291_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt3291_SPCD4.pdf", width = 12)
  print(pt3291_SPCD4$Plot)
  dev.off()
}
if (nrow(pt335_SPCD4$Results) > 0){
  write.table(x = pt335_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt335_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt335_SPCD4.pdf", width = 12)
  print(pt335_SPCD4$Plot)
  dev.off()
}
if (nrow(pt344_SPCD4$Results) > 0){
  write.table(x = pt344_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt344_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt344_SPCD4.pdf", width = 12)
  print(pt344_SPCD4$Plot)
  dev.off()
}
if (nrow(pt4068_SPCD4$Results) > 0){
  write.table(x = pt4068_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt4068_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt4068_SPCD4.pdf", width = 12)
  print(pt4068_SPCD4$Plot)
  dev.off()
}
if (nrow(pt5242_SPCD4$Results) > 0){
  write.table(x = pt5242_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt5242_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt5242_SPCD4.pdf", width = 12)
  print(pt5242_SPCD4$Plot)
  dev.off()
}
if (nrow(pt5438_SPCD4$Results) > 0){
  write.table(x = pt5438_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt5438_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt5438_SPCD4.pdf", width = 12)
  print(pt5438_SPCD4$Plot)
  dev.off()
}
if (nrow(pt5676_SPCD4$Results) > 0){
  write.table(x = pt5676_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt5676_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt5676_SPCD4.pdf", width = 12)
  print(pt5676_SPCD4$Plot)
  dev.off()
}
if (nrow(pt9160_SPCD4$Results) > 0){
  write.table(x = pt9160_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt9160_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt9160_SPCD4.pdf", width = 12)
  print(pt9160_SPCD4$Plot)
  dev.off()
}
if (nrow(pt9175_SPCD4$Results) > 0){
  write.table(x = pt9175_SPCD4$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt9175_SPCD4.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/SPCD4/pt9175_SPCD4.pdf", width = 12)
  print(pt9175_SPCD4$Plot)
  dev.off()
}



# DNearly
if (nrow(Thy_DNearly$Results) > 0){
  write.table(x = Thy_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/Thy_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/Thy_DNearly.pdf", width = 12)
  print(Thy_DNearly$Plot)
  dev.off()
}
if (nrow(pt10138_DNearly$Results) > 0){
  write.table(x = pt10138_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt10138_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt10138_DNearly.pdf", width = 12)
  print(pt10138_DNearly$Plot)
  dev.off()
}
if (nrow(pt1179_DNearly$Results) > 0){
  write.table(x = pt1179_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt1179_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt1179_DNearly.pdf", width = 12)
  print(pt1179_DNearly$Plot)
  dev.off()
}
if (nrow(pt11801_DNearly$Results) > 0){
  write.table(x = pt11801_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt11801_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt11801_DNearly.pdf", width = 12)
  print(pt11801_DNearly$Plot)
  dev.off()
}
if (nrow(pt2283_DNearly$Results) > 0){
  write.table(x = pt2283_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt2283_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt2283_DNearly.pdf", width = 12)
  print(pt2283_DNearly$Plot)
  dev.off()
}
if (nrow(pt2337_DNearly$Results) > 0){
  write.table(x = pt2337_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt2337_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt2337_DNearly.pdf", width = 12)
  print(pt2337_DNearly$Plot)
  dev.off()
}
if (nrow(pt3045_DNearly$Results) > 0){
  write.table(x = pt3045_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt3045_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt3045_DNearly.pdf", width = 12)
  print(pt3045_DNearly$Plot)
  dev.off()
}
if (nrow(pt315_DNearly$Results) > 0){
  write.table(x = pt315_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt315_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt315_DNearly.pdf", width = 12)
  print(pt315_DNearly$Plot)
  dev.off()
}
if (nrow(pt3291_DNearly$Results) > 0){
  write.table(x = pt3291_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt3291_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt3291_DNearly.pdf", width = 12)
  print(pt3291_DNearly$Plot)
  dev.off()
}
if (nrow(pt335_DNearly$Results) > 0){
  write.table(x = pt335_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt335_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt335_DNearly.pdf", width = 12)
  print(pt335_DNearly$Plot)
  dev.off()
}
if (nrow(pt344_DNearly$Results) > 0){
  write.table(x = pt344_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt344_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt344_DNearly.pdf", width = 12)
  print(pt344_DNearly$Plot)
  dev.off()
}
if (nrow(pt4068_DNearly$Results) > 0){
  write.table(x = pt4068_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt4068_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt4068_DNearly.pdf", width = 12)
  print(pt4068_DNearly$Plot)
  dev.off()
}
if (nrow(pt5242_DNearly$Results) > 0){
  write.table(x = pt5242_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt5242_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt5242_DNearly.pdf", width = 12)
  print(pt5242_DNearly$Plot)
  dev.off()
}
if (nrow(pt5438_DNearly$Results) > 0){
  write.table(x = pt5438_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt5438_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt5438_DNearly.pdf", width = 12)
  print(pt5438_DNearly$Plot)
  dev.off()
}
if (nrow(pt5676_DNearly$Results) > 0){
  write.table(x = pt5676_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt5676_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt5676_DNearly.pdf", width = 12)
  print(pt5676_DNearly$Plot)
  dev.off()
}
if (nrow(pt9160_DNearly$Results) > 0){
  write.table(x = pt9160_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt9160_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt9160_DNearly.pdf", width = 12)
  print(pt9160_DNearly$Plot)
  dev.off()
}
if (nrow(pt9175_DNearly$Results) > 0){
  write.table(x = pt9175_DNearly$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt9175_DNearly.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DNearly/pt9175_DNearly.pdf", width = 12)
  print(pt9175_DNearly$Plot)
  dev.off()
}



# DN3
if (nrow(Thy_DN3$Results) > 0){
  write.table(x = Thy_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/Thy_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/Thy_DN3.pdf", width = 12)
  print(Thy_DN3$Plot)
  dev.off()
}
if (nrow(pt10138_DN3$Results) > 0){
  write.table(x = pt10138_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt10138_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt10138_DN3.pdf", width = 12)
  print(pt10138_DN3$Plot)
  dev.off()
}
if (nrow(pt1179_DN3$Results) > 0){
  write.table(x = pt1179_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt1179_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt1179_DN3.pdf", width = 12)
  print(pt1179_DN3$Plot)
  dev.off()
}
if (nrow(pt11801_DN3$Results) > 0){
  write.table(x = pt11801_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt11801_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt11801_DN3.pdf", width = 12)
  print(pt11801_DN3$Plot)
  dev.off()
}
if (nrow(pt2283_DN3$Results) > 0){
  write.table(x = pt2283_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt2283_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt2283_DN3.pdf", width = 12)
  print(pt2283_DN3$Plot)
  dev.off()
}
if (nrow(pt2337_DN3$Results) > 0){
  write.table(x = pt2337_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt2337_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt2337_DN3.pdf", width = 12)
  print(pt2337_DN3$Plot)
  dev.off()
}
if (nrow(pt3045_DN3$Results) > 0){
  write.table(x = pt3045_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt3045_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt3045_DN3.pdf", width = 12)
  print(pt3045_DN3$Plot)
  dev.off()
}
if (nrow(pt315_DN3$Results) > 0){
  write.table(x = pt315_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt315_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt315_DN3.pdf", width = 12)
  print(pt315_DN3$Plot)
  dev.off()
}
if (nrow(pt3291_DN3$Results) > 0){
  write.table(x = pt3291_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt3291_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt3291_DN3.pdf", width = 12)
  print(pt3291_DN3$Plot)
  dev.off()
}
if (nrow(pt335_DN3$Results) > 0){
  write.table(x = pt335_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt335_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt335_DN3.pdf", width = 12)
  print(pt335_DN3$Plot)
  dev.off()
}
if (nrow(pt344_DN3$Results) > 0){
  write.table(x = pt344_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt344_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt344_DN3.pdf", width = 12)
  print(pt344_DN3$Plot)
  dev.off()
}
if (nrow(pt4068_DN3$Results) > 0){
  write.table(x = pt4068_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt4068_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt4068_DN3.pdf", width = 12)
  print(pt4068_DN3$Plot)
  dev.off()
}
if (nrow(pt5242_DN3$Results) > 0){
  write.table(x = pt5242_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt5242_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt5242_DN3.pdf", width = 12)
  print(pt5242_DN3$Plot)
  dev.off()
}
if (nrow(pt5438_DN3$Results) > 0){
  write.table(x = pt5438_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt5438_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt5438_DN3.pdf", width = 12)
  print(pt5438_DN3$Plot)
  dev.off()
}
if (nrow(pt5676_DN3$Results) > 0){
  write.table(x = pt5676_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt5676_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt5676_DN3.pdf", width = 12)
  print(pt5676_DN3$Plot)
  dev.off()
}
if (nrow(pt9160_DN3$Results) > 0){
  write.table(x = pt9160_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt9160_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt9160_DN3.pdf", width = 12)
  print(pt9160_DN3$Plot)
  dev.off()
}
if (nrow(pt9175_DN3$Results) > 0){
  write.table(x = pt9175_DN3$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt9175_DN3.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DN3/pt9175_DN3.pdf", width = 12)
  print(pt9175_DN3$Plot)
  dev.off()
}


# DP
if (nrow(Thy_DP$Results) > 0){
  write.table(x = Thy_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/Thy_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/Thy_DP.pdf", width = 12)
  print(Thy_DP$Plot)
  dev.off()
}
if (nrow(pt10138_DP$Results) > 0){
  write.table(x = pt10138_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt10138_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt10138_DP.pdf", width = 12)
  print(pt10138_DP$Plot)
  dev.off()
}
if (nrow(pt1179_DP$Results) > 0){
  write.table(x = pt1179_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt1179_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt1179_DP.pdf", width = 12)
  print(pt1179_DP$Plot)
  dev.off()
}
if (nrow(pt11801_DP$Results) > 0){
  write.table(x = pt11801_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt11801_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt11801_DP.pdf", width = 12)
  print(pt11801_DP$Plot)
  dev.off()
}
if (nrow(pt2283_DP$Results) > 0){
  write.table(x = pt2283_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt2283_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt2283_DP.pdf", width = 12)
  print(pt2283_DP$Plot)
  dev.off()
}
if (nrow(pt2337_DP$Results) > 0){
  write.table(x = pt2337_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt2337_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt2337_DP.pdf", width = 12)
  print(pt2337_DP$Plot)
  dev.off()
}
if (nrow(pt3045_DP$Results) > 0){
  write.table(x = pt3045_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt3045_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt3045_DP.pdf", width = 12)
  print(pt3045_DP$Plot)
  dev.off()
}
if (nrow(pt315_DP$Results) > 0){
  write.table(x = pt315_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt315_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt315_DP.pdf", width = 12)
  print(pt315_DP$Plot)
  dev.off()
}
if (nrow(pt3291_DP$Results) > 0){
  write.table(x = pt3291_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt3291_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt3291_DP.pdf", width = 12)
  print(pt3291_DP$Plot)
  dev.off()
}
if (nrow(pt335_DP$Results) > 0){
  write.table(x = pt335_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt335_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt335_DP.pdf", width = 12)
  print(pt335_DP$Plot)
  dev.off()
}
if (nrow(pt344_DP$Results) > 0){
  write.table(x = pt344_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt344_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt344_DP.pdf", width = 12)
  print(pt344_DP$Plot)
  dev.off()
}
if (nrow(pt4068_DP$Results) > 0){
  write.table(x = pt4068_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt4068_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt4068_DP.pdf", width = 12)
  print(pt4068_DP$Plot)
  dev.off()
}
if (nrow(pt5242_DP$Results) > 0){
  write.table(x = pt5242_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt5242_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt5242_DP.pdf", width = 12)
  print(pt5242_DP$Plot)
  dev.off()
}
if (nrow(pt5438_DP$Results) > 0){
  write.table(x = pt5438_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt5438_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt5438_DP.pdf", width = 12)
  print(pt5438_DP$Plot)
  dev.off()
}
if (nrow(pt5676_DP$Results) > 0){
  write.table(x = pt5676_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt5676_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt5676_DP.pdf", width = 12)
  print(pt5676_DP$Plot)
  dev.off()
}
if (nrow(pt9160_DP$Results) > 0){
  write.table(x = pt9160_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt9160_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt9160_DP.pdf", width = 12)
  print(pt9160_DP$Plot)
  dev.off()
}
if (nrow(pt9175_DP$Results) > 0){
  write.table(x = pt9175_DP$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt9175_DP.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/DP/pt9175_DP.pdf", width = 12)
  print(pt9175_DP$Plot)
  dev.off()
}


# NK
if (nrow(Thy_NK$Results) > 0){
  write.table(x = Thy_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/Thy_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/Thy_NK.pdf", width = 12)
  print(Thy_NK$Plot)
  dev.off()
}
if (nrow(pt10138_NK$Results) > 0){
  write.table(x = pt10138_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt10138_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt10138_NK.pdf", width = 12)
  print(pt10138_NK$Plot)
  dev.off()
}
if (nrow(pt1179_NK$Results) > 0){
  write.table(x = pt1179_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt1179_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt1179_NK.pdf", width = 12)
  print(pt1179_NK$Plot)
  dev.off()
}
if (nrow(pt11801_NK$Results) > 0){
  write.table(x = pt11801_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt11801_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt11801_NK.pdf", width = 12)
  print(pt11801_NK$Plot)
  dev.off()
}
if (nrow(pt2283_NK$Results) > 0){
  write.table(x = pt2283_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt2283_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt2283_NK.pdf", width = 12)
  print(pt2283_NK$Plot)
  dev.off()
}
if (nrow(pt2337_NK$Results) > 0){
  write.table(x = pt2337_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt2337_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt2337_NK.pdf", width = 12)
  print(pt2337_NK$Plot)
  dev.off()
}
if (nrow(pt3045_NK$Results) > 0){
  write.table(x = pt3045_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt3045_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt3045_NK.pdf", width = 12)
  print(pt3045_NK$Plot)
  dev.off()
}
if (nrow(pt315_NK$Results) > 0){
  write.table(x = pt315_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt315_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt315_NK.pdf", width = 12)
  print(pt315_NK$Plot)
  dev.off()
}
if (nrow(pt3291_NK$Results) > 0){
  write.table(x = pt3291_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt3291_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt3291_NK.pdf", width = 12)
  print(pt3291_NK$Plot)
  dev.off()
}
if (nrow(pt335_NK$Results) > 0){
  write.table(x = pt335_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt335_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt335_NK.pdf", width = 12)
  print(pt335_NK$Plot)
  dev.off()
}
if (nrow(pt344_NK$Results) > 0){
  write.table(x = pt344_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt344_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt344_NK.pdf", width = 12)
  print(pt344_NK$Plot)
  dev.off()
}
if (nrow(pt4068_NK$Results) > 0){
  write.table(x = pt4068_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt4068_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt4068_NK.pdf", width = 12)
  print(pt4068_NK$Plot)
  dev.off()
}
if (nrow(pt5242_NK$Results) > 0){
  write.table(x = pt5242_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt5242_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt5242_NK.pdf", width = 12)
  print(pt5242_NK$Plot)
  dev.off()
}
if (nrow(pt5438_NK$Results) > 0){
  write.table(x = pt5438_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt5438_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt5438_NK.pdf", width = 12)
  print(pt5438_NK$Plot)
  dev.off()
}
if (nrow(pt5676_NK$Results) > 0){
  write.table(x = pt5676_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt5676_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt5676_NK.pdf", width = 12)
  print(pt5676_NK$Plot)
  dev.off()
}
if (nrow(pt9160_NK$Results) > 0){
  write.table(x = pt9160_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt9160_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt9160_NK.pdf", width = 12)
  print(pt9160_NK$Plot)
  dev.off()
}
if (nrow(pt9175_NK$Results) > 0){
  write.table(x = pt9175_NK$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt9175_NK.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/NK/pt9175_NK.pdf", width = 12)
  print(pt9175_NK$Plot)
  dev.off()
}


# gdTcell
if (nrow(Thy_gdTcell$Results) > 0){
  write.table(x = Thy_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/Thy_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/Thy_gdTcell.pdf", width = 12)
  print(Thy_gdTcell$Plot)
  dev.off()
}
if (nrow(pt10138_gdTcell$Results) > 0){
  write.table(x = pt10138_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt10138_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt10138_gdTcell.pdf", width = 12)
  print(pt10138_gdTcell$Plot)
  dev.off()
}
if (nrow(pt1179_gdTcell$Results) > 0){
  write.table(x = pt1179_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt1179_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt1179_gdTcell.pdf", width = 12)
  print(pt1179_gdTcell$Plot)
  dev.off()
}
if (nrow(pt11801_gdTcell$Results) > 0){
  write.table(x = pt11801_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt11801_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt11801_gdTcell.pdf", width = 12)
  print(pt11801_gdTcell$Plot)
  dev.off()
}
if (nrow(pt2283_gdTcell$Results) > 0){
  write.table(x = pt2283_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt2283_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt2283_gdTcell.pdf", width = 12)
  print(pt2283_gdTcell$Plot)
  dev.off()
}
if (nrow(pt2337_gdTcell$Results) > 0){
  write.table(x = pt2337_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt2337_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt2337_gdTcell.pdf", width = 12)
  print(pt2337_gdTcell$Plot)
  dev.off()
}
if (nrow(pt3045_gdTcell$Results) > 0){
  write.table(x = pt3045_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt3045_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt3045_gdTcell.pdf", width = 12)
  print(pt3045_gdTcell$Plot)
  dev.off()
}
if (nrow(pt315_gdTcell$Results) > 0){
  write.table(x = pt315_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt315_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt315_gdTcell.pdf", width = 12)
  print(pt315_gdTcell$Plot)
  dev.off()
}
if (nrow(pt3291_gdTcell$Results) > 0){
  write.table(x = pt3291_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt3291_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt3291_gdTcell.pdf", width = 12)
  print(pt3291_gdTcell$Plot)
  dev.off()
}
if (nrow(pt335_gdTcell$Results) > 0){
  write.table(x = pt335_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt335_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt335_gdTcell.pdf", width = 12)
  print(pt335_gdTcell$Plot)
  dev.off()
}
if (nrow(pt344_gdTcell$Results) > 0){
  write.table(x = pt344_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt344_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt344_gdTcell.pdf", width = 12)
  print(pt344_gdTcell$Plot)
  dev.off()
}
if (nrow(pt4068_gdTcell$Results) > 0){
  write.table(x = pt4068_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt4068_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt4068_gdTcell.pdf", width = 12)
  print(pt4068_gdTcell$Plot)
  dev.off()
}
if (nrow(pt5242_gdTcell$Results) > 0){
  write.table(x = pt5242_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt5242_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt5242_gdTcell.pdf", width = 12)
  print(pt5242_gdTcell$Plot)
  dev.off()
}
if (nrow(pt5438_gdTcell$Results) > 0){
  write.table(x = pt5438_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt5438_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt5438_gdTcell.pdf", width = 12)
  print(pt5438_gdTcell$Plot)
  dev.off()
}
if (nrow(pt5676_gdTcell$Results) > 0){
  write.table(x = pt5676_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt5676_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt5676_gdTcell.pdf", width = 12)
  print(pt5676_gdTcell$Plot)
  dev.off()
}
if (nrow(pt9160_gdTcell$Results) > 0){
  write.table(x = pt9160_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt9160_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt9160_gdTcell.pdf", width = 12)
  print(pt9160_gdTcell$Plot)
  dev.off()
}
if (nrow(pt9175_gdTcell$Results) > 0){
  write.table(x = pt9175_gdTcell$Results[,1:7], file = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt9175_gdTcell.txt", sep = "\t", row.names = FALSE)
  pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PathwaysAnalysis/gdTcell/pt9175_gdTcell.pdf", width = 12)
  print(pt9175_gdTcell$Plot)
  dev.off()
}







