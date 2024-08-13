#!/bin/bash

## Protein DE analysis

## libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(tidyverse)
#BiocManager::install("dittoSeq")
#library(dittoSeq)
library(RColorBrewer)
library(viridis)
library(cowplot)

setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/')
#setwd('/Users/verapoort/hpc/projects/TallClonal/3_Output/Test_objects/')

## load data
srat_pts <- readRDS('srat_pts.rds')
#srat_pts <- readRDS('pts_small.rds')
#prot2gene <- read.csv('../../1_Input/ProteinGeneSymbols.csv', header = TRUE, sep = ",")

#exclude non T-cell types
exclCelltypes <- c("NK", "Monocyte",  "Dendritic",  "Bcell")

# colors for plotting
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A", "#E28028", "#963C25", "grey")
names(pal) <- c("ALLBULK", "DN-like", "DNearly","SPCD4", "DP", "MSCBULK", "MONOBulk","DN3", "iSPCD4", "SPCD8", "gdTcell", "unknown")

# flow markers 
flowmarkers9 <- c("CD7.1","CD5.1", "CD44.1","CD1a", "CD4.1", "CD8","CD3","TCR-A-B","TCR-Y-O")
flowmarkersALL <- c("CD117", "CD34.1", "CD99.1", "CD7.1","CD5.1", "CD2.1","CD44.1","CD127", "CD38.1", "CD10", "CD1a", "CD4.1", "CD8","CD3","TCR-A-B","TCR-Y-O","CD71")
## corresponding genes
flowgenes9 <- c("CD7", "CD5", "CD44", "CD1A", "CD4", "CD8A", "CD3E", 'PTCRA')
flowgenesAll <- c("KIT", "CD34", "CD99","CD7", "CD5", "CD2", "CD44","IL7R", "CD38", "MME", "CD1A", "CD4", "CD8A", "CD3E", "PTCRA", "TFRC")

#ordering
my_levels <- c("DNearly", "DN3", "iSPCD4", "DP", "SPCD4", "SPCD8", "gdTcell","unknown")

#plot_list1 <- list()
#plot_list2 <- list()

## make heatmap
n = 5
for (patient in unique(srat_pts$orig.ident)) {
  #subset per patient
  Idents(srat_pts) <- "orig.ident"
  name = as.character(patient)
  srat_sub <- subset(srat_pts, idents = name)
  #remove the subsets with < ncells for better comparison
  Idents(srat_sub) <- "gating"
  idents.count <- as.data.frame(table(srat_sub@meta.data$gating))
  idents.keep  <- idents.count %>% group_by(Var1) %>% filter(Freq >= n) %>% filter(!Var1 %in% exclCelltypes)
  srat_sub <- subset(srat_sub, idents = idents.keep$Var1)
  # scale all proteins so that all proteins are in the scale.data slot
  DefaultAssay(srat_sub) <- "Protein"
  all.proteins <- rownames(srat_sub) 
  srat_sub <- ScaleData(srat_sub, features = all.proteins)
  # plot heatmap DE proteins
  #markers <- markers[order(match(markers,prot2gene$Protein))] #order markers
  srat_sub$gating <- factor(x = srat_sub$gating, levels = my_levels)
  p1 <- DoHeatmap(srat_sub, features = flowmarkers9, 
                  assay = "Protein", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal )+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "PiYG")) ) +
    labs(title = as.character(name))
  
  ## Plot expression of genes to compare
  DefaultAssay(srat_sub) <- "RNA" # Used RNA because some genes are filtered out in SCT
  srat_sub <- NormalizeData(srat_sub)
  all.genes <- rownames(srat_sub)
  srat_sub <- ScaleData(srat_sub, features = all.genes, do.center = TRUE)
  genes_filt <- intersect(flowgenes9, rownames(srat_sub@assays$RNA))
  
  # plot heatmap genes transcribing DE proteins
  #genes_filt <- genes_filt[order(match(genes_filt, prot2gene$Gene))] #order genes same as markers
  srat_sub$gating <- factor(x = srat_sub$gating, levels = my_levels)
  p2 <- DoHeatmap(srat_sub, features = genes_filt, 
                  assay = "RNA", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")) )+
    labs(title = as.character(name))
  
  p_grid1 <- plot_grid(p1 , p2 , labels = c('Protein', 'RNA'))
  #temp_list1 <- list(grid1=p_grid1)
  ## plot second option with additional markers important for T-development
  DefaultAssay(srat_sub) <- "Protein"
  srat_sub$gating <- factor(x = srat_sub$gating, levels = my_levels)
  p3 <- DoHeatmap(srat_sub, features = flowmarkersALL, 
                  assay = "Protein", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal)+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "PiYG")) ) + 
    labs(title = as.character(name))

  
  ## Plot expression of genes to compare
  DefaultAssay(srat_sub) <- "RNA" # Used RNA because some genes are filtered out in SCT
  genes_filt <- intersect(flowgenesAll, rownames(srat_sub@assays$RNA))
  
  # plot heatmap genes transcribing DE proteins
  #genes_filt <- genes_filt[order(match(genes_filt, prot2gene$Gene))] #order genes same as markers
  srat_sub$gating <- factor(x = srat_sub$gating, levels = my_levels)
  p4 <- DoHeatmap(srat_sub, features = genes_filt, 
                  assay = "RNA", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")) ) +
    labs(title = as.character(name))
  
  p_grid2 <- plot_grid(p3, p4, labels = c('Protein', 'RNA'))
  #temp_list2 <- list(grid2=p_grid2)
  
  pdf(paste0('PerPatient/Pretty/', patient, "_ProtRNA.pdf"), width = 16, height = 8)
  print(p_grid1)
  print(p_grid2)
  #ggsave(p_grid1,  filename = paste0(patient, "_ProtRNA.pdf"), device = "pdf", width = 16, height = 8)
  dev.off()
  saveRDS(p_grid1, paste0('PerPatient/Pretty/', patient, "_ProtRNAplot.rds"))
}

#saveRDS(plot_list1[1:4], 'PerPatient/Pretty/ProtRNA_plotlist1.rds')
#saveRDS(plot_list1[5:8], 'PerPatient/Pretty/ProtRNA_plotlist2.rds')
#saveRDS(plot_list1[9:length(plot_list1)], 'PerPatient/Pretty/ProtRNA_plotlist3.rds')

#plot_list1[1]

#saveRDS(plot_list2, 'PerPatient/Pretty/ProtRNA_ExtMarkers_plotlist.rds')

