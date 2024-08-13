#!/bin/bash

## Protein DE analysis

## libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(tidyverse)
BiocManager::install("dittoSeq")
library(dittoSeq)
library(RColorBrewer)
library(viridis)
library(cowplot)

setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb')

## load data
srat_pts <- readRDS('srat_pts.rds')
prot2gene <- read.csv('../../1_Input/ProteinGeneSymbols.csv', header = TRUE, sep = ",")

#exclude non T-cell types
exclCelltypes <- c("NK", "Monocyte",  "Dendritic",  "Bcell")

# colors for plotting
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A", "#E28028", "#963C25", "grey")
names(pal) <- c("ALLBULK", "DN-like", "DNearly","SPCD4", "DP", "MSCBULK", "MONOBulk","DN3", "iSPCD4", "SPCD8", "gdTcell", "unknown")


## make heatmap
n = 5
for (patient in unique(srat_pts$orig.ident)) {
  #subset per patient
  Idents(srat_pts) <- "orig.ident"
  name = as.character(patient)
  srat_sub <- subset(srat_pts, idents = name)
  #remove the subsets with < ncells for better comparison
  Idents(srat_sub) <- "gating"
  levels(srat_sub) <- c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8", "NK", "Monocyte",  "Dendritic",  "Bcell", "unknown")
  idents.count <- as.data.frame(table(srat_sub@meta.data$gating))
  idents.keep  <- idents.count %>% group_by(Var1) %>% filter(Freq >= n) %>% filter(!Var1 %in% exclCelltypes)
  srat_sub <- subset(srat_sub, idents = idents.keep$Var1)
  #find DE genes per patient for clusters > n cells
  DefaultAssay(srat_sub) <- "Protein"
  srat_sub <- NormalizeData(srat_sub, normalization.method = "CLR", margin = 2)
  Idents(srat_sub) <- "gating"
  DEmarkers <- FindAllMarkers(srat_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
  write.csv(DEmarkers, paste0("PerPatient/DEmarkers_", name, ".csv"))
  ## prepare for heatmap
  top10 <- DEmarkers %>% group_by(cluster) %>% top_n(10, avg_log2FC) 
  markers <- intersect(top10$gene, rownames(srat_sub)) # safety that all markers are in the object
  # scale all proteins so that all proteins are in the scale.data slot
  all.proteins <- rownames(srat_sub) 
  srat_sub <- ScaleData(srat_sub, features = all.proteins)
  print(markers)
  # plot heatmap DE proteins
  markers <- markers[order(match(markers,prot2gene$Protein))] #order markers
  p1 <- DoHeatmap(srat_sub, features = markers, 
                  assay = "Protein", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal) + scale_fill_viridis()
  
  ## convert markers to genes to compare
  DE_df <- filter(prot2gene, Protein %in% markers)
  genes <- DE_df$Gene
  DefaultAssay(srat_sub) <- "RNA" # Used RNA because some genes are filtered out in SCT
  all.genes <- rownames(srat_sub)
  srat_sub <- ScaleData(srat_sub, features = all.genes)
  genes_filt <- intersect(genes, rownames(srat_sub@assays$RNA))
  
  # plot heatmap genes transcribing DE proteins
  genes_filt <- genes_filt[order(match(genes_filt, prot2gene$Gene))] #order genes same as markers
  p2 <- DoHeatmap(srat_sub, features = genes_filt, 
                  assay = "RNA", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) +  guides(color=FALSE)
  
  p_grid <- plot_grid(p1, p2)
  
  pdf(paste0("PerPatient/Pretty/hDEProtRNA_", name, ".pdf"), paper = "a4r", width = 14)
  print(p_grid)
  dev.off()
}
