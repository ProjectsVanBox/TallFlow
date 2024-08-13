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

setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb')

## load data
srat_pts <- readRDS('srat_pts.rds')
prot2gene <- read.csv('../../1_Input/ProteinGeneSymbols.csv', header = TRUE, sep = ";")

#exclude non T-cell types
exclCelltypes <- c("NK", "Monocyte",  "Dendritic",  "Bcell")

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
  ## convert markers to genes to compare
  DE_df <- filter(prot2gene, Protein %in% markers)
  genes <- DE_df$Gene
  genes_filt <- intersect(genes, rownames(srat_sub@assays$RNA))
  srat_sub <- NormalizeData(srat_sub, assay = "RNA")
  # plot Heatmaps
  p1 <- dittoHeatmap(srat_sub, markers, assay = "Protein", 
                     annot.by = c("gating", "orig.ident"),
                     order.by = "gating",
                     scaled.to.max = TRUE,
                     show_colnames = FALSE,
                     show_rownames = TRUE, 
                     complex = TRUE)
  
  p2 <- dittoHeatmap(srat_sub, genes_filt, assay = "RNA",
                       annot.by = c("gating", "orig.ident"),
                       order.by = "gating",
                       scaled.to.max = TRUE,
                       show_colnames = FALSE,
                       show_rownames = TRUE, 
                       complex = TRUE)
  
  pdf(paste0("PerPatient/heatmapDEProt_", name, "_logNorm.pdf"))
  print(p1)
  print(p2)
  dev.off()
}




