#!/bin/bash

## load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
BiocManager::install("dittoSeq")
library(dittoSeq)

setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/')

## load data
srat_pts <- readRDS('srat_pts.rds')
#DE_thymi <- read.csv("DEgenes_thymi.csv")
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
  idents.keep  <- idents.count %>% group_by(Var1) %>% filter(Freq >= n)
  srat_sub <- subset(srat_sub, idents = idents.keep$Var1)
  #find DE genes per patient for clusters > n cells
  Idents(srat_sub) <- "gating"
  DEgenes <- FindAllMarkers(srat_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(DEgenes, paste0("PerPatient/DEgenes_", name, ".csv"))
  ## prepare for heatmap
  top20 <- DEgenes %>% group_by(cluster) %>% top_n(20, avg_log2FC) 
  genes <- intersect(top20$gene, rownames(srat_sub)) # safety that all genes are in the object
  
  p1 <- dittoHeatmap(srat_sub, genes,
                     annot.by = c("gating", "orig.ident"),
                     order.by = "gating",
                     scaled.to.max = TRUE,
                     show_colnames = FALSE,
                     show_rownames = TRUE, 
                     complex = TRUE)
  pdf(paste0("PerPatient/heatmapDEgenes_", name, ".pdf"))
  print(p1)
  dev.off()
}




