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
DE_thymi <- read.csv("DEgenes_thymi.csv")
## make heatmap

#nrow(DE_thymi)
top20 <- DE_thymi %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC) 

# Pick Genes
#exclude genes that are not present in both lists
genes <- intersect(top20$gene, rownames(srat_pts))

# Annotating and ordering cells by some meaningful feature(s):
p1 <- dittoHeatmap(srat_pts, genes,
             annot.by = c("gating", "orig.ident"),
             order.by = "gating",
             scaled.to.max = TRUE,
             show_colnames = FALSE,
             show_rownames = FALSE, 
             complex = TRUE)

pdf('heatmap_celltypes.pdf')
print(p1)
dev.off()

