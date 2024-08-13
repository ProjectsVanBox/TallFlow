#!/bin/bash


## predict identity of Cite annotated files based on Park transcriptional data

## libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(RColorBrewer)
#library(viridis)
library(cowplot)

#setwd('/Users/verapoort/hpc/projects/TallClonal/3_Output/Test_objects/')
setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/')

## load data
srat_pts <- readRDS('srat_pts.rds')
#srat_pts <- readRDS('pts_small.rds')
#thy.query <- readRDS('Thymi_small.rds')
thy.ref <- readRDS('srat_thymi.rds')

#thy.ref <- thy.query

#exclude non T-cell types
exclCelltypes <- c("NK", "Monocyte",  "Dendritic",  "Bcell")

# colors for plotting
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A", "#E28028", "#963C25", "grey")
names(pal) <- c("ALLBULK", "DN-like", "DNearly","SPCD4", "DP", "MSCBULK", "MONOBulk","DN3", "iSPCD4", "SPCD8", "gdTcell", "unknown")

sample_cols <- c("pt1179" = "#339999", "pt2337" = "#124D4D",
             "pt344" = "#CC3366", "pt3291" = "#801438",
             "pt4068" = "#99CC66", "pt5242" = "#44631F",
             "pt2283" = "#FF9933", "pt11801" = "#8E4503",
             "pt9160" = "#FFCC00", "pt315" = "#B29330",
             "pt335" = "#003366", "pt9175" = "#7CE3D8",
             "pt3045" = "#b0c4de", "pt5438" = "#6495ed",
             "pt5676" = "#875692", "pt10138" = "#006666",
             "Thy01" = "#ffa07a", "Thy02" = "#C2B280", "Thy03" = "#C71585")

#ordering
my_levels <- c("DNearly", "DN3", "iSPCD4", "DP", "SPCD4", "SPCD8", "gdTcell","unknown")

## prepare references data
DefaultAssay(thy.ref) <- "RNA" 
thy.ref <- thy.ref %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(dims = 1:30) %>% FindNeighbors() %>% FindClusters()
thy.ref <- RunUMAP(thy.ref, dims = 1:30)
#DimPlot(thy.ref, group.by = "gatingCelltype")

## prepare query data
DefaultAssay(srat_pts) <- "RNA"
srat_pts <- srat_pts %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(dims = 1:30) %>% FindNeighbors() %>% FindClusters() %>% RunUMAP(dims = 1:30)

## Cell type classification
thymus.anchors <- FindTransferAnchors(reference = thy.ref, 
                                      query = srat_pts, 
                                      dims = 1:30, 
                                      reduction = "pcaproject")

predictions <- TransferData(anchorset = thymus.anchors, 
                            refdata = thy.ref$gatingCelltype, 
                            dims = 1:30)

srat_pts <- AddMetaData(srat_pts, metadata = predictions)

srat_pts$prediction.match <- srat_pts$predicted.id == srat_pts$gating
print(table(srat_pts$prediction.match))


p1 <- DimPlot(srat_pts, reduction = "umap", group.by = c("gating", "predicted.id"), cols =pal)
p2 <- DimPlot(srat_pts, reduction = "umap", group.by = "SampleID", cols = sample_cols)


pdf('../../3_Output/10x_ThyPtsComb/celltypePrediction_thymus.pdf')
print(p1)
print(p2)
dev.off()

saveRDS(p1, '../../3_Output/10x_ThyPtsComb/celltypePredictionPlot1.rds')
saveRDS(p2, '../../3_Output/10x_ThyPtsComb/celltypepredictionPlot2.rds')

