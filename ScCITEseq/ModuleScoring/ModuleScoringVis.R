#!/bin/bash

## Module scoring dotplots visualizations from own healhty Thymus reference

## libraries 
## load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

## 
setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/')

# load data -> module scoring is already stored in the object
srat_pts <- readRDS('srat_pts.rds')

dotfeatures <- c("DNearly1", "DN31", "gdTcell1", "DP1", "SPCD41","SPCD81", "NK1", "Monocyte1", "Dendritic1","Bcell1", "unknown1")

# split by patient. 
n = 5  # Set the minimum count threshold
Idents(srat_pts) <- "gating"
for (patient in unique(srat_pts$orig.ident)) {
  Idents(srat_pts) <- "orig.ident"
  name = as.character(patient)
  srat_sub <- subset(srat_pts, idents = name)
  Idents(srat_sub) <- "gating"
  #srat_sub <- subset(srat_pts, idents = c("NK", "Monocyte",  "Dendritic",  "Bcell", "unknown"), invert= TRUE )
  levels(srat_sub) <- c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8","NK", "Monocyte",  "Dendritic",  "Bcell", "unknown" )
  idents.count <- as.data.frame(table(srat_sub@meta.data$gating))
  idents.keep  <- idents.count %>% group_by(Var1) %>% filter(Freq >= n) %>% filter(!Var1 %in% c("NK", "Monocyte",  "Dendritic",  "Bcell", "unknown"))
  p6 <- DotPlot(srat_sub, features = dotfeatures, idents = idents.keep$Var1, cols = "RdYlBu") + theme(axis.text.x = element_text(angle = 90))
  pdf(paste0(name,'_Modulescoring2.pdf'))
  print(p6)
  dev.off()
}


