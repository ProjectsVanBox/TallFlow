#!/bin/bash

## integrating data from previously published Cordes et al.

## libraries 
## load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

## 
setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/')

## load data 
## Cordes data 
srat_ref <- readRDS('/hpc/pmc_vanboxtel/projects/TallClonal/1_Input/ExternalDatasets/Cordes/GSE195812_All_thymus_samples.rds')

# check what is in this data 
#colnames(srat_ref@meta.data)

p1 <- DimPlot(srat_ref, reduction = "UMAP", group.by = 'orig.ident')

pdf('CordesInt/UmapReference.pdf')
p1
dev.off()

## integration methods? Or just extracting a gene score and then using this for the single cell patients. 

# find DE genes
Idents(srat_ref) <- "orig.ident"
DEgenes <- FindAllMarkers(srat_ref, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DEgenes, "CordesInt/DEgenes_cordes.csv")

## module scoring
## load data:
srat_pts <- readRDS('srat_pts.rds')

top20 <- DEgenes %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC) # 10 genes is too little as not all genes are found in my data set

# make sure population names and module score names will not overwrite.
for (pop in unique(DEgenes$cluster)){
  pop_top20 <- top20[top20$cluster == as.character(pop),]
  srat_pts <- AddModuleScore(srat_pts, features = list(pop = pop_top20$gene), name = as.character(pop))
  srat_ref <- AddModuleScore(srat_ref, features = list(pop = pop_top20$gene), name = as.character(pop))
}

#  cell populations
dotfeatures <- c("DN11", "DN21", "DN31", "ISP1", "DP_CD3min1","DP_CD3plus1", "CD41", "CD81")

# plot reference plot
Idents(srat_ref) <- "orig.ident"
levels(srat_ref) <- c("DN1", "DN2", "DN3", "ISP", "DP_CD3min","DP_CD3plus", "CD4", "CD8")
p5 <- DotPlot(srat_ref, features = dotfeatures) + theme(axis.text.x = element_text(angle = 90))
pdf('CordesInt/ModuleScoring_reference.pdf')
print(p5)
dev.off()

# plot patients bulk
Idents(srat_pts) <- "gating"
levels(srat_pts) <- c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8", "NK","Monocyte",  "Dendritic",  "Bcell","unknown" )
p6 <- DotPlot(srat_pts, features = dotfeatures) + theme(axis.text.x = element_text(angle = 90))
pdf('CordesInt/ModuleScoring_ptsbulk.pdf')
print(p6)
dev.off()


# split by patient. 
#Idents(srat_pts) <- "gating"
n = 5
for (patient in unique(srat_pts$orig.ident)) {
  Idents(srat_pts) <- "orig.ident"
  name = as.character(patient)
  srat_sub <- subset(srat_pts, idents = name)
  Idents(srat_sub) <- "gating"
  levels(srat_sub) <- c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8", "NK", "Monocyte",  "Dendritic",  "Bcell", "unknown")
  idents.count <- as.data.frame(table(srat_sub@meta.data$gating))
  idents.keep <- idents.count %>% group_by(Var1) %>% filter(Freq >= n)
  p7 <- DotPlot(srat_sub, features = dotfeatures, idents = idents.keep$Var1) + theme(axis.text.x = element_text(angle = 90))
  pdf(paste0('CordesInt/',name,'_Modulescoring.pdf'))
  print(p7)
  dev.off()
}

#write module scoring to csv file that can later be added to a seurat object using AddMetaData
# this save a lot of storage space!
srat_pts_modulesCordes <- srat_pts[[dotfeatures]]
write.csv(srat_pts_modulesCordes, 'CordesInt/srat_pts_moduleScores.csv')

cordes_moduleScore <- srat_ref[[dotfeatures]]
write.csv(cordes_moduleScore, 'CordesInt/Cordes_modulescores.csv')


#saveRDS(srat_pts, 'CordesInt/srat_ptsCordes.rds')
