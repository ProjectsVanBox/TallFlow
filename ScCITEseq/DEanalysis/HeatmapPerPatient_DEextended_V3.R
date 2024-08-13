#!/bin/bash

## load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)

setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/')

## load data
srat_pts <- readRDS('srat_pts.rds')

# thymi version
#srat_obj <- readRDS('srat_thymi.rds')
#DE_thymi <- read.csv("DEgenes_thymi.csv")
#cellannotation <- read.csv('cellgating.csv', header = TRUE, row.names = 1)
#colnames(cellannotation) <- "gating"

# add meta data of cell annotation to thymi object
#srat_obj <- AddMetaData(srat_obj, cellannotation, col.name = "gating")
#srat_pts <- srat_obj

# colors for plotting
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A", "#E28028", "#963C25", "grey")
names(pal) <- c("ALLBULK", "DN-like", "DNearly","SPCD4", "DP", "MSCBULK", "MONOBulk","DN3", "iSPCD4", "SPCD8", "gdTcell", "unknown")

plot_order <- c("DNearly", "DN3", "iSPCD4", "DP", "SPCD4", "SPCD8", "gdTcell", "unknown")
plot_order_rev <- rev(plot_order)

#exclude non T-cell types
exclCelltypes <- c("NK", "Monocyte",  "Dendritic",  "Bcell")

# exclude genes from DE analysis: sex genes, mitochondrial genes, ribosomal genes. 
exclgenes <- read.csv('../../1_Input/10x_Pts/20210902_all_exclud_genes.csv', header = TRUE)
#exclgenes <- read.csv('~/hpc/projects/TallClonal/1_Input/10x_Pts/20210902_all_exclud_genes.csv', header = TRUE)


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
  ## DEFAULT SETTING = Wilcoxon Rank sum test
  # filter to used only variable features, and no mitochondrial/ ribosomal genes
  Idents(srat_sub) <- "gating"
  srat_sub <- FindVariableFeatures(srat_sub, assay = "SCT")
  gene.list <- srat_sub@assays$SCT@var.features
  gene.list.filt <- as.data.frame(gene.list) %>% filter(!gene.list %in% exclgenes$x)
  #gene.list.filt <- gene.list[grep("^MT", gene.list, invert=TRUE)]
  #gene.list.filt <- gene.list.filt[grep("^RPL", gene.list.filt, invert=TRUE)]
  #gene.list.filt <- gene.list.filt[grep("^RPS", gene.list.filt, invert=TRUE)]
  DEgenes <- FindAllMarkers(srat_sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, features = gene.list.filt$gene.list)
  write.csv(DEgenes, paste0("PerPatient/DEtestingNew/DEgenesWilC_", name, ".csv"))
  ## prepare for heatmap
  top20 <- DEgenes %>% group_by(cluster) %>% top_n(20, avg_log2FC) 
  genes <- intersect(top20$gene, rownames(srat_sub)) # safety that all genes are in the object
  
  p1 <- DoHeatmap(srat_sub, features = genes, 
                  assay = "SCT", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal) + 
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
  
  
  top10 <- DEgenes %>% group_by(cluster) %>% top_n(10, avg_log2FC) 
  genes10 <- intersect(top10$gene, rownames(srat_sub)) # safety that all genes are in the object
  
  Idents(srat_sub) <- "gating"
  levels(srat_sub) <- plot_order
  p2 <- DotPlot(srat_sub, features = genes10,
                assay = "SCT", cols = "RdYlBu") + coord_flip()
  
  top5 <- DEgenes %>% group_by(cluster) %>% top_n(5, avg_log2FC) 
  genes5 <- intersect(top5$gene, rownames(srat_sub)) # safety that all genes are in the object
  
  p3 <- DoHeatmap(srat_sub, features = genes5, 
                  assay = "SCT", group.by = "gating", 
                  group.bar=TRUE, group.colors = pal) + 
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
  
  
  top3 <- DEgenes %>% group_by(cluster) %>% top_n(3, avg_log2FC) 
  genes3 <- intersect(top3$gene, rownames(srat_sub)) # safety that all genes are in the object
  
  Idents(srat_sub) <- "gating"
  levels(srat_sub) <- plot_order
  p4 <- DotPlot(srat_sub, features = genes3,
                assay = "SCT", cols = "RdYlBu") + coord_flip()
  
  
  pdf(paste0("PerPatient/DEtestingNew/DEplots_", name, ".pdf"))
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
}


#print("Heatmaps and Dotplots are finished, start pseudobulk visualization")

#DEgenes.list <- list.files("PerPatient/DEtesting/", recursive = FALSE, pattern = ".csv")




