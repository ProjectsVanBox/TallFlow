#!/bin/bash

## heatmap DE expression genes all patients

#Packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(sctransform)
devtools::install_github("xmc811/Scillus", ref = "development")
library(Scillus)

##
setwd('/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_ThyPtsComb/')

## load data set 
srat_pts <- readRDS('srat_pts.rds')
#srat_pts <- readRDS('../10x_SinglePt/Preprocessing/srat_10x_SinglePt.rds')

### Colors and names for coming analysis 
colors = list(SampleID = c("pt1179" = "#339999", "pt2337" = "#124D4D",
                           "pt344" = "#CC3366", "pt3291" = "#801438",
                           "pt4068" = "#99CC66", "pt5242" = "#44631F",
                           "pt2283" = "#FF9933", "pt11801" = "#8E4503",
                           "pt9160" = "#FFCC00", "pt315" = "#B29330",
                           "pt335" = "#003366", "pt9175" = "#7CE3D8",
                           "pt3045" = "#b0c4de", "pt5438" = "#6495ed",
                           "pt5676" = "#875692", "pt10138" = "#006666",
                           "Thy01" = "#ffa07a", "Thy02" = "#C2B280", "Thy03" = "#C71585"),
              PatientID = c("TALL1" = "#339999", "TALL2" = "#CC3366", 
                            "TALL3" = "#99CC66","TALL4" = "#FF9933", 
                            "TALL5" = "#FFCC00", "TALL6" = "#003366", 
                            "TALL7" = "#7CE3D8", "TALL8" = "#b0c4de",
                            "TALL9" = "#6495ed", "TALL10" = "#875692", 
                            "TALL11" = "#006666", "Thy01" = "ffa07a", 
                            "Thy02" = "#C2B280", "Thy03" = "C71585"),
              Site = c("BM" = "43C5F9", "PB" = "#FF0000"),
              Gender = c("M" = "#003366", "F" = "#CC3366"),
              Type = c("TALL" = "#052955", "THYMUS" = "#C59434"),
              Stage = c("D" = "#C59434", "R" = "#A3B7F9"),
              Clusters = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A",
                           "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD",
                           "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
                           "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"))

# colors for plotting Gating stages
pal = c("#33A02C", "#385725", "#E7298A", "#B17BA6", "#1965B0", "#E28028", "#963C25", "grey") 
names(pal) <- c("DNearly", "DN3", "iSPCD4", "DP", "SPCD4", "SPCD8", "gdTcell", "unknown")
colors_sampleID <- c("#006666", "#339999", "#FF9933", "#003366", "#CC3366", "#6495ed", "#875692","#7CE3D8")

#exclude non T-cell types
exclCelltypes <- c("NK", "Monocyte",  "Dendritic",  "Bcell")

# filter for only the patients that are in these paper. 
# only diagnosis
Idents(srat_pts) <- "SampleID"
srat_pts <- subset(srat_pts, 
                   idents = c("pt10138", "pt1179", "pt2283", "pt335",
                              "pt344", "pt5438", "pt5676", "pt9175"))

## downsample for script test
#srat_pts <- subset(srat_pts, downsample = 100)

## Umap of srat_pts
### RNA log normalization
DefaultAssay(srat_pts) <- "RNA"
srat_pts <- NormalizeData(srat_pts) 
srat_pts <- ScaleData(srat_pts)
srat_pts <- FindVariableFeatures(srat_pts)

### Run umap on the RNA set 
# perform visualization and clustering steps
srat_pts <- srat_pts %>%
  RunPCA(reduction.name = "pca", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", graph.name = "ssn") %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE) 
# Create plot
png("Umap_8pts_log.png")
print(DimPlot(srat_pts, group.by = "SampleID", reduction = "umap", cols = colors$SampleID) + ggtitle("by sample"))
dev.off()

## exclude non-T-cells.
DefaultAssay(srat_pts) <- "RNA"
Idents(srat_pts) <- "gating"
srat_pts <- subset(srat_pts, idents = exclCelltypes, invert = TRUE)

## find variable features
Idents(srat_pts) <- "gating"
DE_genes <- FindAllMarkers(srat_pts, only.pos = TRUE, min.pct = 0.25, 
              logfc.threshold = 0.25)
write.csv(DE_genes, "PerPatient/DEtesting/DEgenesPtsComblog.csv")
#DE_genes <- read.csv2("PerPatient/DEtesting/DEgenesPtsComblog.csv", header = TRUE, sep = ",")

## plot top 20 genes
top20 <- DE_genes %>% group_by(cluster) %>% top_n(20, avg_log2FC) 
genes <- intersect(top20$gene, rownames(srat_pts)) # safety that all genes are in the object

#ordering
my_levels <- c("DNearly", "DN3", "iSPCD4", "DP", "SPCD4", "SPCD8", "gdTcell","unknown")
srat_pts$gating <- factor(x = srat_pts$gating, levels = my_levels)
# Do heatmap with two upper bars. 
p1 <- plot_heatmap(dataset = srat_pts,
             markers = genes,
             sort_var = c("gating", "SampleID"),
             anno_var = c("gating","SampleID"), 
             anno_colors = list(pal, colors_sampleID))

png('DE_ptscombined_logDE.png')
print(p1)
dev.off()

pdf('DE_ptscombined_logDE.pdf')
print(p1)
dev.off()

## save srat_obj
saveRDS(srat_pts,"srat_8pts_log.rds")
## save plot object
saveRDS(p1, "DE_ptscombined_logDE.rds")

