
### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

## load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(randomcoloR)
library(egg)


setwd("~/Desktop")
getwd()

## colors 
cell_types_col <- c( "#124D4D", "#222222", "#F3C300", "#875692", "#F38400", "#A1CAF1",
                     "#BE0032", "#C2B280","#F99379", "#008856", "#E68FAC", "#0067A5",
                     "#999999", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#43C5F9", "cornflowerblue",  "#003366")

## load samples
sratThyB <- readRDS("/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallClonal/2_Code/ThymocyteTrajectory/srat_parkintegrated.rds")
sratThyP <- readRDS("/Users/verapoort/Downloads/srat_10xThyComb_ThymocytesAnnotated.rds")
#1 minuut laad tijd
# srat_parkintegrated is not the correct file
# srat_thyintegrated is not 
#plot 
sratThyP
sratThyP@meta.data$predicted.id
p1 <- DimPlot(sratThyP, reduction = "wnn.umap", group.by = "predicted.id", cols = cell_types_col )
p1

ptsrat <- readRDS("/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallClonal/3_Output/CellTypePrediction/simplify_celltypes/pt1179.RDS")

ptsrat
# the umap is strange, probably the coordinates of the combined analysis
# rerun umap script for 1 patient. 
p2 <- DimPlot(ptsrat, reduction = "wnn.umap")
p2

features = c("CD4.1", "CD19.1", "CD14.1", "CD56", "CD11C", "CD123", "CD8", "CD3", "CD1a", "CD7.1", "TCR-Y-O")
DefaultAssay(sratThyP) <- "Protein"

Idents(sratThyP) <- "predicted.id"
p3 <- RidgePlot(sratThyP, features = features[7:11] , ncol = 3)
p3
pdf('cutoffvaluesProtein.pdf')
p3
dev.off()
getwd()
#rownames(sratThyP@assays$Protein)

## do the asssignment:
type2 <- rep("unknown", nrow(sratThyP@meta.data))
names(type2) <- rownames(sratThyP@meta.data)

type2[WhichCells(sratThyP, expression = CD19.1 > 2)] <- "Bcell"
type2[WhichCells(sratThyP, expression = CD14.1 > 2)] <- "Monocyte"
type2[WhichCells(sratThyP, expression = CD56 > 2 )] <- "NK"
type2[WhichCells(sratThyP, expression = CD11C > 2)] <- "Dendritic"
type2[WhichCells(sratThyP, expression = CD123 > 2)] <- "Dendritic"
type2[WhichCells(sratThyP, expression = CD4.1 > 1 & CD8 > 1)] <- "DP"
type2[WhichCells(sratThyP, expression = CD4.1 > 1 & CD1a > 1.5 & CD3 < 0.5)] <- "iSPCD4"
type2[WhichCells(sratThyP, expression = CD4.1 > 1 & CD8 < 1)] <- "SPCD4"
type2[WhichCells(sratThyP, expression = CD8 > 1 & CD4.1 < 1)] <- "SPCD8"
type2[WhichCells(sratThyP, expression = TCR-Y-O > 1.8)] <- "gdTcell"
type2[WhichCells(sratThyP, expression = CD4.1 < 1 & CD1a > 1.5 & CD3 < 0.5 & CD8 < 1)] <- "DN3"
type2[WhichCells(sratThyP, expression = CD4.1 < 1 & CD1a < 1.5 & CD3 < 0.5 & CD8 < 1)] <- "DNearly"

sratThyP <- AddMetaData(sratThyP, col.name = "CITEanno", metadata = type2)

p4 <- DimPlot(sratThyP, reduction = "wnn.umap", group.by = "CITEanno", label = T, label.box = T)
p4

p5 <- DimPlot(sratThyP, reduction = "wnn.umap", group.by = "predicted.id", label = T, label.box = T)
p5

p6 <- FeaturePlot(sratThyP, reduction = "wnn.umap", features = features, order = T)
p6


# check module score 

# load markers of reference population

# Add module score
# Dotplot 

DE_park <- readRDS("/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallClonal/2_Code/ThymocyteTrajectory/DE_parkcelltypes.rds")
head(DE_park)

unique(DE_park$cluster)

nrow(DE_park)
top10 <- DE_park %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC) # 10 genes is too little as not all genes are found in my data set

DefaultAssay(sratThyP) <- "SCT"
for (pop in unique(DE_park$cluster)){
  pop_top10 <- top10[top10$cluster == as.character(pop),]
  sratThyP <- AddModuleScore(sratThyP, features = list(pop = pop_top10$gene), name = as.character(pop))
}

colnames(sratThyP@meta.data)

dotfeatures <- c("DN.early.1", "DN.P.1", "DN.Q.1", "αβT.entry.1","DP.Q.1", "DP.P.1", "CD8.T1", "CD8αα.II.1","CD8αα.I.1", "CD4.T1" , "γδT1")

Idents(sratThyP) <- "CITEanno"
levels(sratThyP) <- reverse(c("DNearly", "DN3", "iSPCD4", "DP", "SPCD4", "SPCD8", "unknown", "Dendritic", "NK", "Bcell"))
DotPlot(sratThyP, features = dotfeatures) + theme(axis.text.x = element_text(angle = 90))




