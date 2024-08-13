## Run locally
# on merged thymus incl non thymocytes! + T-ALL patient samples
# make sure data is CLR margin 2 normalized

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

## libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(cutoff)
library(SeuratDisk)

#pbmc <- LoadH5Seurat("/Users/verapoort/Downloads/multi.h5seurat")

#rownames(pbmc@assays$ADT)

#rownames(pbmc@assays$RNA)
setwd('/Users/verapoort/hpc/projects/TallClonal/3_Output/10x_ThyPtsComb/')

# load data
srat <- readRDS('thyptscomb_protein.rds')
#srat <- readRDS('/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallClonal/3_Output/CellTypePrediction/simplify_celltypes/pt315.RDS')
srat

srat <- NormalizeData(srat, normalization.method = "CLR", margin = 2)
colnames(pbmc@meta.data)
unique(pbmc@meta.data$celltype.l3)
# plot to check if data is correct
nrow(pbmc@assays$ADT)
DimPlot(srat, reduction = "ab_umap", group.by = "SampleID")
srat

#CITE features
features = c("CD4.1", "CD19.1", "CD14.1", "CD56", "CD11C", "CD123", "CD8", "CD3", "CD1a", "CD7.1", "TCR-Y-O")

features_hao = c("CD4-1", "CD4-2", "CD19", "CD14", "CD11c", "CD123", "CD8a", "CD8", "CD3-1", "CD3-2", "CD1a", "CD7.1", "TCR-Y-O", "CD56-1", "CD56-2")
RidgePlot(pbmc, features = features[7:12] , ncol = 3)

# get cutoff values
Idents(srat) <- "SingleR"
p1 <- RidgePlot(srat, features = features[1:6] , ncol = 3)
p2 <- RidgePlot(srat, features = features[7:11] , ncol = 3)


pdf('histvaluesProtein.pdf')
print(p1)
print(p2)
dev.off()


# annotate using cut-off values

# cut-off annotation:
DefaultAssay(srat) <- "Protein"
type2 <- rep("unknown", nrow(srat@meta.data))
names(type2) <- rownames(srat@meta.data)

type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD19.1") > 1)] )] <- "Bcell"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD14.1") > 0.75)] )] <- "Monocyte"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD56") > 1)] )] <- "NK"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD11C") > 1)] )] <- "Dendritic"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD123") > 2)] )] <- "Dendritic"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD8") > 1 & FetchData(srat, vars = "CD4.1") > 1.5) ] )] <- "DP"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD4.1") > 1.5 & FetchData(srat, vars = "CD1a") > 1.5 & FetchData(srat, vars = "CD3") < 0.1 )] )] <- "iSPCD4"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD8") < 1 & FetchData(srat, vars = "CD4.1") > 1.5) ] )] <- "SPCD4"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD8") > 1 & FetchData(srat, vars = "CD4.1") < 1.5) ] )] <- "SPCD8"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "TCR-Y-O") > 1.7)] )] <- "gdTcell"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD4.1") < 1.5 & FetchData(srat, vars = "CD1a") > 1.5 & FetchData(srat, vars = "CD8") < 1 )] )] <- "DN3"
type2[colnames(srat@assays$Protein[,which(FetchData(srat, vars = "CD4.1") < 1.5 & FetchData(srat, vars = "CD1a") < 1.5 & FetchData(srat, vars = "CD8") < 1 & FetchData(srat, vars = "CD7.1") > 0.5 )] )] <- "DNearly"

srat <- AddMetaData(srat, col.name = "gating", metadata = type2)

write.csv(type2, "cellgating.csv")
getwd()
p2 <- DimPlot(srat, reduction = "umap", group.by = "gating")

png('umapSCTGating.png')
p2
dev.off()



#---- testing using an automatic cut-off generator

# make binary using Finite mixture models: 
#http://marcchoisy.free.fr/fmm/index.html
for (marker in features){
  names <- marker
  marker_out <- em(srat@assays$Protein@data[marker],"normal","normal")
  marker_out
  #confint(marker_out,nb=100,level=.95)
  cut_off <- cutoff(marker_out)
  
  hist(srat@assays$Protein@data["CD4"],breaks = 20, xlab = paste0(marker," expression"))
  #pdf(paste0(marker, "_cutoffHist.png"))
  print(p1)
  print(abline(v=cut_off[-1],lty=2, col = "blue"))
  print(abline(v=cut_off[1],lty=2, col = "blue"))
  #dev.off()
}

hist(srat@assays$Protein@data["CD4"],breaks = 20, xlab = paste0(marker," expression"))

features <- "CD4.1"

for (marker in features){
  names <- marker
  print(paste0(marker, "_cutoffHist.png"))
}

marker_out <- em(srat@assays$Protein@data["CD3",],"normal","normal")
confint(marker_out,nb=100,level=.95)

pdf('...pdf')
hist(srat@assays$Protein@data["3",],breaks = 20, xlab = "CD3expression")
lines(marker_out,lwd=1.5,col="red")
dev.off()

cut_off <- cutoff(marker_out)

hist(srat@assays$Protein@data["CD3",],breaks = 20, xlab = "CD56expression")
abline(v=cut_off[-1],lty=2, col = "blue")
abline(v=cut_off[1],lty=2, col = "blue")


hist(srat@assays$Protein@scale.data["CD56",])
