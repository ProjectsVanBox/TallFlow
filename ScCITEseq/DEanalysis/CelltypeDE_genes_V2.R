

###--- To make cell type labels comparable and get the Differential Expressed genes per cell type for module scoring.
## subset on thymus samples so that the genes taken into account are for healthy conditions.

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

# mem needed, at least 40GB

## libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(cutoff)
#library(SeuratDisk)
library(sctransform)

#load Seurat objects
srat <- readRDS("/Users/verapoort/hpc/projects/TallClonal/3_Output/10x_ThyPtsComb/ThyptscombALL.rds")
# sratCordes <- readRDS('...')

## cell gating was stored as csv.
anno <- read.csv('cellgating.csv', row.names = TRUE)
colnames(anno) <- "gatingCelltype"

# add anno to seurat object
srat <- AddMetaData(srat, anno)

#subset to cells of healthy thymi donors
Idents(srat) <- "orig.ident"
srat_thymi <- subset(srat, idents = c("Thy01", "Thy02", "Thy03"))

#save subset
saveRDS(srat_thymi, "srat_thymi.rds")
#srat_thymi <- readRDS('/Users/verapoort/hpc/projects/TallClonal/3_Output/10x_ThyPtsComb/srat_thymi.rds')
# SCT transform
srat_thymi <- SCTransform(srat_thymi, vst.flavor = "v2") %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)


# find DE genes
Idents(srat_thymi) <- "gatingCelltype"
# FInd All markers
DEgenes <- FindAllMarkers(srat_thymi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DEgenes, "DEgenes_thymi.csv")


p1 <- DimPlot(srat_thymi, reduction = "umap") 
p2 <- DimPlot(srat_thymi, reduction = "umap", group.by = "orig.ident")
pdf("srat_thymi.pdf")
p1
p2
dev.off()
# the gating does not seem perfect on the SCT transformed umap. 
# maybe addin the CITE umap as well will improve this? also the SCTransform umap does not look great. 

DefaultAssay(srat_thymi) <- "Protein"
srat_thymi <- srat_thymi %>%  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA(reduction.name = "apca", verbose = FALSE) %>%
  FindNeighbors(reduction = "apca", graph.name = "Protein_snn") %>%
  RunUMAP(dims = 1:30, reduction = "apca", reduction.name = "ab_umap", reduction.key = "ABUMAP_", verbose = FALSE) 
# CITE umap is sotred as srat_thymiCite -> does look like UMAPS generated in earlier stages of the project. 
# with nice CD4/CD8 distribution. 

# read the seurat patients only object to integrate these
# based on vignette: https://satijalab.org/seurat/articles/sctransform_v2_vignette
Idents(srat) <- "orig.ident"
srat_pts <- subset(srat, idents = c("Thy01", "Thy02", "Thy03"), invert = TRUE)

srat_pts <- SCTransform(srat_pts, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
saveRDS(srat_pts, 'srat_pts.rds')

srat.list <- list(ctrl = srat_thymi, stim = srat_pts)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(srat.list, anchor.features = features)

srat.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
                                         anchor.features = features)
srat.combined.sct <- IntegrateData(anchorset = srat.anchors, normalization.method = "SCT")

## DE analysis

immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)

b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "B_STIM", ident.2 = "B_CTRL",
                                     verbose = FALSE)
head(b.interferon.response, n = 15)




# module scoring: 

## the transformation is not performed yet!! so make sure to do this. 

DE_thymi <- read.csv("DEgenes_thymi.csv")
head(DE_thymi)

unique(DE_thymi$cluster)

nrow(DE_thymi)
top20 <- DE_thymi %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC) # 10 genes is too little as not all genes are found in my data set

# check the refence population vs de module scoring as control
#DefaultAssay(srat_thymi) <- "SCT"
for (pop in unique(DE_thymi$cluster)){
  pop_top20 <- top20[top20$cluster == as.character(pop),]
  srat_thymi <- AddModuleScore(srat_thymi, features = list(pop = pop_top20$gene), name = as.character(pop))
}

colnames(srat_thymi@meta.data)

dotfeatures <- c("DNearly1", "DN31", "gdTcell1", "DP1", "SPCD41","SPCD81", "NK1", "Monocyte1", "Dendritic1","Bcell1", "unknown1")


Idents(srat_thymi) <- "gatingCelltype"
levels(srat_thymi) <- c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8", "NK","Monocyte",  "Dendritic",  "Bcell","unknown" )
p4 <- DotPlot(srat_thymi, features = dotfeatures) + theme(axis.text.x = element_text(angle = 90))
pdf('ModuleScoring_thymicontr.pdf')
p4
dev.off()

# check the patient population vs de module scoring
#srat_pts <- readRDS('srat_pts.rds')

DE_thymi <- read.csv("DEgenes_thymi.csv")
head(DE_thymi)

nrow(DE_thymi)
top20 <- DE_thymi %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC) # 10 genes is too little as not all genes are found in my data set

for (pop in unique(DE_thymi$cluster)){
  pop_top20 <- top20[top20$cluster == as.character(pop),]
  srat_pts <- AddModuleScore(srat_pts, features = list(pop = pop_top20$gene), name = as.character(pop))
}

dotfeatures <- c("DNearly1", "DN31", "gdTcell1", "DP1", "SPCD41","SPCD81", "NK1", "Monocyte1", "Dendritic1","Bcell1", "unknown1")

Idents(srat_pts) <- "gating"
levels(srat_pts) <- c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8", "NK","Monocyte",  "Dendritic",  "Bcell","unknown" )
p5 <- DotPlot(srat_pts, features = dotfeatures) + theme(axis.text.x = element_text(angle = 90))
pdf('ModuleScoring_ptsbulk.pdf')
print(p5)
dev.off()


# split by patient. 
n = 5  # Set the minimum count threshold
Idents(srat_pts) <- "gating"
for (patient in unique(srat_pts$orig.ident)) {
  Idents(srat_pts) <- "orig.ident"
  name = as.character(patient)
  srat_sub <- subset(srat_pts, idents = name)
  Idents(srat_sub) <- "gating"
  levels(srat_sub) <- c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8", "NK", "Monocyte",  "Dendritic",  "Bcell", "unknown")
  idents.count <- as.data.frame(table(srat_sub@meta.data$gating))
  idents.keep  <- idents.count %>% group_by(Var1) %>% filter(Freq >= n)
  p6 <- DotPlot(srat_sub, features = dotfeatures, idents = idents.keep$Var1) + theme(axis.text.x = element_text(angle = 90))
  pdf(paste0(name,'_Modulescoring.pdf'))
  print(p6)
  dev.off()
}

#-> the DF now only contains the values to keep.
# use these to find the cells in the Seurat obj to exclude. 
# OR SPECIFY THESE IDENTS IN THE DOTPLOT !!!!
# use the DF as identiti to plot in the DotPlot. 
p6 <- DotPlot(srat_sub, features = dotfeatures, idents = df) + theme(axis.text.x = element_text(angle = 90))


#WhichCells(object = object, ident.remove = "ident.remove")	WhichCells(object = object, idents = "ident.remove", invert = TRUE)




#with standard log-normalization we get the exact same Dotplots

p7 <- DoHeatmap(subset(srat_sub, downsample = 100), features = dotfeatures, size = 3)




