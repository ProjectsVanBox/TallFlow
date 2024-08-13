#### To do trajectory analysis for thymocytes

#Load packages
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(DBI)
library(AnnotationDbi)
library(SingleCellExperiment)
library(tidyverse)
library(scales)
library(cowplot)
library(RCurl)
library(org.Hs.eg.db)
library(SeuratDisk)
library(matrixStats)
library(ggsci)
library(SingleR)
#library(pals)
library(batchelor)
library(BiocParallel)
library(BiocNeighbors)
library(celldex)

setwd("/hpc/pmc_vanboxtel/projects/TallClonal/2_Code/PreprocessingDataset/")
#setwd("/Users/ricohagelaar/Documents/Thymus_Project/SingleCellSequencing/TallClonal/2_Code/PreprocessingDataset/")

# prepare data and create seurat object
#dir.create("Output")
#dir.create("Plots")
inputdir = "../../1_Input/10x_SinglePt/"
outputdir = "../../3_Output/10x_SinglePt/Preprocessing/"
plotdir = paste0(outputdir, "Plots/")
project = "10x_SinglePt"

#Metadata
meta_data = as.data.frame(readxl::read_excel(paste0(inputdir,"CohortInfoScStudy.xlsx")))
genesToExclude <- read.csv(paste0(inputdir,"20210902_all_exclud_genes.csv"))

#Colors 
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
              Clusters = c("#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5", 
                           "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", 
                           
                           "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", "#99600F", "#B3823E", "#CCAA7A", "#E6D2B8", "#54990F", "#78B33E", 
                           "#A3CC7A", "#CFE6B8", "#0F8299", "#3E9FB3", "#7ABECC", "#B8DEE6", "#3D0F99", "#653EB3", "#967ACC", "#C7B8E6", 
                           "#333333", "#666666", "#999999", "#CCCCCC"))
# #ffa07a = "lightsalmon",  #C71585 = "mediumvioletred", #b0c4de = "lightsteelblue4",  #6495ed = "cornflowerblue", #FF0000= "red"
# kelly(n=22)[-1][-1] = "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26"
#
#In inputdir store cell ranger output filtered folder with count matrix files 
# Rename folder to pt"1234"  this name will be use in the analysis 
mtx_dirs = list.dirs(path= inputdir, 
                     recursive = F, 
                     full.names = TRUE)

count_files <- list()
srat_objects <- list()
for( i in 1:length(mtx_dirs)){
  filename <-gsub(paste0(inputdir, "/"), "", mtx_dirs[i])
  count_files[[filename]] <- append(count_files, Read10X(mtx_dirs[i]))
  name <- gsub(paste0(inputdir, "/"), "", mtx_dirs[i])
  srat_objects[[name]] <- CreateSeuratObject(counts = count_files[[i]]$`Gene Expression`, 
                                             assay = "RNA", strip.suffix = T)
  srat_objects[[name]][["Protein"]] <- CreateAssayObject(counts = count_files[[i]]$`Antibody Capture`)
}

# merge data sets to srat 
cell_id <- as.list(names(srat_objects))
srat_all <- srat_objects[[1]]#merge(x = srat_objects[[1]], y = srat_objects[-1], 
            #      add.cell.ids = cell_id, project = project)
srat_all$orig.ident <- sapply(X = strsplit(colnames(srat_all), split = "_"), FUN = "[", 1)
srat_all$SampleID <- sapply(X = strsplit(colnames(srat_all), split = "_"), FUN = "[", 1)

srat = srat_all
#str(srat)


#clean up global env
rm(count_files)
rm(srat_objects)
rm(srat_all)
gc()

# Add meta data
m <- match(srat$SampleID, meta_data$Sample_ID)
srat$PatientID = meta_data$Patient_ID[m]
srat$Stage = meta_data$Stage[m]
srat$Site = meta_data$Site[m]
srat$Gender = meta_data$Gender[m]
srat$Type = meta_data$Type[m]

#mitochondrial stress
srat$percent.mt <- PercentageFeatureSet(object = srat, pattern = "^MT-")

#UPR stress
GO_genes <- read.csv(file = paste0(inputdir, "GO_UPR.csv"), sep = ";")
GO_gene_list <- unique(GO_genes$Symbol)
GO_gene_list <- GO_gene_list[GO_gene_list %in% row.names(srat)]   
srat$percent.UPR <- PercentageFeatureSet(object = srat, features = GO_gene_list, assay = "RNA")

#likely some RBCs in the data set with low n genes but high n transcripts
hb_genes <- rownames(srat[["RNA"]]@data)[grep(pattern = "^HB", rownames(srat[["RNA"]]@data))]
hb_counts <- Matrix::colSums(srat@assays$RNA@counts[hb_genes,])
srat <- AddMetaData(srat, col.name = "log2hb_genes", 
                    metadata = log2(1 +hb_counts))
srat <- AddMetaData(srat, col.name = "pct_hemo", metadata = 100 *
                      hb_counts/srat@meta.data$nCount_RNA)

#quantile(srat$pct_hemo)
#median(srat$pct_hemo)

#colnames(srat@meta.data)

# Stats
# cell counts
#table(srat@meta.data$PatientID)
Idents(srat) <- "SampleID"

p1 <- plot_grid(VlnPlot(srat, features=c('nCount_RNA'), pt.size = 0, y.max = 30000, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"), 
                VlnPlot(srat, features=c('nFeature_RNA'), pt.size = 0, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),
                VlnPlot(srat, features=c('percent.mt'), pt.size = 0, y.max = 40, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),
                VlnPlot(srat, features=c('percent.UPR'), pt.size = 0, y.max = 10, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),  
                VlnPlot(srat, features=c('nCount_Protein'), pt.size = 0, y.max = 25000,cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),
                VlnPlot(srat, features=c('nFeature_Protein'), pt.size = 0, y.max = 220, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"))

#p1 

png(paste0(plotdir,"Features.png"))
p1 
dev.off()
### Testing area 
#png(paste0(plotdir,"Test_Features.png"))
#t1 
#dev.off()
#is.na(srat$percent.mt)
#packageVersion("Seurat")

data = srat@meta.data

p2 = ggplot(data, aes(x=nCount_RNA, y= percent.mt))+
  geom_point(shape=16, size=1)+
  scale_fill_igv()+
  coord_trans(x = "log10")+
  scale_y_continuous(breaks= scales::pretty_breaks(n=10))+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_light()

p3 = ggplot(data, aes(x=nCount_RNA, y= nFeature_RNA))+
  scale_fill_igv()+
  geom_point(shape=16, size=1)+
  coord_trans(x = "log10", y = "log10")+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_light()
#p2 + p3

png(paste0(plotdir,"Quality.png"))
plot_grid(p2, p3) 
dev.off()

v <- VlnPlot(srat, features = "pct_hemo", cols = colors$SampleID) & NoLegend()
#v

# Filter low quality cells
srat$drop = !c(srat$nCount_RNA <= 500 | srat$nFeature_RNA <= 250 | srat$percent.mt >= 30
               | srat$percent.UPR >= 4 | srat$nCount_RNA >= 40000 | srat$nFeature_RNA >= 7000
               | srat$nCount_Protein >= 20000 | srat$pct_hemo >= 2)


data = srat@meta.data

p4 = ggplot(data, aes(x=nCount_RNA, y= percent.mt, color = drop))+
  geom_point(shape=16, size=1)+
  scale_fill_igv()+
  coord_trans(x = "log10")+
  scale_y_continuous(breaks= scales::pretty_breaks(n=10))+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_color_manual(values=c("red", "black")) +
  theme_light() & NoLegend()

p5 = ggplot(data, aes(x=nCount_RNA, y= nFeature_RNA, color = drop))+
  scale_fill_igv()+
  geom_point(shape=16, size=1)+
  coord_trans(x = "log10", y = "log10")+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_color_manual(values=c("red", "black")) +
  theme_light() & NoLegend()
#p4 + p5

png(paste0(plotdir,"Filtered.png"))
plot_grid(p4, p5) 
dev.off()

srat_filtered = subset(srat, subset = drop)

# cell counts before and after
#ncol(srat)
#ncol(srat_filtered)

# filter low quality genes -> This is only done for the lognormalization
# for SCT normalization it is not needed 
nexprs = rowSums(srat_filtered@assays$RNA@counts > 0)
keep = rownames(srat_filtered)[nexprs >= 3]

srat_filtered2 = subset(srat_filtered, features = keep)
srat_filtered[["RNA"]] <- CreateAssayObject(counts = srat_filtered2@assays$RNA@counts)

#clean up env
rm(keep)
rm(srat_filtered2)
rm(data)

# gene counts before and after
#nrow(srat)
#nrow(srat_filtered)

srat <- srat_filtered

# Prepare Protein data
DefaultAssay(srat) <- "Protein"
# perform visualization and clustering steps
# Alter resolution of FindClusters if too many or too few clusters
srat <- srat %>% NormalizeData(normalization.method = "CLR") %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA(reduction.name = "apca", verbose = FALSE) %>%
  FindNeighbors(reduction = "apca", graph.name = "Protein_snn") %>%
  FindClusters(reduction = "apca", resolution = 0.8, graph.name = "Protein_snn", verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = "apca", reduction.name = "ab_umap", reduction.key = "ABUMAP_", verbose = FALSE) 


# FastMNN batch effect correction for Protein assay
#MNN = reducedMNN(srat@reductions$apca@cell.embeddings,
#                 batch = srat$orig.ident,
#                 BPPARAM = MulticoreParam(workers=8),
#                 BNPARAM = HnswParam())
#srat[["MNNProtein"]] = CreateDimReducObject(embeddings = MNN$corrected,
#                                            assay="Protein",
#                                            key="ProtMNN_")
#srat <- RunUMAP(srat, dims = 1:15, 
#                verbose=FALSE, 
#                reduction = "MNNProtein", 
#                reduction.key="ab_UMAPMNN_", 
#                reduction.name = "ab_UMAPMNN")

protein_features <- c("CD8", "CD4.1", "CD3", "CD1a", "CD20", "CD19.1", "CD14.1",
                      "CD56", "CD7.1", "CD5.1", "CD45", "CD44.1", "CD34.1", "CD117", "CD2.1",
                      "CD38.1", "CD95")


p6b <- FeaturePlot(srat, features = protein_features, reduction = "ab_umap")
p6 <- DimPlot(srat, label = TRUE, reduction = "ab_umap", cols  = colors$Clusters) + ggtitle("by cluster")
p6a <- DimPlot(srat, group.by = "SampleID", reduction = "ab_umap", cols = colors$SampleID) + ggtitle("by sample")
#p6a + p6 

png(paste0(plotdir,"UmapProtein.png"))
plot_grid(p6a, p6) 
dev.off()

png(paste0(plotdir,"UmapProteinFeatures.png"))
plot_grid(p6b) 
dev.off()

# Prepare RNA data
DefaultAssay(srat) <- "RNA"
#Normalization and scaling
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat)
srat <- ScaleData(srat)
srat <- RunPCA(srat, reduction.name = "pca", verbose = FALSE)
srat <- FindNeighbors(srat, reduction = "pca", graph.name = "ssn")
srat <- FindClusters(srat, reduction = "pca", resolution = 0.8, graph.name = "ssn", verbose = FALSE)
var.features <- setdiff(VariableFeatures(srat), genesToExclude$x)
srat <- RunUMAP(srat, features=var.features, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)

#srat <- RunUMAP(srat, umap.method = 'umap-learn', metric = 'correlation', features=var.features, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)

#srat <- srat %>% NormalizeData() %>% 
#  FindVariableFeatures() %>% 
#  ScaleData() %>%
#  RunPCA(reduction.name = "pca", verbose = FALSE) %>%
#  FindNeighbors(reduction = "pca", graph.name = "ssn") 
#features=setdiff(VariableFeatures(srat), genesToExclude$x),
#srat@assays$RNA@var.features

#srat <- srat %>%
#  FindClusters(reduction = "pca", resolution = 0.8, graph.name = "ssn", verbose = FALSE) %>%
#  RunUMAP(features=setdiff(VariableFeatures(srat), genesToExclude$x), reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE) 

p8 <- DimPlot(srat, label = TRUE, reduction = "umap", cols = colors$Clusters) + ggtitle("by cluster")
p8a <- DimPlot(srat, group.by = "SampleID", reduction = "umap", cols = colors$SampleID) + ggtitle("by sample")
#p8a + p8
p8b <- FeaturePlot(srat, features = protein_features)

png(paste0(plotdir,"UmapRNA.png"))
plot_grid(p8a, p8) 
dev.off()

png(paste0(plotdir,"UmapRNAFeatures.png"))
plot_grid(p8b) 
dev.off()

#SCTransform RNA data 
srat <- suppressWarnings(SCTransform(srat, vars.to.regres = NULL,
                                     verbose = FALSE, variable.features.n = 3000))
srat <- srat %>% FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA(reduction.name = "pca", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", graph.name = "ssn") %>%
  FindClusters(reduction = "pca", resolution = 0.8, graph.name = "ssn", verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE) 

p9 <- DimPlot(srat, label = TRUE, reduction = "umap", cols = colors$Clusters) + ggtitle("by cluster")
p9a <- DimPlot(srat, group.by = "SampleID", reduction = "umap", cols = colors$SampleID) + ggtitle("by sample")
#p9a + p9
p9b <- FeaturePlot(srat, features = protein_features)
#p9b


#Batch correct between thymus samples using FastMNN


DefaultAssay(srat) <- "SCT"

#MNN = reducedMNN(srat@reductions$pca@cell.embeddings,
#                 batch = srat$orig.ident,
#                 BPPARAM = MulticoreParam(workers=8),
#                 BNPARAM = HnswParam())
#srat[["MNN"]] = CreateDimReducObject(embeddings = MNN$corrected,
#                                     assay="RNA",
#                                     key="MNN_")
#srat <- RunUMAP(srat, dims = 1:15, 
#                verbose=FALSE, 
#                reduction = "MNN", 
#                reduction.key="UMAPMNN_", 
#                reduction.name = "UMAPMNN")

#srat <- FindNeighbors(srat, dims = 1:15, reduction = "MNN", graph.name = "MNN_snn")
#srat <- FindClusters(srat, reduction = "MNN", graph.name = "MNN_snn")

#p9 = DimPlot(srat, group.by="SampleID", reduction = "UMAPMNN", cols = colors$SampleID) + ggtitle("Batch corrected")
#p9a = DimPlot(srat, group.by="SampleID", reduction = "umap", cols = colors$SampleID) + ggtitle("not corrected")
#p10 = DimPlot(srat, group.by="MNN_snn_res.0.8", reduction = "UMAPMNN", cols = colors$Clusters) + ggtitle("Batch corrected")
#p9 + p9a
#p10
#png(paste0(plotdir,"UmapRNAFeatures.png"))
#plot_grid(p9a, p9) 
#dev.off()





# WNN 
srat <- FindMultiModalNeighbors(srat, 
                                reduction.list = list("pca","apca"), 
                                dims.list = list(1:20,1:15),
                                modality.weight.name = "RNA.weight")

srat <- RunUMAP(srat, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap", 
                reduction.key = "wnnUMAP_")

srat <- FindClusters(srat, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
srat <- FindClusters(srat, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)
srat <- FindClusters(srat, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)



# save data
saveRDS(srat, paste0(outputdir,"srat_10x_SinglePt.rds"))


