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
inputdir = "../../1_Input/10x_Pts/"
outputdir = "../../3_Output/10x_Pts_Seperate/Preprocessing/"
plotdir = paste0(outputdir, "Plots/")
project = "Pts_Seperate"


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


# Load in previous set
#srat <- readRDS("/Users/ricohagelaar/Documents/Thymus_Project/SingleCellSequencing/TallClonal/2_Code/PreprocessingDataset/srat_thymocytes.rds")
srat <- readRDS("/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_Pts_30dims/Preprocessing/srat_10xPts.rds")
split_seurat<- SplitObject(srat, split.by = "orig.ident")


# Prepare data 
for (i in 1:length(split_seurat)) {
  # First dimensional reduction
  split_seurat[[i]] <- RunPCA(split_seurat[[i]])
  split_seurat[[i]] <- RunUMAP(split_seurat[[i]], dims = 1:30, verbose = F)
  split_seurat[[i]] <- FindNeighbors(split_seurat[[i]], verbose = F)
  split_seurat[[i]] <- FindClusters(split_seurat[[i]], resolution = 0.8, verbose = F)
  
  File.Name <- unique(split_seurat[[i]]$orig.ident)
  dir.create(paste0(plotdir,  File.Name))
  # First plots
  p1 <- DimPlot(split_seurat[[i]], 
                reduction = "umap"# , group.by = "Phase")
  ) + ggtitle(File.Name) 
  png(paste0(plotdir, File.Name, "/", File.Name, "_umap_30dims.png"))
  print(p1)
  dev.off()
  
  ### Skip the samples that make the mysterious error: disk quota exceeded 122. 
  # And the samples that already finished :D
  if (File.Name %in% c("pt10138","pt1179","pt11801","pt2283","pt2337","pt3045","pt315","pt3291","pt335","pt344","pt4068","pt5242","pt5438","pt5676")){
    print(paste0("Skipping: ", File.Name))
  }else{
    # Prepare Protein data
    DefaultAssay(split_seurat[[i]]) <- "Protein"
    # perform visualization and clustering steps
    # Alter resolution of FindClusters if too many or too few clusters
    split_seurat[[i]] <- split_seurat[[i]] %>% NormalizeData(normalization.method = "CLR") %>% 
      FindVariableFeatures() %>% 
      ScaleData() %>%
      RunPCA(reduction.name = "apca", verbose = FALSE) %>%
      FindNeighbors(reduction = "apca", graph.name = "Protein_snn") %>%
      FindClusters(reduction = "apca", resolution = 0.8, graph.name = "Protein_snn", verbose = FALSE) %>%
      RunUMAP(dims = 1:30, reduction = "apca", reduction.name = "ab_umap", reduction.key = "ABUMAP_", verbose = FALSE) 
    
    protein_features <- c("CD8", "CD4.1", "CD3", "CD1a", "CD20", "CD19.1", "CD14.1",
                          "CD56", "CD7.1", "CD5.1", "CD45", "CD44.1", "CD34.1", "CD117", "CD2.1",
                          "CD38.1", "CD95")
    p2b <- FeaturePlot(split_seurat[[i]], features = protein_features, reduction = "ab_umap")
    p2 <- DimPlot(split_seurat[[i]], label = TRUE, reduction = "ab_umap", cols  = colors$Clusters) + ggtitle("by cluster")
    p2a <- DimPlot(split_seurat[[i]], group.by = "SampleID", reduction = "ab_umap", cols = colors$SampleID) + ggtitle("by sample")
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_UmapProtein.png"))
    print(plot_grid(p2a, p2) )
    dev.off()
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_UmapProteinFeatures.png"))
    print(plot_grid(p2b) )
    dev.off()
    
    # Prepare RNA data
    DefaultAssay(split_seurat[[i]]) <- "RNA"
    #Normalization and scaling
    split_seurat[[i]] <- NormalizeData(split_seurat[[i]])
    split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]])
    split_seurat[[i]] <- ScaleData(split_seurat[[i]])
    split_seurat[[i]] <- RunPCA(split_seurat[[i]], reduction.name = "pca", verbose = FALSE)
    split_seurat[[i]] <- FindNeighbors(split_seurat[[i]], reduction = "pca", graph.name = "ssn")
    split_seurat[[i]] <- FindClusters(split_seurat[[i]], reduction = "pca", resolution = 0.8, graph.name = "ssn", verbose = FALSE)
    var.features <- setdiff(VariableFeatures(split_seurat[[i]]), genesToExclude$x)
    split_seurat[[i]] <- RunUMAP(split_seurat[[i]], features=var.features, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)
    #features=setdiff(VariableFeatures(split_seurat[[i]]), genesToExclude$x),
    
    p3 <- DimPlot(split_seurat[[i]], label = TRUE, reduction = "umap", cols = colors$Clusters) + ggtitle("by cluster")
    p3a <- DimPlot(split_seurat[[i]], group.by = "SampleID", reduction = "umap", cols = colors$SampleID) + ggtitle("by sample")
    #p8a + p8
    p3b <- FeaturePlot(split_seurat[[i]], features = protein_features)
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_UmapRNA.png"))
    print(plot_grid(p3a, p3) )
    dev.off()
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_UmapRNAFeatures.png"))
    print(plot_grid(p3b) )
    dev.off()
    
    #SCTransform RNA data 
    split_seurat[[i]] <- suppressWarnings(SCTransform(split_seurat[[i]], vars.to.regres = NULL,
                                                      verbose = FALSE, variable.features.n = 3000))
    split_seurat[[i]] <- split_seurat[[i]] %>% FindVariableFeatures() %>% 
      ScaleData() %>%
      RunPCA(reduction.name = "pca", verbose = FALSE) %>%
      FindNeighbors(reduction = "pca", graph.name = "ssn") %>%
      FindClusters(reduction = "pca", resolution = 0.8, graph.name = "ssn", verbose = FALSE) %>%
      RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE) 
    
    p4 <- DimPlot(split_seurat[[i]], label = TRUE, reduction = "umap", cols = colors$Clusters) + ggtitle("by cluster")
    p4a <- DimPlot(split_seurat[[i]], group.by = "SampleID", reduction = "umap", cols = colors$SampleID) + ggtitle("by sample")
    p4b <- FeaturePlot(split_seurat[[i]], features = protein_features)
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_UmapRNA_SCT.png"))
    print(plot_grid(p4a, p4) )
    dev.off()
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_UmapRNAFeatures_SCT.png"))
    print(plot_grid(p4b) )
    dev.off()
    
    # WNN 
    split_seurat[[i]] <- FindMultiModalNeighbors(split_seurat[[i]], 
                                                 reduction.list = list("pca","apca"), 
                                                 dims.list = list(1:20,1:15),
                                                 modality.weight.name = "RNA.weight")
    
    split_seurat[[i]] <- RunUMAP(split_seurat[[i]], nn.name = "weighted.nn", 
                                 reduction.name = "wnn.umap", 
                                 reduction.key = "wnnUMAP_")
    
    split_seurat[[i]] <- FindClusters(split_seurat[[i]], graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
    split_seurat[[i]] <- FindClusters(split_seurat[[i]], graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)
    split_seurat[[i]] <- FindClusters(split_seurat[[i]], graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)
    
    p5 <- DimPlot(split_seurat[[i]], label = TRUE, reduction = "wnn.umap", cols = colors$Clusters) + ggtitle("by cluster")
    p5a <- DimPlot(split_seurat[[i]], label = TRUE, reduction = "wnn.umap", cols = colors$Clusters, group.by = "wsnn_res.0.2") + ggtitle("wsnn_res.0.2")
    p5b <- DimPlot(split_seurat[[i]], label = TRUE, reduction = "wnn.umap", cols = colors$Clusters, group.by = "wsnn_res.0.5") + ggtitle("wsnn_res.0.5")
    p5c <- DimPlot(split_seurat[[i]], label = TRUE, reduction = "wnn.umap", cols = colors$Clusters, group.by = "wsnn_res.1") + ggtitle("wsnn_res.1")
    
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_Umap_wsnn.png"))
    print(plot_grid(p5) )
    dev.off()
    
    png(paste0(plotdir, File.Name, "/", File.Name, "_Umap_wsnn_res.png"))
    print(plot_grid(p5a, p5b, p5c) )
    dev.off()
    
    
    
    # save data
    saveRDS(split_seurat[[i]], paste0(outputdir, File.Name, ".rds"))
  }
}













#for (i in 1:length(split_seurat)) {
#  File.Name <- unique(split_seurat[[i]]$orig.ident)
#  print(File.Name)
#  if (File.Name %in% c("Thy01", "Thy03")){
#    print(paste0("Skipping: ", File.Name))
#  }else{
#    print(paste0("Processing: ", File.Name))
#  }
#}












