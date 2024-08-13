#!/usr/bin/env Rscript

### Empty the variables
rm(list= ls(all=TRUE))
options(stringsAsFactors = FALSE)


## compare CITE seq data with 17c panel

### --- convert preprocessed seurat objects to reduce size, performed on the cluster

## libraries
library(Seurat)
library(flowCore)

setwd("/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_Pts_Seperate/Preprocessing/")

## load data
inputDir <- "/hpc/pmc_vanboxtel/projects/TallClonal/3_Output/10x_Pts_Seperate/Preprocessing/"
files <- list.files(path=inputDir, pattern="*.rds", full.names=TRUE, recursive=FALSE)

for (file in files) {
  srat <- readRDS(file)
  DefaultAssay(srat) <- "Protein"
  srat1 <- DietSeurat(srat, assays = c("Protein"))
  srat1 <- NormalizeData(srat1, assay = "Protein", normalization.method = "CLR", margin = 2)
  
  #df_srat <- as.data.frame(t(srat1@assays$Protein@data))
  SampleID <- srat1@meta.data$orig.ident[1]
  ## make fcs file for the 10X data
  #fr <- new("flowFrame", exprs = as.matrix(df_srat))
  #head(fr)
  
  #write.FCS(x = fr, filename = paste0(SampleID, "_CITE.fcs"))
  saveRDS(srat1, paste0(SampleID, "_sObject_prot.rds"))
}




##------ plotting (done locally):
library(ggplot2)
library(cowplot)
library(gridExtra)
library(egg)
library(grid)

## load data
inputDir <- "~/hpc/projects/TallClonal/3_Output/10x_Pts_Seperate/Preprocessing/"
files <- list.files(path=inputDir, pattern="*prot.rds", full.names=TRUE, recursive=FALSE)
outputdir <- "/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallFlow/4_Lab/CITE/"

for (file in files){
  srat <- readRDS(file)
  sampleID <- srat@meta.data$orig.ident[1]
  print(sampleID)

  p1 <- FeatureScatter(srat, feature1="CD5.1", feature2="CD7.1", group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend()
  
  p2 <- FeatureScatter(srat,  feature1="CD45", feature2="CD7.1", group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend()
  
  p3 <- FeatureScatter(srat, feature1="CD3",  feature2="CD7.1",group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend()
  
  p4 <- FeatureScatter(srat, feature1="CD1a", feature2="CD7.1", group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend()
  
  p5 <- FeatureScatter(srat, feature1="CD4.1", feature2="CD7.1", group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend()
  
  p6 <- FeatureScatter(srat, feature1="CD8", feature2="CD7.1", group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend()
  
  p7 <- FeatureScatter(srat, feature1="TCR-A-B", feature2="CD7.1", group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend()
  
  p8 <- FeatureScatter(srat, feature1="TCR-Y-O", feature2="CD7.1", group.by = "SampleID", cols = "white") + 
    geom_bin_2d(bins = 80) + scale_fill_distiller(palette = "Spectral") + labs(title = NULL) & NoLegend() 
  
  
  #grid1 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 1)
  
  pdf(paste0(outputdir, sampleID, "_CITE.pdf"), width = 7, height = 7)
  grid.arrange(grobs=lapply(list(
    p1, p2, p3, p4, p5, p6, p7, p8),
    set_panel_size), ncol = 2, top = textGrob(sampleID, gp=gpar(fontsize=20, font=3)))
  dev.off()

}

### ploting density plots of CITE data
library(ggridges)
library(reshape2)


features = c("TCR-Y-O", "TCR-A-B", "CD8", "CD4.1", "CD3","CD1a", "CD7.1", "CD45")

for (file in files){
  srat <- readRDS(file)
  sampleID <- srat@meta.data$orig.ident[1]
  print(sampleID)
  
  data <- as.data.frame(t(srat@assays$Protein@data))
  data1 <- data %>% select(all_of(features))
  data2 <- melt(data1)
  
  p1 <-ggplot(data2, aes(x=value, y = variable)) + 
    geom_density_ridges(fill="#69b3a2", color="#e9ecef", alpha=0.8, scale = 1) + 
    theme_classic() +
    ggtitle(sampleID) + xlab("Intensity")
  
  pdf(paste0(outputdir, sampleID, "_CITEDens.pdf"), width = 7, height = 7)
  print(p1)
  dev.off()
}

p1 <-ggplot(data2, aes(x=value, y = variable)) + 
  geom_density_ridges(fill="#69b3a2", color="#e9ecef", alpha=0.8, scale = 1) + 
  theme_classic() +
  ggtitle(sampleID) + xlab("Intensity")

