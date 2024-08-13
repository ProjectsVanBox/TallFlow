### Warning! 
# RunUmap script needs to be run before starting this script


# Load packages
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
library(batchelor)
library(BiocParallel)
library(BiocNeighbors)
library(celldex)
library(viridis)
library(pheatmap)


### Change directory 
setwd("/hpc/pmc_vanboxtel/projects/TallClonal/2_Code/CellTypeIdent/")


### Prepare data structure and read in info files 
inputdir = "../../1_Input/10x_Pts/"
outputdir = "../../3_Output/10x_PtsComb/CellTypeIdent/"
plotdir = paste0(outputdir, "Plots/")
project = "Pts_Combined"
# Metadata
meta_data = as.data.frame(readxl::read_excel(paste0(inputdir,"CohortInfoScStudy.xlsx")))


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
              Clusters = c("#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5", 
                           "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", 
                           
                           "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", "#99600F", "#B3823E", "#CCAA7A", "#E6D2B8", "#54990F", "#78B33E", 
                           "#A3CC7A", "#CFE6B8", "#0F8299", "#3E9FB3", "#7ABECC", "#B8DEE6", "#3D0F99", "#653EB3", "#967ACC", "#C7B8E6", 
                           "#333333", "#666666", "#999999", "#CCCCCC"))


### Read in data 
srat <- readRDS("../../3_Output/10x_PtsComb/RunUmap/srat_10xPtsComb_RunUmap.rds")


### Load in healthy reference sets
#ref_mon <- MonacoImmuneData()
#ref_nov <- NovershternHematopoieticData()
#ref_blue <- BlueprintEncodeData()
#saveRDS(ref_mon, "ref_mon.rds")
#saveRDS(ref_nov, "ref_nov.rds")
#saveRDS(ref_blue, "ref_blue.rds")
ref_mon <- readRDS("ref_mon.rds")
ref_nov <- readRDS("ref_nov.rds")
ref_blue <- readRDS("ref_blue.rds")


### Perform SingleR
singler_mon <- SingleR(test = GetAssayData(srat, assay = "SCT",slot = "data"), ref = ref_mon, labels = ref_mon$label.main)
singler_nov <- SingleR(test = GetAssayData(srat, assay = "SCT",slot = "data"), ref = ref_nov, labels = ref_nov$label.main)
singler_blue <- SingleR(test = GetAssayData(srat, assay = "SCT",slot = "data"), ref = ref_blue, labels = ref_blue$label.main)


### Create plots 
pdf(paste0(plotdir,"SingleR_SCT_Monaco.pdf"))
plotScoreHeatmap(singler_mon, show.labels = TRUE, max.labels = 100, show.pruned = FALSE, order.by = "clusters", clusters = srat@meta.data[["ssn_res.1.2"]])
dev.off()
pdf(paste0(plotdir,"SingleR_SCT_Novershtern.pdf"))
plotScoreHeatmap(singler_nov, show.labels = TRUE, max.labels = 100, show.pruned = FALSE, order.by = "clusters", clusters = srat@meta.data[["ssn_res.1.2"]])
dev.off()
pdf(paste0(plotdir,"SingleR_SCT_BlueprintEncode.pdf"))
plotScoreHeatmap(singler_blue, show.labels = TRUE, max.labels = 100, show.pruned = FALSE, order.by = "clusters", clusters = srat@meta.data[["ssn_res.1.2"]])
dev.off()


# Resolution is probably between 0.8 and 1.4. 1 looks good




