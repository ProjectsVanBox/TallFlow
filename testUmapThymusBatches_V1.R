#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(flowWorkspace)
library(flowCore)
library(DescTools)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(data.table)
library(dendextend)
library(uwot)
library(magrittr)
library(Rtsne)
library(ggforce)
library(ggridges)
library(RColorBrewer)
library(FNN)
library(igraph)
#library(ConsensusClusterPlus)
#library(FlowSOM)


### function for winsorizing pct_level can be adjusted to alter winsorizing treshold (we checked for this dataset and 1% looks much better than 5%)
wins_vars <- function(x, pct_level = 0.005){
  if(is.numeric(x)){
    Winsorize(x, probs = c(pct_level, 1-pct_level), na.rm = T)
  } else {x}
}


### Function to plot the histograms for the expression of markers 
plot_clustering_distr_wrapper <- function(expr, cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(list(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering, 
                            labels = paste0(levels(cell_clustering), " (", freq_clust, "%)"))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr), 
              id.vars = "cluster", value.name = "expression", 
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster, 
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster, 
                                             color = reference, fill = reference), alpha = 0.3, scale = 0.9) +
    facet_grid_paginate( ~ antigen, scales = "free_y", nrow = 2, page = 2) +
    theme_bw() +
    theme(axis.text = element_text(size = 5, angle = 45), 
          strip.text = element_text(size = 7), legend.position = "none") 
}


### create flow Frame from Data frame 
# from: https://rdrr.io/github/ssmpsn2/flowAssist/src/R/DFtoFF.R
DFtoFF<-function(DF){
  if(class(DF) == "data.frame"){
    return(flowFrame(as.matrix(DF)))
  }
  if(class(DF) == "list"){
    frameList<-as.list(NULL)
    length(frameList)<-length(DF)
    for(i in 1:length(DF)){
      if(class(DF[[i]]) == "data.frame"){
        frameList[[i]]<-flowFrame(as.matrix(DF[[i]]))
        names(frameList)[[i]]<-names(DF)[[i]]
      }
      else{
        warning(paste("Object at index",i,"not of type data.frame"))
      }
    }
    return(frameList)
  }
  else {
    stop(" Object is not of type data.frame")
  }
}


#### function for Biex transformation
trans <- flowjo_biexp(channelRange = 4096, maxValue = 300000, 
                      pos =  4.33, neg = 1, widthBasis = -10, inverse = FALSE)

### Local variables 
n = 20000
n_iter = 20000
perplex = 50
lying_iter = 5000
ncell = 5000 
lineage_markers <- c("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7",
                     "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "TCRg_d")

### Colour selection 
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN-like", "ISP-like", "gdTcell-like","DP-like", "CD4-like", "CD8-like", "cDC-like", 
                "B-cell-like", "Monocyte-like", "NK-cell-like", "pDC-like", 
                "Unknown")

### read in data  
#file_dir <- "/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PythonGating/Patients/"
file_dir <- "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/"
setwd(file_dir)

csvfiles <- c(#"/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A4 Thy02_depleted_stained_CD45+.csv", 
              "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A6 Thy071_thawed_stained_CD45+.csv",
              "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A6 Thy072_70Mdepleted_stained_CD45+.csv",
              "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A6 Thy075_10Mdepleted_stained_CD45+.csv",
              "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A4 Thy076_70Mdepleted_stained_CD45+.csv",
              "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A12 Thy076_thawed_stained_CD45+.csv",
              "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A4 Thy079_depleted_stained_CD45+.csv",
              "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A4 Thy085_Depleted_stained_CD45+.csv")



#"/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/"

all.df <- {}
for (f in csvfiles){
  # Get the sample names
  SampleName = gsub(".*/", "", gsub("_stained_CD45+.csv", "", f))
  SampleNameShort <- gsub(".* ", "", SampleName)
  print(SampleName)
  print(SampleNameShort)
  
  data.file <- read.table(file = f, sep = ",", header = TRUE, row.names = NULL)
  data.file$SampleID <- SampleName
  data.file$SampleIDshort <- SampleNameShort
  
  ### Select data 
  df <- data.file %>% select(matches('__'))
  colnames(df) <- gsub(x = gsub(colnames(df), pattern = ".*__", replacement = ""), pattern = ".A", replacement = "")
  df <- df %>% select(-starts_with("Celltype")) %>% select(-starts_with("DeadCells")) %>% 
    select(- contains("invNK_T")) %>% select(- contains("DeadCells")) %>% select(- contains("CD127")) #%>% select(- contains("CD45"))
  #df$Celltype <- data.file$Celltype
  df$SampleID <- data.file$SampleID
  df$SampleIDshort <- data.file$SampleIDshort
  
  ### Winsorize dataset in bulk (seems better if winsorizing is done before transformation!)
  df[] <- bind_cols(lapply(df, wins_vars))
  win.all.df <- df
  
  ### Biex transform data
  trans.all.df <- win.all.df[,1:17]
  trans.all.df[] <- bind_cols(lapply(trans.all.df, trans))
  trans.all.df[,18:19] <- df[,18:19] 
  
  ### Downsample 
  if(nrow(trans.all.df) < ncell){
    data.sub <- trans.all.df[ sample(nrow(trans.all.df), nrow(trans.all.df)), ]
  }else{
    data.sub <- trans.all.df[ sample(nrow(trans.all.df), ncell), ]
  }
  ### Combine files 
  all.df <- rbind(all.df, data.sub)
}



### Run umaps
umap.input <- all.df[,1:17] %>% select(- contains("CD45"))
umap.data <- umap(umap.input, n_neighbors = perplex, n_threads = 12, n_epochs = n_iter)
all.df %<>% mutate(umap1 = umap.data[, 1],  umap2 = umap.data[, 2])
#p1 <- ggplot(all.df,  aes(x = umap1, y = umap2, color = Celltype)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw()+ scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
p2 <- ggplot(all.df,  aes(x = umap1, y = umap2, color = SampleID)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw()+  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
p2




ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/TestThymusBatches/Thymusbatch_Umap_SampleID.pdf", p2)
saveRDS(all.df, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/TestThymusBatches/Thymusbatch.rds")




all.df <- readRDS("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/TestThymusBatches/Thymusbatch.rds")


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/TestThymusBatches/Thymusbatch_Umap_SampleID.pdf", height = 7, width = 14)
ggplot(all.df,  aes(x = umap1, y = umap2, color = SampleID)) +geom_point(size = 0.1, show.legend = TRUE) + theme_bw()+  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
dev.off()
