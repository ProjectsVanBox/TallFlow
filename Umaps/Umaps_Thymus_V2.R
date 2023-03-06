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
library(ConsensusClusterPlus)
library(FlowSOM)


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


#### function for Biex transformation
trans <- flowjo_biexp(channelRange = 4096, maxValue = 300000, 
                      pos =  4.33, neg = 1, widthBasis = -10, inverse = FALSE)


### Local variables 
n = 20000
n_iter = 20000
perplex = 50
lying_iter = 5000
ncell = 5000 
k = 100
lineage_markers <- c("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7",
                     "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "TCRg_d")

### Colour selection 
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                "B-cell","Monocyte", "NK-cell", "pDC", 
                "Unknown")


### Read in data  
#file_dir <- "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/ThymusMore/"
file_dir <- "/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PythonGating/ThymusMore/"
setwd(file_dir)
csvfiles <- dir(".", pattern="*csv$")

all.df <- {}
for (f in csvfiles){
  # Get the sample names
  SampleName = gsub("_keymarkers.csv", "", f)
  SampleNameShort <- gsub("_.*", "", gsub(".*pt", "pt", SampleName))
  print(SampleName)
  print(SampleNameShort)
  
  ## Read in data 
  data.file <- read.table(file = f, sep = ",", header = TRUE, row.names = NULL)
  data.file$SampleID <- SampleName
  data.file$SampleIDshort <- SampleNameShort
  
  ### Select data 
  df <- data.file %>% select(matches('__'))
  colnames(df) <- gsub(x = gsub(colnames(df), pattern = ".*__", replacement = ""), pattern = ".A", replacement = "")
  df <- df %>% select(-starts_with("Celltype")) %>% select(-starts_with("DeadCells")) %>% 
    select(- contains("invNK_T")) %>% select(- contains("DeadCells")) %>% select(- contains("CD127")) #%>% select(- contains("CD45"))
  df$Celltype <- data.file$Celltype
  df$SampleID <- data.file$SampleID
  df$SampleIDshort <- data.file$SampleIDshort
  
  ### Winsorize dataset in bulk (seems better if winsorizing is done before transformation!)
  df[] <- bind_cols(lapply(df, wins_vars))
  win.all.df <- df
  
  ### Biex transform data
  trans.all.df <- win.all.df[,1:17]
  trans.all.df[] <- bind_cols(lapply(trans.all.df, trans))
  trans.all.df[,18:20] <- df[,18:20] 
  
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

### Plot
p1 <- ggplot(all.df,  aes(x = umap1, y = umap2, color = Celltype)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw()+ scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
p2 <- ggplot(all.df,  aes(x = umap1, y = umap2, color = SampleID)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw()+  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "../../../3_Output/Umaps/Thymus/CombinedUmap_Celltype.pdf", p1)
ggsave(filename = "../../../3_Output/Umaps/Thymus/CombinedUmap_SampleID.pdf", p2)


### Louvain clustering for umap
knn.umap = get.knn(as.matrix(umap.data), k = k)
knn.umap = data.frame(from = rep(1:nrow(knn.umap$nn.index), 
                                 k), to = as.vector(knn.umap$nn.index), weight = 1/(1 + as.vector(knn.umap$nn.dist)))
nw.umap = graph_from_data_frame(knn.umap, directed = FALSE)
nw.umap = simplify(nw.umap)
lc.umap = cluster_louvain(nw.umap)
all.df$umap.louvain = as.factor(membership(lc.umap))
lc.umap.cent = all.df %>% group_by(umap.louvain) %>% select(umap1, umap2) %>% summarize_all(mean)
plot(all.df$umap1, all.df$umap2, pch=16, cex=0.3, col=all.df$umap.louvain)
pdf("../../../3_Output/Umaps/Thymus/CombinedUmap_louvainclusters.pdf")
ggplot(all.df,  aes(x = umap1, y = umap2, color = umap.louvain)) +
  geom_point(size = 0.2, show.legend = TRUE) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
dev.off()


### Plot the expression of the markers 
for (marker in lineage_markers){
  print(marker)
  p = ggplot(all.df,  aes(x = umap1, y = umap2, color = all.df[, marker])) + geom_point(size = 0.3) + theme_bw() +
    scale_color_gradientn(marker, colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
  png(paste0("../../../3_Output/Umaps/Thymus/MarkerPlots/Combined_umap_", marker, ".png"))
  print(p)
  dev.off()
}


### Save data
saveRDS(all.df, "../../../3_Output/Umaps/Thymus/CombinedThymi.rds")








