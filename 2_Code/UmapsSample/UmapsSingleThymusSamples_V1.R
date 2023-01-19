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


### Load data  
setwd("/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Umaps/Thymus/")
all.df <- readRDS("./CombinedThymi.rds")


### Loop over samples and perform single sample umaps
sampleNames <- unique(all.df$SampleIDshort)
for (s in sampleNames){
  print(s)
  ### Select data 
  single.df <- all.df[all.df$SampleIDshort == s, ][1:20]
  #colnames(single.df)
  
  ### Create directories
  dir.create(paste0("./SingleThymusUmap/", s))
  dir.create(paste0("./SingleThymusUmap/", s, "/MarkerPlots"))
  
  ### Run umaps
  umap.input <- single.df[,1:17] %>% select(- contains("CD45"))
  umap.data <- umap(umap.input, n_neighbors = perplex, n_threads = 12, n_epochs = n_iter)
  single.df %<>% mutate(umap1 = umap.data[, 1],  umap2 = umap.data[, 2])
  
  ### Plot
  p1 <- ggplot(single.df,  aes(x = umap1, y = umap2, color = Celltype)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw()+ scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  p2 <- ggplot(single.df,  aes(x = umap1, y = umap2, color = SampleID)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw()+  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  ggsave(filename = paste0("./SingleThymusUmap/", s, "/", s, "_Umap_Celltype.pdf"), p1)
  ggsave(filename = paste0("./SingleThymusUmap/", s, "/", s, "_Umap_SampleID.pdf"), p2)
  
  ### Louvain clustering for umap
  knn.umap = get.knn(as.matrix(umap.data), k = k)
  knn.umap = data.frame(from = rep(1:nrow(knn.umap$nn.index), k), to = as.vector(knn.umap$nn.index), weight = 1/(1 + as.vector(knn.umap$nn.dist)))
  nw.umap = graph_from_data_frame(knn.umap, directed = FALSE)
  nw.umap = simplify(nw.umap)
  lc.umap = cluster_louvain(nw.umap)
  single.df$umap.louvain = as.factor(membership(lc.umap))
  lc.umap.cent = single.df %>% group_by(umap.louvain) %>% select(umap1, umap2) %>% summarize_all(mean)
  p3 <- ggplot(single.df,  aes(x = umap1, y = umap2, color = umap.louvain)) +
    geom_point(size = 0.2, show.legend = TRUE) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
  ggsave(filename = paste0("./SingleThymusUmap/", s, "/", s, "_Umap_louvainclusters.pdf"), p3)
  
  ### Plot the expression of the markers 
  for (marker in lineage_markers){
    #print(marker)
    p = ggplot(single.df,  aes(x = umap1, y = umap2, color = single.df[, marker])) + geom_point(size = 0.3) + theme_bw() +
      scale_color_gradientn(marker, colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
    png(paste0("./SingleThymusUmap/", s, "/MarkerPlots/", s, "_umap_", marker, ".png"))
    print(p)
    dev.off()
  }
  
  ### Save data
  saveRDS(single.df, paste0("./SingleThymusUmap/", s, "/", s, "_data.rds"))
}




