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


#rsync -var --progress ~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/Umaps/Thymus/CombinedThymi.rds /Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/
all.df <- readRDS("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/CombinedThymi.rds")
all.df.sub <- all.df[all.df$Celltype %in% c("DN", "ISP", "gdTcell","DP", "CD4", "CD8"), ]
### Run umaps
umap.input <- all.df.sub[,1:17] %>% select(- contains("CD45"))
umap.data <- umap(umap.input, n_neighbors = perplex, n_threads = 12, n_epochs = n_iter)
all.df.sub %<>% mutate(umap1 = umap.data[, 1],  umap2 = umap.data[, 2])

p1 <- ggplot(all.df.sub,  aes(x = umap1, y = umap2, color = Celltype)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw()+ scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
p1
ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/Umap_Thymus_noHealthy.pdf", p1)

p2 <- ggplot(all.df.sub,  aes(x = umap1, y = umap2, color = SampleIDshort)) +geom_point(size = 0.2, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
p2
ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/Umap_Thymus_noHealthy_SampleID.pdf", p2)

for (s in unique(all.df.sub$SampleIDshort)){
  print(s)
  p.sample <- ggplot(subset(all.df.sub,SampleIDshort %in% c(s)),  aes(x = umap1, y = umap2)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  ggsave(filename = paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/SamplePlots_NoHealthy/Umap_Thymus_noHealthy_", s, ".pdf"), p.sample)
}
p.sample <- ggplot(subset(all.df.sub,SampleIDshort %in% c(s)),  aes(x = umap1, y = umap2)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
p.sample



#all.df$TissueType <- ifelse(grepl("Thy",all.df$SampleIDshort),'Thymus','Patient')
saveRDS(all.df.sub, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/Thymus_noHealthy.rds")
