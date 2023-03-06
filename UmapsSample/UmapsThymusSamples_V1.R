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
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/")
all.df <- readRDS("./CombinedThymi.rds")


### Subselect data 
thy070 <- all.df[all.df$SampleID %in% c("Thy070_depleted", "Thy070_undepleted"), ]
thy072 <- all.df[all.df$SampleID %in% c("Thy072_depleted", "Thy072_undepleted"), ]
thy073 <- all.df[all.df$SampleID %in% c("Thy073_depleted", "Thy073_undepleted"), ]


### Clean up old variables in order to reconstruct new ones based on the subseleted data 
thy070 <- thy070 %>% select(-starts_with("umap1")) %>% select(-starts_with("umap2"))  %>% select(-starts_with("umap.louvain"))
thy072 <- thy072 %>% select(-starts_with("umap1")) %>% select(-starts_with("umap2"))  %>% select(-starts_with("umap.louvain"))
thy073 <- thy073 %>% select(-starts_with("umap1")) %>% select(-starts_with("umap2"))  %>% select(-starts_with("umap.louvain"))


### Run umap thy070
dir.create("./Thy070")
thy070.umap.input <- thy070[,1:17] %>% select(- contains("CD45"))
thy070.umap.data <- umap(thy070.umap.input, n_neighbors = perplex, n_threads = 12, n_epochs = n_iter)
thy070 %<>% mutate(umap1 = thy070.umap.data[, 1],  umap2 = thy070.umap.data[, 2])
thy070.p1 <- ggplot(thy070, aes(x = umap1, y = umap2, color = Celltype)) + geom_point(size = 0.2, show.legend = TRUE) + theme_bw() + scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
thy070.p2 <- ggplot(thy070, aes(x = umap1, y = umap2, color = SampleID)) + geom_point(size = 0.2, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "./Thy070/Thy070_Umap_Celltype.pdf", thy070.p1)
ggsave(filename = "./Thy070/Thy070_Umap_SampleID.pdf", thy070.p2)
### Louvain clustering for umap thy070
thy070.knn.umap = get.knn(as.matrix(thy070.umap.data), k = k)
thy070.knn.umap = data.frame(from = rep(1:nrow(thy070.knn.umap$nn.index), 
                                 k), to = as.vector(thy070.knn.umap$nn.index), weight = 1/(1 + as.vector(thy070.knn.umap$nn.dist)))
thy070.nw.umap = graph_from_data_frame(thy070.knn.umap, directed = FALSE)
thy070.nw.umap = simplify(thy070.nw.umap)
thy070.lc.umap = cluster_louvain(thy070.nw.umap)
thy070$umap.louvain = as.factor(membership(thy070.lc.umap))
thy070.lc.umap.cent = thy070 %>% group_by(umap.louvain) %>% select(umap1, umap2) %>% summarize_all(mean)
plot(thy070$umap1, thy070$umap2, pch=16, cex=0.3, col=thy070$umap.louvain)
pdf("./Thy070/Thy070_Umap_louvainclusters.pdf")
ggplot(thy070,  aes(x = umap1, y = umap2, color = umap.louvain)) +
  geom_point(size = 0.2, show.legend = TRUE) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
dev.off()
### Save data 
saveRDS(thy070, "./Thy070/Thy070.rds")


### Run umap thy072
dir.create("./Thy072")
thy072.umap.input <- thy072[,1:17] %>% select(- contains("CD45"))
thy072.umap.data <- umap(thy072.umap.input, n_neighbors = perplex, n_threads = 12, n_epochs = n_iter)
thy072 %<>% mutate(umap1 = thy072.umap.data[, 1],  umap2 = thy072.umap.data[, 2])
thy072.p1 <- ggplot(thy072, aes(x = umap1, y = umap2, color = Celltype)) + geom_point(size = 0.2, show.legend = TRUE) + theme_bw() + scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
thy072.p2 <- ggplot(thy072, aes(x = umap1, y = umap2, color = SampleID)) + geom_point(size = 0.2, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "./Thy072/Thy072_Umap_Celltype.pdf", thy072.p1)
ggsave(filename = "./Thy072/Thy072_Umap_SampleID.pdf", thy072.p2)
### Louvain clustering for umap thy072
thy072.knn.umap = get.knn(as.matrix(thy072.umap.data), k = k)
thy072.knn.umap = data.frame(from = rep(1:nrow(thy072.knn.umap$nn.index), 
                                        k), to = as.vector(thy072.knn.umap$nn.index), weight = 1/(1 + as.vector(thy072.knn.umap$nn.dist)))
thy072.nw.umap = graph_from_data_frame(thy072.knn.umap, directed = FALSE)
thy072.nw.umap = simplify(thy072.nw.umap)
thy072.lc.umap = cluster_louvain(thy072.nw.umap)
thy072$umap.louvain = as.factor(membership(thy072.lc.umap))
thy072.lc.umap.cent = thy072 %>% group_by(umap.louvain) %>% select(umap1, umap2) %>% summarize_all(mean)
plot(thy072$umap1, thy072$umap2, pch=16, cex=0.3, col=thy072$umap.louvain)
pdf("./Thy072/Thy072_Umap_louvainclusters.pdf")
ggplot(thy072,  aes(x = umap1, y = umap2, color = umap.louvain)) +
  geom_point(size = 0.2, show.legend = TRUE) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
dev.off()
### Save data 
saveRDS(thy072, "./Thy072/Thy072.rds")


### Run umap thy073
dir.create("./Thy073")
thy073.umap.input <- thy073[,1:17] %>% select(- contains("CD45"))
thy073.umap.data <- umap(thy073.umap.input, n_neighbors = perplex, n_threads = 12, n_epochs = n_iter)
thy073 %<>% mutate(umap1 = thy073.umap.data[, 1],  umap2 = thy073.umap.data[, 2])
thy073.p1 <- ggplot(thy073, aes(x = umap1, y = umap2, color = Celltype)) + geom_point(size = 0.2, show.legend = TRUE) + theme_bw() + scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
thy073.p2 <- ggplot(thy073, aes(x = umap1, y = umap2, color = SampleID)) + geom_point(size = 0.2, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "./Thy073/Thy073_Umap_Celltype.pdf", thy073.p1)
ggsave(filename = "./Thy073/Thy073_Umap_SampleID.pdf", thy073.p2)
### Louvain clustering for umap thy073
thy073.knn.umap = get.knn(as.matrix(thy073.umap.data), k = k)
thy073.knn.umap = data.frame(from = rep(1:nrow(thy073.knn.umap$nn.index), 
                                        k), to = as.vector(thy073.knn.umap$nn.index), weight = 1/(1 + as.vector(thy073.knn.umap$nn.dist)))
thy073.nw.umap = graph_from_data_frame(thy073.knn.umap, directed = FALSE)
thy073.nw.umap = simplify(thy073.nw.umap)
thy073.lc.umap = cluster_louvain(thy073.nw.umap)
thy073$umap.louvain = as.factor(membership(thy073.lc.umap))
thy073.lc.umap.cent = thy073 %>% group_by(umap.louvain) %>% select(umap1, umap2) %>% summarize_all(mean)
plot(thy073$umap1, thy073$umap2, pch=16, cex=0.3, col=thy073$umap.louvain)
pdf("./Thy073/Thy073_Umap_louvainclusters.pdf")
ggplot(thy073,  aes(x = umap1, y = umap2, color = umap.louvain)) +
  geom_point(size = 0.2, show.legend = TRUE) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
dev.off()
### Save data 
saveRDS(thy073, "./Thy073/Thy073.rds")












