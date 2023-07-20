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
pal2 = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal2) <- c("DN-like", "ISP-like", "gdTcell-like","DP-like", "CD4-like", "CD8-like", "cDC-like", 
                "B-cell-like", "Monocyte-like", "NK-cell-like", "pDC-like", 
                "Unknown")
### Colour selection 
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                "B-cell","Monocyte", "NK-cell", "pDC", 
                "Unknown")




all.df <- readRDS("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/ThyPtsCombined.rds")
all.df$TissueType <- ifelse(grepl("Thy",all.df$SampleIDshort),'Thymus','Patient')


p1 <- ggplot(all.df,  aes(x = umap1, y = umap2, color = TissueType)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/ThyPtsCombined_TissueType.pdf", p1)
p2 <- ggplot(subset(all.df,TissueType %in% c("Thymus")),  aes(x = umap1, y = umap2, color = TissueType)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/ThyPtsCombined_TissueType_Thymus.pdf", p2)
p3 <- ggplot(subset(all.df,TissueType %in% c("Patient")),  aes(x = umap1, y = umap2, color = TissueType)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/ThyPtsCombined_TissueType_Patient.pdf", p3)
p4 <- ggplot(subset(all.df,TissueType %in% c("Thymus")),  aes(x = umap1, y = umap2, color = Celltype)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/ThyPtsCombined_TissueType_Thymus_CT.pdf", p4)
p5 <- ggplot(subset(all.df,TissueType %in% c("Patient")),  aes(x = umap1, y = umap2, color = Celltype)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + scale_color_manual(values = pal2) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave(filename = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/ThyPtsCombined_TissueType_Thymus_CT_pts.pdf", p5)



unique(all.df$SampleIDshort)

ggplot(subset(all.df,SampleIDshort %in% c("pt335")),  aes(x = umap1, y = umap2)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + scale_color_manual(values = pal) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggplot(subset(all.df,SampleIDshort %in% c("pt335")),  aes(x = umap1, y = umap2)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 

sample.list <- c(unique(all.df$SampleIDshort))
sample.list <- sample.list[! grepl("Thy",sample.list)]


for (s in sample.list){
  print(s)
  p <- ggplot(subset(all.df,SampleIDshort %in% c(s)),  aes(x = umap1, y = umap2)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + xlim(-21.438626, 9.383204)+ ylim(-5.733191, 9.602134)
  ggsave(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/PerPatientPlot/", s, "_umap.pdf"), p)
}


ggplot(subset(all.df,SampleIDshort %in% c(s)),  aes(x = umap1, y = umap2), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 
ggplot(subset(all.df,SampleIDshort %in% c("pt2723")),  aes(x = umap1, y = umap2)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() + xlim(-21.438626, 9.383204)+ ylim(-5.733191, 9.602134)

#[1] "pt5438"
#Saving 7.04 x 3.85 in image
#Warning messages:
#  1: Removed 1 rows containing missing values (geom_point). 
#  2: Removed 1 rows containing missing values (geom_point). 
p1$plot.margin
layer_scales(p1)$y$range$range
layer_scales(p1)$x$range$range





test.df <- subset(all.df,SampleIDshort %in% c(s))
#test.df$Kmeans <- kmeans(test.df[,c("umap1", "umap2")], centers = 10)
#knn(test.df[c("umap1", "umap2")])
#test.df$Knn <- knn.dist(test.df[c("umap1", "umap2")])
#ggplot(test.df,  aes(x = umap1, y = umap2), color = Knn)
#knn.tsne.trans = get.knn(as.matrix(raw.tsne.trans$Y), k = k)
test.knn <- get.knn(as.matrix(test.df[c("umap1", "umap2")]), k = 10)
test.knn = data.frame(from = rep(1:nrow(test.knn$nn.index), 
                                      k), to = as.vector(test.knn$nn.index), weight = 1/(1 + as.vector(test.knn$nn.dist)))
test.knn.nw = graph_from_data_frame(test.knn, directed = FALSE)
test.knn.nw = simplify(test.knn.nw)
test.knn.lc = cluster_louvain(test.knn.nw)
test.df$umap.louvain = as.factor(membership(test.knn.lc))

library(Seurat)
#test.srat <- FindClusters(asS4(test.df) )
#test.srat <- FindClusters(asS4(test.knn) )
#BiocManager::install("clusterExperiment")
library(clusterExperiment)
library(leiden)
#clusterExpiment()
#test.ce <- ClusterExperiment(as.data.frame(test.knn))
test.leiden <- leiden(test.knn.nw)
table(test.leiden)
table(test.df$umap.louvain)




ggplot(test.df,  aes(x = umap1, y = umap2), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 
ggplot(test.df,  aes(x = umap1, y = umap2, color = umap.louvain), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 



test.knn.1k <- get.knn(as.matrix(test.df[c("umap1", "umap2")]), k = 1000)
test.knn.1k = data.frame(from = rep(1:nrow(test.knn.1k$nn.index), 
                                 k), to = as.vector(test.knn.1k$nn.index), weight = 1/(1 + as.vector(test.knn.1k$nn.dist)))
test.knn.1k.nw = graph_from_data_frame(test.knn.1k, directed = FALSE)
test.knn.1k.nw = simplify(test.knn.1k.nw)
test.knn.1k.lc = cluster_louvain(test.knn.1k.nw)
test.df$umap.1k.louvain = as.factor(membership(test.knn.1k.lc))

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/testUmapClusters/louvain_1k.pdf")
ggplot(test.df,  aes(x = umap1, y = umap2, color = umap.1k.louvain), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 
dev.off()

test.knn.0.5k <- get.knn(as.matrix(test.df[c("umap1", "umap2")]), k = 500)
test.knn.0.5k = data.frame(from = rep(1:nrow(test.knn.0.5k$nn.index), 
                                      k), to = as.vector(test.knn.0.5k$nn.index), weight = 1/(1 + as.vector(test.knn.0.5k$nn.dist)))
test.knn.0.5k.nw = graph_from_data_frame(test.knn.0.5k, directed = FALSE)
test.knn.0.5k.nw = simplify(test.knn.0.5k.nw)
test.knn.0.5k.lc = cluster_louvain(test.knn.0.5k.nw)
test.df$umap.0.5k.louvain = as.factor(membership(test.knn.0.5k.lc))

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/testUmapClusters/louvain_0.5k.pdf")
ggplot(test.df,  aes(x = umap1, y = umap2, color = umap.0.5k.louvain), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 
dev.off()

test.knn.0.1k <- get.knn(as.matrix(test.df[c("umap1", "umap2")]), k = 100)
test.knn.0.1k = data.frame(from = rep(1:nrow(test.knn.0.1k$nn.index), 
                                      k), to = as.vector(test.knn.0.1k$nn.index), weight = 1/(1 + as.vector(test.knn.0.1k$nn.dist)))
test.knn.0.1k.nw = graph_from_data_frame(test.knn.0.1k, directed = FALSE)
test.knn.0.1k.nw = simplify(test.knn.0.1k.nw)
test.knn.0.1k.lc = cluster_louvain(test.knn.0.1k.nw)
test.df$umap.0.1k.louvain = as.factor(membership(test.knn.0.1k.lc))

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/testUmapClusters/Louvain_0.1k.pdf")
ggplot(test.df,  aes(x = umap1, y = umap2, color = umap.0.1k.louvain), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 
dev.off()



table(all.df$umap.louvain)
ggplot(all.df,  aes(x = umap1, y = umap2, color = umap.louvain), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 

BiocManager::install("scDataviz")


umap.input <- test.df[,1:17] %>% select(- contains("CD45"))
library(scDataviz)
test.datavix <- clusKNN(umap.input,
                        k.param = 20,
                        prune.SNN = 1/15,
                        resolution = 0.01,
                        algorithm = 2)
test.datavix <- clusKNN(umap.input)
test.df$Dataviz <- as.factor(test.datavix)

#all.df$Dataviz <- clusKNN(all.df[1:17])
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/testUmapClusters/Dataviz_exprs.pdf")
ggplot(test.df,  aes(x = umap1, y = umap2, color = Dataviz), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 
dev.off()
#ggplot(all.df,  aes(x = umap1, y = umap2, color = Dataviz), xlim = c(-21.438626, 9.383204),
#      ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.1, show.legend = TRUE) + theme_bw() 

test.datavix.umap <- clusKNN(test.df[c("umap1", "umap2")],
                             k.param = 20,
                             prune.SNN = 1/15,
                             resolution = 0.01,
                             algorithm = 2)
test.df$Dataviz.umap <- as.factor(test.datavix.umap)


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/testUmapClusters/Dataviz_umap.pdf")
ggplot(test.df,  aes(x = umap1, y = umap2, color = Dataviz.umap), xlim = c(-21.438626, 9.383204),
       ylim = c(-5.733191, 9.602134)) + geom_point(size = 0.2, show.legend = TRUE) + theme_bw() 
dev.off()


library(M3C)
m3c.data <- M3C(umap.input)
r <- M3C(t(umap.input),method=2)
umap(pollen$data,labels=as.factor(r$assignments),printres = TRUE,printwidth = 24)

m3c.data$assignments
r$assignments


