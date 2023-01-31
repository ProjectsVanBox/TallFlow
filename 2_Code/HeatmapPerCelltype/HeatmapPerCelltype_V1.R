#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(data.table)
library(dendextend)
library(RColorBrewer)
library(pheatmap)


### Function to create a heatmap per cluster
plot_clustering_heatmap_wrapper <- function(expr, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% 
    summarize_all(list(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% 
    summarize_all(list(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (", 
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster 
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
# Colors for the heatmap
color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  pheatmap(expr_heat, color = color, 
           cluster_cols = TRUE, cluster_rows = cluster_rows, 
           labels_col = labels_col, labels_row = labels_row, 
           display_numbers = TRUE, number_color = "black", 
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors, 
           annotation_legend = annotation_legend)
  
}


### Function to get distinct colours 
getDistinctColors <- function(n) {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unique (unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))));
  stopifnot (n <= length(col_vector));
  xxx <- col2rgb(col_vector);
  dist_mat <- as.matrix(dist(t(xxx)));
  diag(dist_mat) <- 1e10;
  while (length(col_vector) > n) {
    minv <- apply (dist_mat,1,function(x)min(x));
    idx <- which(minv==min(minv))[1];
    dist_mat <- dist_mat[-idx, -idx];
    col_vector <- col_vector[-idx]
  }
  col_vector
}


### Z-score 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}





setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow")
### Loop over file and select columns, add sample name 
file.list <-  list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Patients/", pattern = "*.csv$", full.names = TRUE, recursive = TRUE)
combine.list <- list()
for (f in file.list){
  data.file <- read.table(file = f, sep = ",", header = TRUE, row.names = NULL)
  sample.name <- gsub(x = gsub(x = f, pattern = "/.*/", replacement = ""), pattern = "_keymarkers.csv", replacement = "")
  SampleNameShort <- gsub("_.*", "", gsub(".*pt", "pt", sample.name))
  data.file$SampleID <- sample.name
  data.file$SampleIDshort <- SampleNameShort
  print(sample.name)
  
  ### Select data 
  df <- data.file %>% select(matches('__'))
  colnames(df) <- gsub(x = gsub(colnames(df), pattern = ".*__", replacement = ""), pattern = ".A", replacement = "")
  df <- df %>% select(-starts_with("Celltype")) %>% select(-starts_with("DeadCells")) %>% 
    select(- contains("invNK_T")) %>% select(- contains("DeadCells")) %>% select(- contains("CD127")) #%>% select(- contains("CD45"))
  df$Celltype <- data.file$Celltype
  df$SampleID <- data.file$SampleID
  df$SampleIDshort <- data.file$SampleIDshort
  #data.sub <- data.file[c("Index", "Celltype", "SampleID")]
  combine.list[[sample.name]] <- df
}
### Combine and filter data 
combine.file <- do.call("rbind", combine.list)
dim(combine.file)
combine.file<- combine.file[combine.file$SampleIDshort != "pt5716",]
combine.file.sub <- combine.file[! combine.file$Celltype %in% c("cDC-like", "B-cell-like","Monocyte-like", "NK-cell-like", "pDC-like", "Unknown"), ]


### Calculate the average expression of each replicate 
combine.average <- combine.file.sub %>% group_by(Celltype , SampleID) %>% summarise_each(funs(mean=mean(., na.rm=TRUE)))
input.heat <- as.matrix(combine.average[3:19])
rownames(input.heat) <- paste0(combine.average$Celltype, "_", combine.average$SampleID)




### Perform the z-score part and see if it works better/worse
pheatmap(input.heat)
input.heat.z <- t(apply(input.heat, 1, cal_z_score))
pheatmap(input.heat.z)




### Get the annotation and colours 
#my_row_col <- read.table("./1_Input/PatientCellTypeClustering_Strict.txt", sep = "\t", header = T, row.names = 1)
my_row_col <- read.table("./1_Input/PatientCellTypeClustering_Strict_k10.txt", sep = "\t", header = T, row.names = 1)
my_row_col <- my_row_col[c(1,4)]
cell.color <- getDistinctColors(length(unique(my_row_col$CellType)))
names(cell.color) <- unique(my_row_col$CellType)
#sample.color <- getDistinctColors(length(unique(my_row_col$SampleIDshort)))
#names(sample.color) <- unique(my_row_col$SampleIDshort)
Cluster.color <- getDistinctColors(length(unique(my_row_col$Cluster)))
names(Cluster.color) <- unique(my_row_col$Cluster)
# List the colours
my_colour = list(
  CellType = cell.color,
  #SampleIDshort = sample.color,
  Cluster = Cluster.color
)

pdf("./3_Output/HeatmapPerCelltype/HeatmapReplicateSamples_TPop_zscore.pdf")
pheatmap(input.heat.z, annotation_colors = my_colour, annotation_row = my_row_col, show_rownames = F)
dev.off()



my_row_col.order <- my_row_col[order( my_row_col[,"CellType"], my_row_col[,"Cluster"] ),]
input.heat.z.order <- input.heat.z[rownames(my_row_col.order),]
input.heat.order <- input.heat[rownames(my_row_col.order),]

pdf("./3_Output/HeatmapPerCelltype/HeatmapReplicateSamples_TPop_zscore_NoRowClust.pdf")
pheatmap(input.heat.z.order, annotation_colors = my_colour, annotation_row = my_row_col, show_rownames = F, cluster_rows = F)
dev.off()



celltype.heat <- celltype.heat[marker.order]
celltype.heat.z <- celltype.heat.z[marker.order]
### Get a part per Cell type 
for (ct in unique(my_row_col.order$CellType)){
  print(ct)
  # Subselect from zscore 
  celltype.heat <- input.heat.z.order[rownames(input.heat.z.order) %in% rownames(my_row_col.order[my_row_col.order$CellType == ct, ]), ]
  pdf(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/HeatmapPerCelltype/subselect_zscore/heatmap_", ct, "_NoClust.pdf"))
  pheatmap(celltype.heat, annotation_colors = my_colour, annotation_row = my_row_col, cluster_rows = F, cluster_cols = F)
  dev.off()
  pdf(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/HeatmapPerCelltype/subselect_zscore/heatmap_", ct, ".pdf"))
  pheatmap(celltype.heat, annotation_colors = my_colour, annotation_row = my_row_col, cluster_cols = F)
  dev.off()
  # Recalculate zscore
  celltype.heat.z <- input.heat.order[rownames(input.heat.order) %in% rownames(my_row_col.order[my_row_col.order$CellType == ct, ]), ]
  celltype.heat.z <- t(apply(celltype.heat.z, 1, cal_z_score))
  pdf(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/HeatmapPerCelltype/ReCalculate_zscore/heatmap_", ct, "_NoClust.pdf"))
  pheatmap(celltype.heat.z, annotation_colors = my_colour, annotation_row = my_row_col, cluster_rows = F, cluster_cols = F)
  dev.off()
  pdf(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/HeatmapPerCelltype/ReCalculate_zscore/heatmap_", ct, ".pdf"))
  pheatmap(celltype.heat.z, annotation_colors = my_colour, annotation_row = my_row_col, cluster_cols = F)
  dev.off()
  # Normalised value 
  celltype.norm.heat <- norm.all.df.order[rownames(norm.all.df.order) %in% rownames(my_row_col.order[my_row_col.order$CellType == ct, ]), ]
  pdf(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/HeatmapPerCelltype/Normalised/heatmap_", ct, "_NoClust.pdf"))
  pheatmap(celltype.norm.heat, cluster_cols = F, cluster_rows = F, annotation_row = my_row_col, fontsize_row = 4, annotation_colors = my_colour)
  dev.off()
  pdf(paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/HeatmapPerCelltype/Normalised/heatmap_", ct, ".pdf"))
  pheatmap(celltype.norm.heat, cluster_cols = F, annotation_row = my_row_col, fontsize_row = 4, annotation_colors = my_colour)
  dev.off()
}


celltype.heat <- input.heat.z.order[rownames(input.heat.z.order) %in% rownames(my_row_col.order[my_row_col.order$CellType == "DP-like", ]), ]
pheatmap(celltype.heat, annotation_colors = my_colour, annotation_row = my_row_col, cluster_rows = F)

celltype.heat.z <- input.heat.order[rownames(input.heat.order) %in% rownames(my_row_col.order[my_row_col.order$CellType == "DP-like", ]), ]
celltype.heat.z <- t(apply(celltype.heat.z, 1, cal_z_score))
pheatmap(celltype.heat.z, annotation_colors = my_colour, annotation_row = my_row_col, cluster_rows = F)
pheatmap(celltype.heat.z, annotation_colors = my_colour, annotation_row = my_row_col)

celltype.norm.heat <- norm.all.df.order[rownames(norm.all.df.order) %in% rownames(my_row_col.order[my_row_col.order$CellType == "DP-like", ]), ]
pheatmap(celltype.norm.heat, cluster_cols = F, cluster_rows = F, annotation_row = my_row_col, fontsize_row = 4, annotation_colors = my_colour)
pheatmap(celltype.norm.heat, cluster_cols = F, annotation_row = my_row_col, fontsize_row = 4, annotation_colors = my_colour)




### Put all values between 0 and 1 
# find min and max of all columns in total dataset
options(scipen=999)
min.v <- as.matrix(apply(input.heat,2,min))
min.v <- as.data.frame(min.v)
min.v <- as.numeric(gsub(",","",min.v$V1))
max.v <- as.matrix(apply(input.heat,2,max))
max.v <- as.data.frame(max.v)
max.v <- as.numeric(gsub(",","",max.v$V1))
### nomalize dataframe containing all files to min and max of each cell type 
norm.all.df <- t((t(input.heat) - min.v) / (max.v - min.v))
norm.all.df[norm.all.df < 0] <- 0
norm.all.df[norm.all.df > 1] <- 1
norm.all.df <- as.data.frame(norm.all.df)


### Change names and order 
marker.order <- c("CD7_mean", "CD5_mean", "CD117_mean", "CD45_mean", "CD44_mean", "CD1a_mean", "CD3_mean", "CD4_mean", "CD8a_mean", 
                  "TCRa_b_mean", "TCRg_d_mean", "CD25_mean", "CD56_mean", "CD123_mean", "CD11c_mean", "CD14_mean", "CD19_mean" )
colnames(norm.all.df)
norm.all.df <- norm.all.df[marker.order]

#colnames(norm.all.df) <- gsub("_mean", "", colnames(norm.all.df))




### Get the annotation and colours 
#my_row_col <- read.table("./1_Input/PatientCellTypeClustering_Strict.txt", sep = "\t", header = T, row.names = 1)
my_row_col <- read.table("./1_Input/PatientCellTypeClustering_Strict_k10.txt", sep = "\t", header = T, row.names = 1)
#my_row_col <- data.frame(Celltype = gsub("_.*", "", rownames(norm.all.df)))
#rownames(my_row_col) <- rownames(norm.all.df)


my_row_col <- my_row_col[c(1,3,4)]

cell.color <- getDistinctColors(length(unique(my_row_col$CellType)))
names(cell.color) <- unique(my_row_col$CellType)
sample.color <- getDistinctColors(length(unique(my_row_col$SampleIDshort)))
names(sample.color) <- unique(my_row_col$SampleIDshort)
Cluster.color <- getDistinctColors(length(unique(my_row_col$Cluster)))
names(Cluster.color) <- unique(my_row_col$Cluster)

norm.all.df.order <- norm.all.df[rownames(my_row_col.order),]

my_colour = list(
  CellType = cell.color,
  SampleIDshort = sample.color,
  Cluster = Cluster.color
)

pdf("./3_Output/HeatmapPerCelltype/HeatmapReplicateSamples_TPop_NoColClust.pdf")
pheatmap(norm.all.df, cluster_cols = F, annotation_row = my_row_col, fontsize_row = 4, annotation_colors = my_colour)
dev.off()





my_row_col
my_row_col[
  with(my_row_col, order("CellType", "Cluster")),
]


my_row_col.order <- my_row_col[order( my_row_col[,"CellType"], my_row_col[,"Cluster"] ),]


norm.all.df.order <- norm.all.df[rownames(my_row_col.order),]
pdf("./3_Output/HeatmapPerCelltype/HeatmapReplicateSamples_TPop_CustomOrder.pdf")
pheatmap(norm.all.df.order, cluster_cols = F, cluster_rows = F, annotation_row = my_row_col, fontsize_row = 4, annotation_colors = my_colour)
dev.off()



##################################################### Old code will follow






















### Calculate the average expression of each sample 
combine.average <- combine.file.sub %>% group_by(Celltype , SampleIDshort) %>% summarise_each(funs(mean=mean(., na.rm=TRUE)))
input.heat <- as.matrix(combine.average[3:19])
rownames(input.heat) <- paste0(combine.average$Celltype, "_", combine.average$SampleIDshort)



# find min and max of all columns in total dataset
options(scipen=999)
min.v <- as.matrix(apply(input.heat,2,min))
min.v <- as.data.frame(min.v)
min.v <- as.numeric(gsub(",","",min.v$V1))
max.v <- as.matrix(apply(input.heat,2,max))
max.v <- as.data.frame(max.v)
max.v <- as.numeric(gsub(",","",max.v$V1))
### nomalize dataframe containing all files to min and max of each cell type 
norm.all.df <- t((t(input.heat) - min.v) / (max.v - min.v))
norm.all.df[norm.all.df < 0] <- 0
norm.all.df[norm.all.df > 1] <- 1
norm.all.df <- as.data.frame(norm.all.df)


marker.order <- c("CD7_mean", "CD5_mean", "CD117_mean", "CD45_mean", "CD44_mean", "CD1a_mean", "CD3_mean", "CD4_mean", "CD8a_mean", 
                  "TCRa_b_mean", "TCRg_d_mean", "CD25_mean", "CD56_mean", "CD123_mean", "CD11c_mean", "CD14_mean", "CD19_mean" )

colnames(norm.all.df)
norm.all.df <- norm.all.df[marker.order]

colnames(norm.all.df) <- gsub("_mean", "", colnames(norm.all.df))


my_row_col <- read.table("./1_Input/PatientCellTypeClustering.txt", sep = "\t", header = T, row.names = 1)


#my_row_col <- data.frame(Celltype = gsub("_.*", "", rownames(norm.all.df)))
#rownames(my_row_col) <- rownames(norm.all.df)



pdf("./3_Output/HeatmapPerCelltype/HeatmapMergedSamples_TPop_NoClust.pdf")
pheatmap(norm.all.df, cluster_cols = F, cluster_rows = F, annotation_row = my_row_col, fontsize_row = 4)
dev.off()
pdf("./3_Output/HeatmapPerCelltype/HeatmapMergedSamples_TPop_NoRowClust.pdf")
pheatmap(norm.all.df, cluster_cols = F, annotation_row = my_row_col, fontsize_row = 4)
dev.off()
pdf("./3_Output/HeatmapPerCelltype/HeatmapMergedSamples_TPop_AllClust.pdf")
pheatmap(norm.all.df, annotation_row = my_row_col, fontsize_row = 4)
dev.off()



for (rn in rownames(my_row_col)){
  print(rn)
}

dev.off()















########################## Old code 


#test.average <- combine.file %>% group_by(Celltype , SampleIDshort) %>% summarise_each(funs(mean=mean(., na.rm=TRUE)))
test.average <- combine.file %>% group_by(Celltype , SampleID) %>% summarise_each(funs(mean=mean(., na.rm=TRUE)))
#test.average <- combine.file %>% group_by(SampleIDshort) %>% summarise_each(funs(mean=mean(., na.rm=TRUE)))


input.heat <- as.matrix(test.average[3:19])
rownames(input.heat) <- paste0(test.average$Celltype, test.average$SampleID)

pheatmap(t(input.heat))


test.average <- combine.file %>% group_by(Celltype , SampleID) %>% summarise_each(funs(mean=mean(., na.rm=TRUE)))
input.heat.dn <- test.average[test.average$Celltype == "DN-like", ]
#input.heat.dn <- as.matrix(input.heat.dn)
rownames.heat  <- paste0(input.heat.dn$Celltype, input.heat.dn$SampleID)
input.heat.dn <- as.matrix(input.heat.dn[3:19])
rownames(input.heat.dn) <- rownames.heat
#pheatmap(t(input.heat.dn))


# find min and max of all columns in total dataset
options(scipen=999)
min.v <- as.matrix(apply(input.heat.dn,2,min))
min.v <- as.data.frame(min.v)
min.v <- as.numeric(gsub(",","",min.v$V1))
max.v <- as.matrix(apply(input.heat.dn,2,max))
max.v <- as.data.frame(max.v)
max.v <- as.numeric(gsub(",","",max.v$V1))
### nomalize dataframe containing all files to min and max of healthy samples. 
#exprs.trim <- trans.all.df %>% dplyr::select(- contains("File")) %>% dplyr::select(- contains("PatientID"))
norm.all.df <- t((t(input.heat.dn) - min.v) / (max.v - min.v))
norm.all.df[norm.all.df < 0] <- 0
norm.all.df[norm.all.df > 1] <- 1
norm.all.df <- as.data.frame(norm.all.df)
#norm.all.df <- cbind(norm.all.df,trans.all.df[,c("File", "PatientID")])
#names(norm.all.df)[18:19] <- c("File", "PatientID")


pdf("./3_Output/HeatmapPerCelltype/HeatmapTest.pdf")
plot_clustering_heatmap_wrapper(expr = norm.all.df, 
                                expr01 = norm.all.df, 
                                cell_clustering = rownames(norm.all.df) , color_clusters = getDistinctColors(70))
dev.off()



# find min and max of all columns in total dataset
options(scipen=999)
min.v <- as.matrix(apply(input.heat,2,min))
min.v <- as.data.frame(min.v)
min.v <- as.numeric(gsub(",","",min.v$V1))
max.v <- as.matrix(apply(input.heat,2,max))
max.v <- as.data.frame(max.v)
max.v <- as.numeric(gsub(",","",max.v$V1))
### nomalize dataframe containing all files to min and max of healthy samples. 
#exprs.trim <- trans.all.df %>% dplyr::select(- contains("File")) %>% dplyr::select(- contains("PatientID"))
norm.all.df <- t((t(input.heat) - min.v) / (max.v - min.v))
norm.all.df[norm.all.df < 0] <- 0
norm.all.df[norm.all.df > 1] <- 1
norm.all.df <- as.data.frame(norm.all.df)


pdf("./3_Output/HeatmapPerCelltype/HeatmapTest_AllCelltypes.pdf")
plot_clustering_heatmap_wrapper(expr = norm.all.df, 
                                expr01 = norm.all.df, 
                                cell_clustering = rownames(norm.all.df) , color_clusters = getDistinctColors(70))
dev.off()
#test.med <- combine.file %>% group_by(Celltype , SampleIDshort) %>% summarise_each(funs(median(., na.rm=TRUE)))


pheatmap(norm.all.df)
