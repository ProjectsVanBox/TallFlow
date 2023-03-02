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



setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow")
### Loop over file and select columns, add sample name 
file.list <-  list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/PatientStrict/", pattern = "*.csv$", full.names = TRUE, recursive = TRUE)
combine.list <- list()
for (f in file.list){
  data.file <- read.table(file = f, sep = ",", header = TRUE, row.names = NULL)
  sample.name <- gsub(x = gsub(x = f, pattern = "/.*/", replacement = ""), pattern = "_keymarkers.csv", replacement = "")
  print(sample.name)
  data.file$SampleID <- sample.name
  data.sub <- data.file[c("Index", "Celltype", "SampleID")]
  combine.list[[sample.name]] <- data.sub
}
combine.file <- do.call("rbind", combine.list)




############################################################################# Add cell count filter

dim(combine.file)
combine.file <- combine.file[! combine.file$SampleID %in% c("C2_pt5716","C3_pt5716"), ]
#unique(combine.file$SampleID)


combine.file.sub <- combine.file[! combine.file$Celltype %in% c("cDC-like", "B-cell-like","Monocyte-like", "NK-cell-like", "pDC-like", "Unknown"), ]
combine.file.sub$MergedSampleName <-  gsub(".*_", "", combine.file.sub$SampleID)

### Count the occurence of all labels per file
LabelCounts <- aggregate(groupcount ~ Celltype + SampleID, transform(combine.file.sub, groupcount = 1), length)
LabelCountsTable <- dcast(data = LabelCounts, formula = Celltype ~ SampleID)
write.table(LabelCountsTable, file = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/LabelCountsTable_Patients_TPop.txt", sep = "\t", row.names = FALSE)

### Calculate percentages
PercentagesLabel <- LabelCountsTable[1:ncol(LabelCountsTable)]
rownames(PercentagesLabel) <- LabelCountsTable$Celltype
PercentagesLabel<-apply(PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(PercentagesLabel, file = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/PercentagesLabel_Patients_TPop.txt", sep = "\t", row.names = T, col.names = NA)

### Reorder the labels
PercentagesPerLabel <- melt(as.matrix(PercentagesLabel))
colnames(PercentagesPerLabel) <- c("Celltype", "SampleID", "Value")
PercentagesPerLabel$Celltype <- factor(PercentagesPerLabel$Celltype, 
                                       levels = rev(c("DN-like", "ISP-like", "gdTcell-like","DP-like", "CD4-like", "CD8-like")))
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00")
names(pal) <- c("DN-like", "ISP-like", "gdTcell-like","DP-like", "CD4-like", "CD8-like")

### Convert NA to 0 
PercentagesLabel <- as.data.frame(PercentagesLabel)
PercentagesLabel[is.na(PercentagesLabel)] <- 0


### Create dendrograms
dend <- as.dendrogram(hclust(as.dist(1-cor(PercentagesLabel))))
dend <- set(dend, "labels_cex", 0.7) # gives error 
#PercentagesPerLabel$SampleID <- factor(PercentagesPerLabel$SampleID, levels = labels(dend))
PercentagesPerLabel$MergedSampleName <-  gsub(".*_", "", PercentagesPerLabel$SampleID)


#ordered by EGIL subtype
PercentagesPerLabel$MergedSampleName <- factor(PercentagesPerLabel$MergedSampleName, 
                                               levels = c("pt8173", "pt540", "pt2322", "pt913", "pt335", "pt419", "pt2229", "pt1946", 
                                                          "pt633", "pt1179", "pt5676", "pt7675", "pt8711", "pt3976", "pt1524", "pt2789", 
                                                          "pt2723", "pt5438", "pt14643","pt7267","pt1949","pt1950","pt258","pt530",
                                                          "pt9175","ptvu10138","pt4564", "pt8148"))

df2 <- PercentagesPerLabel[order(PercentagesPerLabel$MergedSampleName),]
tail(df2)
df2$SampleID <- factor(df2$SampleID, levels = rev(unique(df2$SampleID)))
p.bar <- ggplot(df2, aes(x = SampleID, y = Value, fill=Celltype)) +
  geom_bar(stat='identity')+ coord_flip() + theme_classic() + scale_fill_manual(values =pal)
p.bar
ggsave(filename = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/EgilOrderBargraph_Patients_TPop.pdf", p.bar, width = 7, height = 9)

#ordered by oncogene
PercentagesPerLabel$MergedSampleName <- factor(PercentagesPerLabel$MergedSampleName, 
                                               levels = c("pt8173","pt540","pt913","pt1524","pt8639","pt1946","pt633","pt1179","pt5676",
                                                          "pt3976","pt2723","pt419","pt2229","ptvu10138","pt2322","pt335","pt2789","pt1949",
                                                          "pt1950","pt258","pt530","pt9175","pt7267","pt7675","pt8711","pt5438","pt14643",
                                                          "pt4564","pt8148"))
PercentagesPerLabel$MergedSampleName
df2 <- PercentagesPerLabel[order(PercentagesPerLabel$MergedSampleName),]
tail(df2)
df2$SampleID <- factor(df2$SampleID, levels = rev(unique(df2$SampleID)))
p.bar <- ggplot(df2, aes(x = SampleID, y = Value, fill=Celltype)) +
  geom_bar(stat='identity')+ coord_flip() + theme_classic() + scale_fill_manual(values =pal)
p.bar
ggsave(filename = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/OncogeneOrderBargraph_Patients_TPop.pdf", p.bar, width = 7, height = 9)



### Older code: 

### Create plots 
p.den <- ggplot(as.ggdend(dend), horiz = TRUE, theme = theme_classic()) 
p.bar <- ggplot(PercentagesPerLabel, aes(x = SampleID, y = Value, fill=Celltype)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/Dendrogram_Patients_TPop.pdf", p.den, width = 7, height = 9)
ggsave(filename = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/OrderedBargraph_Patients_TPop.pdf", p.bar, width = 7, height = 9)






### Get clustering and plot it 
#https://www.datacamp.com/tutorial/hierarchical-clustering-R
avg_col_dend <- color_branches(dend, k = 8)
#plot(avg_col_dend)
p.den.col <- ggplot(as.ggdend(avg_col_dend), horiz = TRUE, theme = theme_classic()) 
ggsave(filename = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/DendrogramCol_Patients_TPop.pdf", p.den.col, width = 7, height = 9)


avg_col_dend8 <- color_branches(dend, k = 8)
p.den8.col <- ggplot(as.ggdend(avg_col_dend8), horiz = TRUE, theme = theme_classic()) 


avg_col_dend9 <- color_branches(dend, k = 9)
p.den9.col <- ggplot(as.ggdend(avg_col_dend9), horiz = TRUE, theme = theme_classic()) 

avg_col_dend10 <- color_branches(dend, k = 10)
p.den10.col <- ggplot(as.ggdend(avg_col_dend10), horiz = TRUE, theme = theme_classic()) 

avg_col_dend11 <- color_branches(dend, k = 11)
p.den11.col <- ggplot(as.ggdend(avg_col_dend11), horiz = TRUE, theme = theme_classic()) 
avg_col_dend12 <- color_branches(dend, k = 12)
p.den12.col <- ggplot(as.ggdend(avg_col_dend12), horiz = TRUE, theme = theme_classic()) 

p.den8.col
p.den9.col
p.den10.col
p.den11.col
p.den12.col

ggsave(filename = "./3_Output/Barcharts/PatientStrict/OnlyT_Populations/DendrogramCol_Patients_TPop_k10.pdf", p.den10.col, width = 7, height = 9)



