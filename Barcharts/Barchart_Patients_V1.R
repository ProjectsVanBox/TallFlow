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
file.list <-  list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Patients/", pattern = "*.csv$", full.names = TRUE, recursive = TRUE)
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


### Count the occurence of all labels per file
LabelCounts <- aggregate(groupcount ~ Celltype + SampleID, transform(combine.file, groupcount = 1), length)
LabelCountsTable <- dcast(data = LabelCounts, formula = Celltype ~ SampleID)
write.table(LabelCountsTable, file = "./3_Output/Barcharts/Patients/LabelCountsTable_Patients.txt", sep = "\t", row.names = FALSE)

### Calculate percentages
PercentagesLabel <- LabelCountsTable[1:ncol(LabelCountsTable)]
rownames(PercentagesLabel) <- LabelCountsTable$Celltype
PercentagesLabel<-apply(PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(PercentagesLabel, file = "./3_Output/Barcharts/Patients/PercentagesLabel_Patients.txt", sep = "\t", row.names = T, col.names = NA)

### Reorder the labels
PercentagesPerLabel <- melt(as.matrix(PercentagesLabel))
colnames(PercentagesPerLabel) <- c("Celltype", "SampleID", "Value")
PercentagesPerLabel$Celltype <- factor(PercentagesPerLabel$Celltype, 
                                       levels = rev(c("DN-like", "ISP-like", "gdTcell-like","DP-like", "CD4-like", "CD8-like", "cDC-like", 
                                                      "B-cell-like","Monocyte-like", "NK-cell-like", "pDC-like", 
                                                      "Unknown")))
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN-like", "ISP-like", "gdTcell-like","DP-like", "CD4-like", "CD8-like", "cDC-like", 
                "B-cell-like","Monocyte-like", "NK-cell-like", "pDC-like", 
                "Unknown")

### Convert NA to 0 
PercentagesLabel <- as.data.frame(PercentagesLabel)
PercentagesLabel[is.na(PercentagesLabel)] <- 0


### Create dendrograms
dend <- as.dendrogram(hclust(as.dist(1-cor(PercentagesLabel))))
dend <- set(dend, "labels_cex", 0.7) # gives error 
PercentagesPerLabel$SampleID <- factor(PercentagesPerLabel$SampleID, levels = labels(dend))


### Create plots 
p.den <- ggplot(as.ggdend(dend), horiz = TRUE, theme = theme_classic()) 
p.bar <- ggplot(PercentagesPerLabel, aes(x = SampleID, y = Value, fill=Celltype)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/Patients/Dendrogram_Patients.pdf", p.den)
ggsave(filename = "./3_Output/Barcharts/Patients/OrderedBargraph_Patients.pdf", p.bar)








