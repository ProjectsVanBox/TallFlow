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
file.list <-  list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/ThymusManGated/", pattern = "*.csv$", full.names = TRUE, recursive = TRUE)
combine.list <- list()
for (f in file.list){
  data.file <- read.table(file = f, sep = ",", header = TRUE, row.names = NULL)
  sample.name <- gsub(x = gsub(x = f, pattern = "/.*/", replacement = ""), pattern = "_keymarkers.csv", replacement = "")
  print(sample.name)
  data.file$SampleID <- sample.name
  data.sub <- data.file[c("Index", "Celltype", "SampleID", "ManualGate")]
  combine.list[[sample.name]] <- data.sub
}
combine.file <- do.call("rbind", combine.list)


### Seperate un- and depleted samples 
undepleted_samples <- c("A2_ptThy072_undepleted", "A2_ptThy073_undepleted", "A2_ptThy075_undepleted", "A2_ptThy076_undepleted", 
                        "A2_ptThy079_undepleted", "A2_ptThy083_undepleted", "A2_ptThy085_Undepleted", "A6_ptThy070_undepleted",
                        "A6_ptThy071_thawed", "A10_ptThy075_thawed", "A12_ptThy076_thawed", "A2_ptThy71")
depleted_samples <- c("A4_ptThy073_70Mdepleted", "A4_ptThy075_70Mdepleted", "A4_ptThy076_70Mdepleted", "A4_ptThy079_depleted", 
                      "A4_ptThy083_depleted", "A4_ptThy085_Depleted", "A4_ptThy71_depleted", "A6_ptThy072_70Mdepleted", 
                      "A6_ptThy073_10Mdepleted", "A6_ptThy075_10Mdepleted", "A8_ptThy070_depleted")
undepleted_file <- combine.file[combine.file$SampleID %in% undepleted_samples, ]
depleted_file <- combine.file[combine.file$SampleID %in% depleted_samples, ]


### Count the occurence of all labels per file
undepleted_LabelCounts <- aggregate(groupcount ~ Celltype + SampleID, transform(undepleted_file, groupcount = 1), length)
undepleted_LabelCountsTable <- dcast(data = undepleted_LabelCounts, formula = Celltype ~ SampleID)
write.table(undepleted_LabelCountsTable, file = "./3_Output/Barcharts/Thymus_pygated/undepleted_LabelCountsTable_Thymus.txt", sep = "\t", row.names = FALSE)
depleted_LabelCounts <- aggregate(groupcount ~ Celltype + SampleID, transform(depleted_file, groupcount = 1), length)
depleted_LabelCountsTable <- dcast(data = depleted_LabelCounts, formula = Celltype ~ SampleID)
write.table(depleted_LabelCountsTable, file = "./3_Output/Barcharts/Thymus_pygated/depleted_LabelCountsTable_Thymus.txt", sep = "\t", row.names = FALSE)


### Calculate percentages
undepleted_PercentagesLabel <- undepleted_LabelCountsTable[1:ncol(undepleted_LabelCountsTable)]
rownames(undepleted_PercentagesLabel) <- undepleted_LabelCountsTable$Celltype
undepleted_PercentagesLabel<-apply(undepleted_PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(undepleted_PercentagesLabel, file = "./3_Output/Barcharts/Thymus_pygated/undepleted_PercentagesLabel_Thymus.txt", sep = "\t", row.names = T, col.names = NA)
depleted_PercentagesLabel <- depleted_LabelCountsTable[1:ncol(depleted_LabelCountsTable)]
rownames(depleted_PercentagesLabel) <- depleted_LabelCountsTable$Celltype
depleted_PercentagesLabel<-apply(depleted_PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(depleted_PercentagesLabel, file = "./3_Output/Barcharts/Thymus_pygated/depleted_PercentagesLabel_Thymus.txt", sep = "\t", row.names = T, col.names = NA)


### Reorder the labels
undepleted_PercentagesPerLabel <- melt(as.matrix(undepleted_PercentagesLabel))
colnames(undepleted_PercentagesPerLabel) <- c("Celltype", "SampleID", "Value")
undepleted_PercentagesPerLabel$Celltype <- factor(undepleted_PercentagesPerLabel$Celltype, 
                                       levels = rev(c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                      "B-cell","Monocyte", "NK-cell", "pDC", 
                                                      "Unknown")))
depleted_PercentagesPerLabel <- melt(as.matrix(depleted_PercentagesLabel))
colnames(depleted_PercentagesPerLabel) <- c("Celltype", "SampleID", "Value")
depleted_PercentagesPerLabel$Celltype <- factor(depleted_PercentagesPerLabel$Celltype, 
                                                levels = rev(c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                               "B-cell","Monocyte", "NK-cell", "pDC", 
                                                               "Unknown")))


### Colour variables 
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                "B-cell","Monocyte", "NK-cell", "pDC", 
                "Unknown")


### Convert NA to 0 
undepleted_PercentagesLabel <- as.data.frame(undepleted_PercentagesLabel)
undepleted_PercentagesLabel[is.na(undepleted_PercentagesLabel)] <- 0
depleted_PercentagesLabel <- as.data.frame(depleted_PercentagesLabel)
depleted_PercentagesLabel[is.na(depleted_PercentagesLabel)] <- 0


### Create dendrograms
undepleted_dend <- as.dendrogram(hclust(as.dist(1-cor(undepleted_PercentagesLabel))))
undepleted_dend <- set(undepleted_dend, "labels_cex", 0.7) 
undepleted_PercentagesPerLabel$SampleID <- factor(undepleted_PercentagesPerLabel$SampleID, levels = labels(undepleted_dend))
depleted_dend <- as.dendrogram(hclust(as.dist(1-cor(depleted_PercentagesLabel))))
depleted_dend <- set(depleted_dend, "labels_cex", 0.7) 
depleted_PercentagesPerLabel$SampleID <- factor(depleted_PercentagesPerLabel$SampleID, levels = labels(depleted_dend))


### Create plots 
undepleted_p.den <- ggplot(as.ggdend(undepleted_dend), horiz = TRUE, theme = theme_classic()) 
undepleted_p.bar <- ggplot(undepleted_PercentagesPerLabel, aes(x = SampleID, y = Value, fill=Celltype)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/Thymus_pygated/undepleted_Dendrogram_Thymus.pdf", undepleted_p.den, width = 7, height = 5)
ggsave(filename = "./3_Output/Barcharts/Thymus_pygated/undepleted_OrderedBargraph_Thymus.pdf", undepleted_p.bar, width = 7, height = 5)
depleted_p.den <- ggplot(as.ggdend(depleted_dend), horiz = TRUE, theme = theme_classic()) 
depleted_p.bar <- ggplot(depleted_PercentagesPerLabel, aes(x = SampleID, y = Value, fill=Celltype)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/Thymus_pygated/depleted_Dendrogram_Thymus.pdf", depleted_p.den, width = 7, height = 5)
ggsave(filename = "./3_Output/Barcharts/Thymus_pygated/depleted_OrderedBargraph_Thymus.pdf", depleted_p.bar, width = 7, height = 5)


### Save enviroment
saveRDS(object = undepleted_PercentagesLabel, file = "./3_Output/Barcharts/Thymus_pygated/undepleted_PercentagesLabel.RDS")
saveRDS(object = undepleted_PercentagesPerLabel, file = "./3_Output/Barcharts/Thymus_pygated/undepleted_PercentagesPerLabel.RDS")
saveRDS(object = depleted_PercentagesLabel, file = "./3_Output/Barcharts/Thymus_pygated/depleted_PercentagesLabel.RDS")
saveRDS(object = depleted_PercentagesPerLabel, file = "./3_Output/Barcharts/Thymus_pygated/depleted_PercentagesPerLabel.RDS")






################################################################## Manually gated ##################################################################
### Count the occurence of all labels per file
undepleted_LabelCounts <- aggregate(groupcount ~ ManualGate + SampleID, transform(undepleted_file, groupcount = 1), length)
undepleted_LabelCountsTable <- dcast(data = undepleted_LabelCounts, formula = ManualGate ~ SampleID)
write.table(undepleted_LabelCountsTable, file = "./3_Output/Barcharts/Thymus_mangated/undepleted_LabelCountsTable_Thymus.txt", sep = "\t", row.names = FALSE)
depleted_LabelCounts <- aggregate(groupcount ~ ManualGate + SampleID, transform(depleted_file, groupcount = 1), length)
depleted_LabelCountsTable <- dcast(data = depleted_LabelCounts, formula = ManualGate ~ SampleID)
write.table(depleted_LabelCountsTable, file = "./3_Output/Barcharts/Thymus_mangated/depleted_LabelCountsTable_Thymus.txt", sep = "\t", row.names = FALSE)


### Calculate percentages
undepleted_PercentagesLabel <- undepleted_LabelCountsTable[1:ncol(undepleted_LabelCountsTable)]
rownames(undepleted_PercentagesLabel) <- undepleted_LabelCountsTable$ManualGate
undepleted_PercentagesLabel<-apply(undepleted_PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(undepleted_PercentagesLabel, file = "./3_Output/Barcharts/Thymus_mangated/undepleted_PercentagesLabel_Thymus.txt", sep = "\t", row.names = T, col.names = NA)
depleted_PercentagesLabel <- depleted_LabelCountsTable[1:ncol(depleted_LabelCountsTable)]
rownames(depleted_PercentagesLabel) <- depleted_LabelCountsTable$ManualGate
depleted_PercentagesLabel<-apply(depleted_PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(depleted_PercentagesLabel, file = "./3_Output/Barcharts/Thymus_mangated/depleted_PercentagesLabel_Thymus.txt", sep = "\t", row.names = T, col.names = NA)


### Reorder the labels
undepleted_PercentagesPerLabel <- melt(as.matrix(undepleted_PercentagesLabel))
colnames(undepleted_PercentagesPerLabel) <- c("ManualGate", "SampleID", "Value")
undepleted_PercentagesPerLabel$ManualGate <- factor(undepleted_PercentagesPerLabel$ManualGate, 
                                         levels = rev(c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                        "B-cell","Monocyte", "NK-cell", "pDC", 
                                                        "Unknown")))
depleted_PercentagesPerLabel <- melt(as.matrix(depleted_PercentagesLabel))
colnames(depleted_PercentagesPerLabel) <- c("ManualGate", "SampleID", "Value")
depleted_PercentagesPerLabel$ManualGate <- factor(depleted_PercentagesPerLabel$ManualGate, 
                                                  levels = rev(c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                                 "B-cell","Monocyte", "NK-cell", "pDC", 
                                                                 "Unknown")))


### Colour variables 
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                "B-cell","Monocyte", "NK-cell", "pDC", 
                "Unknown")


### Convert NA to 0 
undepleted_PercentagesLabel <- as.data.frame(undepleted_PercentagesLabel)
undepleted_PercentagesLabel[is.na(undepleted_PercentagesLabel)] <- 0
depleted_PercentagesLabel <- as.data.frame(depleted_PercentagesLabel)
depleted_PercentagesLabel[is.na(depleted_PercentagesLabel)] <- 0


### Create dendrograms
undepleted_dend <- as.dendrogram(hclust(as.dist(1-cor(undepleted_PercentagesLabel))))
undepleted_dend <- set(undepleted_dend, "labels_cex", 0.7) 
undepleted_PercentagesPerLabel$SampleID <- factor(undepleted_PercentagesPerLabel$SampleID, levels = labels(undepleted_dend))
depleted_dend <- as.dendrogram(hclust(as.dist(1-cor(depleted_PercentagesLabel))))
depleted_dend <- set(depleted_dend, "labels_cex", 0.7)
depleted_PercentagesPerLabel$SampleID <- factor(depleted_PercentagesPerLabel$SampleID, levels = labels(depleted_dend))


### Create plots 
undepleted_p.den <- ggplot(as.ggdend(undepleted_dend), horiz = TRUE, theme = theme_classic()) 
undepleted_p.bar <- ggplot(undepleted_PercentagesPerLabel, aes(x = SampleID, y = Value, fill=ManualGate)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/Thymus_mangated/undepleted_Dendrogram_Thymus.pdf", undepleted_p.den, width = 7, height = 5)
ggsave(filename = "./3_Output/Barcharts/Thymus_mangated/undepleted_OrderedBargraph_Thymus.pdf", undepleted_p.bar, width = 7, height = 5)
depleted_p.den <- ggplot(as.ggdend(depleted_dend), horiz = TRUE, theme = theme_classic()) 
depleted_p.bar <- ggplot(depleted_PercentagesPerLabel, aes(x = SampleID, y = Value, fill=ManualGate)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/Thymus_mangated/depleted_Dendrogram_Thymus.pdf", depleted_p.den, width = 7, height = 5)
ggsave(filename = "./3_Output/Barcharts/Thymus_mangated/depleted_OrderedBargraph_Thymus.pdf", depleted_p.bar, width = 7, height = 5)







