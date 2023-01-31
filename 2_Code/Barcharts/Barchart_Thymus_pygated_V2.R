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


### Count the occurence of all labels per file
LabelCounts <- aggregate(groupcount ~ Celltype + SampleID, transform(combine.file, groupcount = 1), length)
LabelCountsTable <- dcast(data = LabelCounts, formula = Celltype ~ SampleID)
write.table(LabelCountsTable, file = "./3_Output/Barcharts/Thymus_pygated/LabelCountsTable_Thymus.txt", sep = "\t", row.names = FALSE)

### Calculate percentages
PercentagesLabel <- LabelCountsTable[1:ncol(LabelCountsTable)]
rownames(PercentagesLabel) <- LabelCountsTable$Celltype
PercentagesLabel<-apply(PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(PercentagesLabel, file = "./3_Output/Barcharts/Thymus_pygated/PercentagesLabel_Thymus.txt", sep = "\t", row.names = T, col.names = NA)

### Reorder the labels
PercentagesPerLabel <- melt(as.matrix(PercentagesLabel))
colnames(PercentagesPerLabel) <- c("Celltype", "SampleID", "Value")
PercentagesPerLabel$Celltype <- factor(PercentagesPerLabel$Celltype, 
                                       levels = rev(c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                      "B-cell","Monocyte", "NK-cell", "pDC", 
                                                      "Unknown")))
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                "B-cell","Monocyte", "NK-cell", "pDC", 
                "Unknown")

### Convert NA to 0 
PercentagesLabel <- as.data.frame(PercentagesLabel)
PercentagesLabel[is.na(PercentagesLabel)] <- 0


### Create dendrograms
dend <- as.dendrogram(hclust(as.dist(1-cor(PercentagesLabel))))
dend <- set(dend, "labels_cex", 0.7) 
PercentagesPerLabel$SampleID <- factor(PercentagesPerLabel$SampleID, levels = labels(dend))


### Create plots 
p.den <- ggplot(as.ggdend(dend), horiz = TRUE, theme = theme_classic()) 
p.bar <- ggplot(PercentagesPerLabel, aes(x = SampleID, y = Value, fill=Celltype)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/Thymus_pygated/Dendrogram_Thymus.pdf", p.den, width = 7, height = 5)
ggsave(filename = "./3_Output/Barcharts/Thymus_pygated/OrderedBargraph_Thymus.pdf", p.bar, width = 7, height = 5)


saveRDS(object = PercentagesLabel, file = "./3_Output/Barcharts/Thymus_pygated/PercentagesLabel.RDS")
saveRDS(object = PercentagesPerLabel, file = "./3_Output/Barcharts/Thymus_pygated/PercentagesPerLabel.RDS")



################################################################## Manually gated ##################################################################
### Count the occurence of all labels per file
LabelCounts <- aggregate(groupcount ~ ManualGate + SampleID, transform(combine.file, groupcount = 1), length)
LabelCountsTable <- dcast(data = LabelCounts, formula = ManualGate ~ SampleID)
write.table(LabelCountsTable, file = "./3_Output/Barcharts/Thymus_mangated/LabelCountsTable_Thymus.txt", sep = "\t", row.names = FALSE)

### Calculate percentages
PercentagesLabel <- LabelCountsTable[1:ncol(LabelCountsTable)]
rownames(PercentagesLabel) <- LabelCountsTable$ManualGate
PercentagesLabel<-apply(PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})
write.table(PercentagesLabel, file = "./3_Output/Barcharts/Thymus_mangated/PercentagesLabel_Thymus.txt", sep = "\t", row.names = T, col.names = NA)

### Reorder the labels
PercentagesPerLabel <- melt(as.matrix(PercentagesLabel))
colnames(PercentagesPerLabel) <- c("ManualGate", "SampleID", "Value")
PercentagesPerLabel$ManualGate <- factor(PercentagesPerLabel$ManualGate, 
                                         levels = rev(c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                        "B-cell","Monocyte", "NK-cell", "pDC", 
                                                        "Unknown")))
pal = c("#33A02C", "#E7298A", "#8b0000", "#B17BA6", "#1965B0", "#FF7F00", "#ffff00", 
        "#606060", "#A0A0A0", "#C8C8C8", "#E8E8E8",
        "#000000")
names(pal) <- c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                "B-cell","Monocyte", "NK-cell", "pDC", 
                "Unknown")

### Convert NA to 0 
PercentagesLabel <- as.data.frame(PercentagesLabel)
PercentagesLabel[is.na(PercentagesLabel)] <- 0


### Create dendrograms
dend <- as.dendrogram(hclust(as.dist(1-cor(PercentagesLabel))))
dend <- set(dend, "labels_cex", 0.7) 
PercentagesPerLabel$SampleID <- factor(PercentagesPerLabel$SampleID, levels = labels(dend))


### Create plots 
p.den <- ggplot(as.ggdend(dend), horiz = TRUE, theme = theme_classic()) 
p.bar <- ggplot(PercentagesPerLabel, aes(x = SampleID, y = Value, fill=ManualGate)) +
  geom_bar(stat='identity')+coord_flip() + theme_classic() + scale_fill_manual(values =pal)
ggsave(filename = "./3_Output/Barcharts/Thymus_mangated/Dendrogram_Thymus.pdf", p.den, width = 7, height = 5)
ggsave(filename = "./3_Output/Barcharts/Thymus_mangated/OrderedBargraph_Thymus.pdf", p.bar, width = 7, height = 5)





################################################################## Get the distribution table ##################################################################
cont_dir <- "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Barcharts/Contingency/"
dir.create(cont_dir)
dist_table <- as.data.frame.matrix(table(combine.file$ManualGate, combine.file$Celltype))
write.table(dist_table, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Barcharts/Contingency/ContingencyTableManVsPy.txt", sep = "\t", col.names = NA)


for (s in unique(combine.file$SampleID)){
  print(s)
  sub_file <- combine.file[combine.file$SampleID == s,  ]
  dist_table <- as.data.frame.matrix(table(sub_file$ManualGate, sub_file$Celltype))
  write.table(dist_table, paste0(cont_dir, s, "_contingencyTableManVsPy.txt"), sep = "\t", col.names = NA)
}





######## Compare the gating strategies 
library(lsa)

dist_table <- as.data.frame.matrix(table(combine.file$ManualGate, combine.file$Celltype))

cos.table <- cosine(as.matrix(dist_table))
write.table(cos.table, paste0(cont_dir, "Cosine_sim_combined.txt"), sep = "\t", col.names = NA)
pdf(paste0(cont_dir, "Cosine_sim_combined_HM.pdf"))
pheatmap(cos.table)
dev.off()

chisq.test(dist_table)
library("gplots")
balloonplot(as.table(dist_table))


library(corrplot)
chisq <- chisq.test(dist_table)
pdf(paste0(cont_dir, "Cor_chisq_residuals.pdf"))
corrplot(chisq$residuals, is.cor = FALSE)
dev.off()
contrib <- 100*chisq$residuals^2/chisq$statistic
pdf(paste0(cont_dir, "Cor_chisq_contribution.pdf"))
corrplot(contrib, is.cor = FALSE)
dev.off()


#library(ISLR)
library(tidyverse)
#library(Rfast)
library(MASS)
loglm( ~ "DP" + "DP", data = dist_table) 
loglm(dist_table[,"DP"], dist_table["DP",])







### Perform the correct statistic
library(caret)
combine.file <- read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Barcharts/Contingency/ContingencyTableManVsPy.txt", 
                           sep = "\t", header = T, row.names = 1)

test <- combine.file
test$Celltype <- factor(test$Celltype, levels = c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                 "B-cell","Monocyte", "NK-cell", "pDC", 
                                                 "Unknown"))
test$ManualGate <- factor(test$ManualGate, levels = c("DN", "ISP", "gdTcell","DP", "CD4", "CD8", "cDC", 
                                                  "B-cell","Monocyte", "NK-cell", "pDC", 
                                                  "Unknown"))
results <- confusionMatrix(test$Celltype, reference = test$ManualGate)
write.table(results$table, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Barcharts/Contingency/Stats/Stats_table.txt", sep = "\t", col.names = NA)
write.table(results$overall, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Barcharts/Contingency/Stats/Stats_overall.txt", sep = "\t", col.names = NA)
write.table(results$byClass, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Barcharts/Contingency/Stats/Stats_byClass.txt", sep = "\t", col.names = NA)








