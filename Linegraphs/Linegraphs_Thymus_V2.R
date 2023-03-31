#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
#BiocManager::install("")

### Load data  
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/")
all.df <- readRDS("./CombinedThymi.rds")

sub.df <- all.df[! all.df$Celltype %in% c("B-cell", "Monocyte", "cDC", "NK-cell", "pDC", "Unknown"),]

### Redo everything and get the proportions 
LabelCounts <- aggregate(groupcount ~ Celltype + SampleID, transform(sub.df, groupcount = 1), length)
LabelCountsTable <- dcast(data = LabelCounts, formula = Celltype ~ SampleID)


prob_table <- apply(LabelCountsTable[,-1],2,function(x){x/sum(x, na.rm = T)})
rownames(prob_table) <- LabelCountsTable$Celltype
prob_table <- as.data.frame(prob_table)
prob_table$Celltype <- LabelCountsTable$Celltype
#plot_file <- as.data.frame(t(prob_table))

prob_table <- prob_table[! rownames(prob_table) %in% c("B-cell", "Monocyte"),]
prob_table[is.na(prob_table)] <- 0


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Thymus_FirstTry.pdf")
ggplot() + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A10_ptThy075_thawed), group = 1, color = "BLUE") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A12_ptThy076_thawed), group = 1, color = "BLUE") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy072_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy073_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy075_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy076_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy079_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy083_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy085_Undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A2_ptThy71_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A4_ptThy073_70Mdepleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A4_ptThy075_70Mdepleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A4_ptThy076_70Mdepleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A4_ptThy079_depleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A4_ptThy083_depleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A4_ptThy085_Depleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A4_ptThy71_depleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A6_ptThy070_undepleted), group = 1, color = "RED") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A6_ptThy071_thawed), group = 1, color = "BLUE") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A6_ptThy072_70Mdepleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A6_ptThy073_10Mdepleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A6_ptThy075_10Mdepleted), group = 1, color = "GREEN") + 
  geom_line(prob_table, mapping=aes(rownames(prob_table), A8_ptThy070_depleted), group = 1, color = "GREEN") +
  theme_bw()
dev.off()











pts.df <- readRDS("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Patients/CombinedPatients.rds")
sub.pts.df <- pts.df[! pts.df$Celltype %in% c("B-cell-like", "Monocyte-like", "cDC-like", "NK-cell-like", "pDC-like", "Unknown"),]

### Redo everything and get the proportions 
pts_LabelCounts <- aggregate(groupcount ~ Celltype + SampleID, transform(sub.pts.df, groupcount = 1), length)
pts_LabelCountsTable <- dcast(data = pts_LabelCounts, formula = Celltype ~ SampleID)


pts_prob_table <- apply(pts_LabelCountsTable[,-1],2,function(x){x/sum(x, na.rm = T)})
rownames(pts_prob_table) <- pts_LabelCountsTable$Celltype
pts_prob_table <- as.data.frame(pts_prob_table)
pts_prob_table$Celltype <- pts_LabelCountsTable$Celltype
#plot_file <- as.data.frame(t(prob_table))

#pts_prob_table <- pts_prob_table[! rownames(pts_prob_table) %in% c("B-cell-like", "Monocyte"),]
pts_prob_table[is.na(pts_prob_table)] <- 0







##########################################################################################
library(matrixStats)

C1_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("F12_pt2789" ,"F11_pt2789" ,"A9_pt8148" ,"A8_pt8148" ,"E3_pt540" ,"E2_pt540" ,"G12_pt1946" ,"G11_pt1946" ,"C12_pt1179" ,"C11_pt1179")]
C2_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("D12_pt633" ,"D11_pt633")]
C3_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("A6_pt8173", "A5_pt8173", "E9_pt419", "E8_pt419")]
C4_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("H6_pt5438", "H5_pt5438", "H2_pt1524", "H3_pt1524", "B3_pt14643", "B2_pt14643", "A12_pt7675", "A11_pt7675", "G9_pt1949", "G8_pt1949", "G3_pt2723", "G2_pt2723", "D6_pt2229", "D5_pt2229", "A2_pt8711")]
C5_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("E12_pt335", "E11_pt335", "C6_pt5676", "C5_pt5676", "D9_pt913", "D8_pt913")]
C6_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("E6_pt530", "E5_pt530")]
C7_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("D3_pt2322", "D2_pt2322", "A3_pt8711")]
C8_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("F3_pt258", "F2_pt258", "B12_pt7267", "B11_pt7267", "G6_pt1950", "G5_pt1950")]
C9_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("F6_ptvu10138", "F5_ptvu10138", "B9_pt9175", "B8_pt9175", "A12_pt4564", "A11_pt4564")]
C10_pts <- pts_prob_table[colnames(pts_prob_table) %in% c("C9_pt3976", "C8_pt3976")]

figure_order <- c("DN-like", "ISP-like", "DP-like", "CD4-like", "CD8-like", "gdTcell-like")


mat <-matrix(ncol=length(rownames(C1_pts)),nrow=3)
colnames(mat) <- rownames(C1_pts)
mat[1,] <- rowSds(as.matrix(C1_pts))
mat[2,] <- rowMedians(as.matrix(C1_pts))
mat <- mat[,figure_order]
mat[3,] <- figure_order
rownames(mat) <- c("SD", "Median", "CellTypes")
mat_df <- as.data.frame(t(mat))
mat_df$SD <- as.numeric(mat_df$SD)
mat_df$Median <- as.numeric(mat_df$Median)
mat_df$CellTypes <- factor(mat_df$CellTypes, levels = figure_order)
#mat2 
mat2 <-matrix(ncol=length(rownames(C2_pts)),nrow=3)
colnames(mat2) <- rownames(C2_pts)
mat2[1,] <- rowSds(as.matrix(C2_pts))
mat2[2,] <- rowMedians(as.matrix(C2_pts))
mat2 <- mat2[,figure_order]
mat2[3,] <- figure_order
rownames(mat2) <- c("SD", "Median", "CellTypes")
mat2_df <- as.data.frame(t(mat2))
mat2_df$SD <- as.numeric(mat2_df$SD)
mat2_df$Median <- as.numeric(mat2_df$Median)
mat2_df$CellTypes <- factor(mat2_df$CellTypes, levels = figure_order)
#mat3 
mat3 <-matrix(ncol=length(rownames(C3_pts)),nrow=3)
colnames(mat3) <- rownames(C3_pts)
mat3[1,] <- rowSds(as.matrix(C3_pts))
mat3[2,] <- rowMedians(as.matrix(C3_pts))
mat3 <- mat3[,figure_order]
mat3[3,] <- figure_order
rownames(mat3) <- c("SD", "Median", "CellTypes")
mat3_df <- as.data.frame(t(mat3))
mat3_df$SD <- as.numeric(mat3_df$SD)
mat3_df$Median <- as.numeric(mat3_df$Median)
mat3_df$CellTypes <- factor(mat3_df$CellTypes, levels = figure_order)
#mat4 
mat4 <-matrix(ncol=length(rownames(C4_pts)),nrow=3)
colnames(mat4) <- rownames(C4_pts)
mat4[1,] <- rowSds(as.matrix(C4_pts))
mat4[2,] <- rowMedians(as.matrix(C4_pts))
mat4 <- mat4[,figure_order]
mat4[3,] <- figure_order
rownames(mat4) <- c("SD", "Median", "CellTypes")
mat4_df <- as.data.frame(t(mat4))
mat4_df$SD <- as.numeric(mat4_df$SD)
mat4_df$Median <- as.numeric(mat4_df$Median)
mat4_df$CellTypes <- factor(mat4_df$CellTypes, levels = figure_order)
#mat5 
mat5 <-matrix(ncol=length(rownames(C5_pts)),nrow=3)
colnames(mat5) <- rownames(C5_pts)
mat5[1,] <- rowSds(as.matrix(C5_pts))
mat5[2,] <- rowMedians(as.matrix(C5_pts))
mat5 <- mat5[,figure_order]
mat5[3,] <- figure_order
rownames(mat5) <- c("SD", "Median", "CellTypes")
mat5_df <- as.data.frame(t(mat5))
mat5_df$SD <- as.numeric(mat5_df$SD)
mat5_df$Median <- as.numeric(mat5_df$Median)
mat5_df$CellTypes <- factor(mat5_df$CellTypes, levels = figure_order)
#mat6 
mat6 <-matrix(ncol=length(rownames(C6_pts)),nrow=3)
colnames(mat6) <- rownames(C6_pts)
mat6[1,] <- rowSds(as.matrix(C6_pts))
mat6[2,] <- rowMedians(as.matrix(C6_pts))
mat6 <- mat6[,figure_order]
mat6[3,] <- figure_order
rownames(mat6) <- c("SD", "Median", "CellTypes")
mat6_df <- as.data.frame(t(mat6))
mat6_df$SD <- as.numeric(mat6_df$SD)
mat6_df$Median <- as.numeric(mat6_df$Median)
mat6_df$CellTypes <- factor(mat6_df$CellTypes, levels = figure_order)
#mat7 
mat7 <-matrix(ncol=length(rownames(C7_pts)),nrow=3)
colnames(mat7) <- rownames(C7_pts)
mat7[1,] <- rowSds(as.matrix(C7_pts))
mat7[2,] <- rowMedians(as.matrix(C7_pts))
mat7 <- mat7[,figure_order]
mat7[3,] <- figure_order
rownames(mat7) <- c("SD", "Median", "CellTypes")
mat7_df <- as.data.frame(t(mat7))
mat7_df$SD <- as.numeric(mat7_df$SD)
mat7_df$Median <- as.numeric(mat7_df$Median)
mat7_df$CellTypes <- factor(mat7_df$CellTypes, levels = figure_order)
#mat8 
mat8 <-matrix(ncol=length(rownames(C8_pts)),nrow=3)
colnames(mat8) <- rownames(C8_pts)
mat8[1,] <- rowSds(as.matrix(C8_pts))
mat8[2,] <- rowMedians(as.matrix(C8_pts))
mat8 <- mat8[,figure_order]
mat8[3,] <- figure_order
rownames(mat8) <- c("SD", "Median", "CellTypes")
mat8_df <- as.data.frame(t(mat8))
mat8_df$SD <- as.numeric(mat8_df$SD)
mat8_df$Median <- as.numeric(mat8_df$Median)
mat8_df$CellTypes <- factor(mat8_df$CellTypes, levels = figure_order)
#mat9 
mat9 <-matrix(ncol=length(rownames(C9_pts)),nrow=3)
colnames(mat9) <- rownames(C9_pts)
mat9[1,] <- rowSds(as.matrix(C9_pts))
mat9[2,] <- rowMedians(as.matrix(C9_pts))
mat9 <- mat9[,figure_order]
mat9[3,] <- figure_order
rownames(mat9) <- c("SD", "Median", "CellTypes")
mat9_df <- as.data.frame(t(mat9))
mat9_df$SD <- as.numeric(mat9_df$SD)
mat9_df$Median <- as.numeric(mat9_df$Median)
mat9_df$CellTypes <- factor(mat9_df$CellTypes, levels = figure_order)
#mat10 
mat10 <-matrix(ncol=length(rownames(C10_pts)),nrow=3)
colnames(mat10) <- rownames(C10_pts)
mat10[1,] <- rowSds(as.matrix(C10_pts))
mat10[2,] <- rowMedians(as.matrix(C10_pts))
mat10 <- mat10[,figure_order]
mat10[3,] <- figure_order
rownames(mat10) <- c("SD", "Median", "CellTypes")
mat10_df <- as.data.frame(t(mat10))
mat10_df$SD <- as.numeric(mat10_df$SD)
mat10_df$Median <- as.numeric(mat10_df$Median)
mat10_df$CellTypes <- factor(mat10_df$CellTypes, levels = figure_order)


#c("#88CBED", "#43AA98", "#137734", "#322A80", "#DDCC77", "#999830", "#CB6677", "#872255", "#A74592", "#DCDCDC") order colors from: 1-10
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Extended_Combined.pdf")
ggplot() + 
  geom_line(mat_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#88CBED") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat_df, 
                width=.2, position=position_dodge(0.05), color = "#88CBED") +
  geom_line(mat2_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#43AA98") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat2_df, 
                width=.2, position=position_dodge(0.05), color = "#43AA98") +
  geom_line(mat3_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#137734") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat3_df, 
                width=.2, position=position_dodge(0.05), color = "#137734") +
  geom_line(mat4_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#322A80") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat4_df, 
                width=.2, position=position_dodge(0.05), color = "#322A80") +
  geom_line(mat5_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#DDCC77") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat5_df, 
                width=.2, position=position_dodge(0.05), color = "#DDCC77") +
  geom_line(mat6_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#999830") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat6_df, 
                width=.2, position=position_dodge(0.05), color = "#999830") +
  geom_line(mat7_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#CB6677") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat7_df, 
                width=.2, position=position_dodge(0.05), color = "#CB6677") +
  geom_line(mat8_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#872255") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat8_df, 
                width=.2, position=position_dodge(0.05), color = "#872255") +
  geom_line(mat9_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#A74592") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat9_df, 
                width=.2, position=position_dodge(0.05), color = "#A74592") +
  geom_line(mat10_df, mapping=aes(x = CellTypes, y = Median), group = 1, color = "#DCDCDC") + 
  geom_errorbar(aes(x = CellTypes, y = Median, ymin=Median-SD, ymax=Median+SD), data = mat10_df, 
                width=.2, position=position_dodge(0.05), color = "#DCDCDC") +
  theme_bw()
dev.off()

typeof(mat2_df$SD)
typeof(mat2_df$Median)
ggplot() + 
  geom_line(mat_df, mapping=aes(x = rownames(mat_df), y = Median), group = 1, color = "#88CBED") + 
  geom_errorbar(aes(x = rownames(mat_df), y = Median, ymin=Median-SD, ymax=Median+SD), data = mat_df, 
                width=.2, position=position_dodge(0.05), color = "#88CBED") +
  geom_line(mat2_df, mapping=aes(x = rownames(mat2_df), y = Median), group = 1, color = "#43AA98") + 
  geom_errorbar(aes(x = rownames(mat2_df), y = Median, ymin=Median-SD, ymax=Median+SD), data = mat2_df, 
                width=.2, position=position_dodge(0.05), color = "#43AA98") +
  theme_bw()



pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Extended_Patients_Cluster1.pdf")
ggplot() + 
  geom_line(plot_input, mapping=aes(x = rownames(plot_input), y = Median), group = 1) + 
  geom_errorbar(aes(x = rownames(plot_input), y = Median, ymin=Median-SD, ymax=Median+SD), data = plot_input, 
                width=.2, position=position_dodge(0.05)) +
  theme_bw()
dev.off()


ggplot(plot_input) + geom_errorbar(aes(ymin=Median-SD, ymax=Median+SD),data = plot_input, width=.2,
                         position=position_dodge(0.05))

##########################################################################################

# Plot in order of clusters 
#c("#88CBED", "#43AA98", "#137734", "#322A80", "#DDCC77", "#999830", "#CB6677", "#872255", "#A74592", "#DCDCDC") op volgorde van 1-10
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster1.pdf")
ggplot() + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F12_pt2789), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F11_pt2789), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A9_pt8148), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A8_pt8148), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E3_pt540), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E2_pt540), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G12_pt1946), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G11_pt1946), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C12_pt1179), group = 1, color = "#88CBED") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C11_pt1179), group = 1, color = "#88CBED") + 
  theme_bw()
dev.off()

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster2.pdf")
ggplot() + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D12_pt633), group = 1, color = "#43AA98") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D11_pt633), group = 1, color = "#43AA98") + 
  theme_bw()
dev.off()




pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster3.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A6_pt8173	), group = 1, color = "#137734") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A5_pt8173	), group = 1, color = "#137734") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E9_pt419	), group = 1, color = "#137734") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E8_pt419	), group = 1, color = "#137734") + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster4.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H6_pt5438	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H5_pt5438	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H2_pt1524	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H3_pt1524	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B3_pt14643	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B2_pt14643	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A12_pt7675	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A11_pt7675	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G9_pt1949	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G8_pt1949	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G3_pt2723	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G2_pt2723	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D6_pt2229	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D5_pt2229	), group = 1, color = "#322A80") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A2_pt8711	), group = 1, color = "#322A80") + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster5.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E12_pt335	), group = 1, color = "#DDCC77") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E11_pt335	), group = 1, color = "#DDCC77") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C6_pt5676	), group = 1, color = "#DDCC77") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C5_pt5676	), group = 1, color = "#DDCC77") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D9_pt913	), group = 1, color = "#DDCC77") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D8_pt913	), group = 1, color = "#DDCC77") + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster6.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E6_pt530	), group = 1, color = "#999830") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E5_pt530	), group = 1, color = "#999830") + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster7.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D3_pt2322	), group = 1, color = "#CB6677") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D2_pt2322	), group = 1, color = "#CB6677") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A3_pt8711	), group = 1, color = "#CB6677") + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster8.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F3_pt258	), group = 1, color = "#872255") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F2_pt258	), group = 1, color = "#872255") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B12_pt7267	), group = 1, color = "#872255") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B11_pt7267	), group = 1, color = "#872255") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G6_pt1950	), group = 1, color = "#872255") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G5_pt1950	), group = 1, color = "#872255") + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster9.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F6_ptvu10138	), group = 1, color = "#A74592") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F5_ptvu10138	), group = 1, color = "#A74592") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B9_pt9175	), group = 1, color = "#A74592") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B8_pt9175	), group = 1, color = "#A74592") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A12_pt4564	), group = 1, color = "#A74592") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A11_pt4564	), group = 1, color = "#A74592") + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_Cluster10.pdf")
ggplot() +
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C9_pt3976	), group = 1, color = "#DCDCDC") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C8_pt3976	), group = 1, color = "#DCDCDC") + 
  theme_bw()
dev.off()


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/Patients_FirstTry.pdf")
ggplot() + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A11_pt4564), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A11_pt7675), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A12_pt4564), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A12_pt7675), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A2_pt8711), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A3_pt8711), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A5_pt8173), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A6_pt8173), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A8_pt8148), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), A9_pt8148), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B11_pt7267), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B12_pt7267), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B2_pt14643), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B3_pt14643), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B8_pt9175), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), B9_pt9175), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C11_pt1179), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C12_pt1179), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C5_pt5676), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C6_pt5676), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C8_pt3976), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), C9_pt3976), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D11_pt633), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D12_pt633), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D2_pt2322), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D3_pt2322), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D5_pt2229), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D6_pt2229), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D8_pt913), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), D9_pt913), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E11_pt335), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E12_pt335), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E2_pt540), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E3_pt540), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E5_pt530), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E6_pt530), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E8_pt419), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), E9_pt419), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F11_pt2789), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F12_pt2789), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F2_pt258), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F3_pt258), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F5_ptvu10138), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), F6_ptvu10138), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G11_pt1946), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G12_pt1946), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G2_pt2723), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G3_pt2723), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G5_pt1950), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G6_pt1950), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G8_pt1949), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), G9_pt1949), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H2_pt1524), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H3_pt1524), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H5_pt5438), group = 1, color = "BLUE") + 
  geom_line(pts_prob_table, mapping=aes(rownames(pts_prob_table), H6_pt5438), group = 1, color = "BLUE") + 
  theme_bw()
dev.off()
colnames(pts_prob_table)

######### Old 

colnames(prob_table)


observations <- observations %>% 
  pivot_longer(-Individuals, names_to = "name", values_to = "value")

ggplot(observations, aes(name, Celltype, color = value)) +
  geom_point() #+
  #scale_color_manual(values = c("0" = "white", "1" = "black"))


PercentagesLabel <- LabelCountsTable[1:ncol(LabelCountsTable)]
rownames(PercentagesLabel) <- LabelCountsTable$Celltype
PercentagesLabel<-apply(PercentagesLabel[,-1],2,function(x){x/sum(x, na.rm = T)})



rowSums(PercentagesLabel)





#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



median_persample <- all.df %>% group_by(Celltype, SampleID) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), median)
median_persample <- median_persample[! median_persample$Celltype %in% c("B-cell", "Monocyte"),]
persample_data <- as.data.frame(sapply(median_persample[3:19], function(x) ( (x-min(x))/(max(x)-min(x))) ))
persample_data$Celltype <- median_persample$Celltype
persample_data$SampleID <- median_persample$SampleID

sd.data <- persample_data %>% group_by(Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), sd)
plot_data <- persample_data %>% group_by(Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), median)
#colnames(sd.data) <- paste0("sd_", colnames(sd.data))

colors = c("CD45"="Green", 
           "CD7"="Yellow", 
           "CD5"="Red", 
           "CD1a"="Orange", 
           "CD3"="Pink",
           "TCRg_d"="Gold",
           "CD4"="Black",
           "CD8a"="steelblue",
           "TCRa_b"="Purple")

ggplot() +
  geom_point(persample_data, mapping=aes(Celltype, CD4), color = "Black") +
  geom_line(plot_data, mapping=aes(Celltype, CD4), group = 1, color = "Black") +
  geom_point(persample_data, mapping=aes(Celltype, CD8a), color = "steelblue") +
  geom_line(plot_data, mapping=aes(Celltype, CD8a), group = 1, color = "steelblue")+
  geom_point(persample_data, mapping=aes(Celltype, CD3), color = "Pink") +
  geom_line(plot_data, mapping=aes(Celltype, CD3), group = 1, color = "Pink")#+ 
  #geom_errorbar(plot_data, mapping=aes(ymin=CD4-sd.data$CD4, ymax=CD4+sd.data$CD4), width=.2,
   #             position=position_dodge(0.05))
  

ggplot(persample_data, aes(Celltype, CD4)) + 
  geom_point() +
  geom_line(plot_data$Celltype, plot_data$CD4)
  #geom_line(plot_data, aes(Celltype, CD4))


ggplot(plot_data, aes(Celltype, CD4)) + 
  #geom_point() +
  geom_line( group = 1)
  #geom_line(aes(group = SampleID), stat = "summary", fun = median)# + 
  #
  #geom_errorbar(aes(ymin=CD4-sd.data$sd_CD4, ymax=CD4+sd.data$sd_CD4), width=.2,
   #             position=position_dodge(0.05))
max(persample_data$CD4)





median.data <- all.df %>% group_by(Celltype, SampleID) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), median)
median.data <- median.data[! median.data$Celltype %in% c("B-cell", "Monocyte"),]
plot.data <- as.data.frame(sapply(median.data[3:19], function(x) ( (x-min(x))/(max(x)-min(x))) ))
plot.data$Celltype <- median.data$Celltype
plot.data$SampleID <- median.data$SampleID

sd.data <- plot.data %>% group_by(Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), sd)

CD4.data <- plot.data[,c("CD4", "SampleID", "Celltype")]
CD4.data$Celltype <- as.factor(CD4.data$Celltype)
#CD4.data_long <- melt(CD4.data, id="Celltype") 

typeof(CD4.data)
CD4.data$CD4 <- as.numeric(unlist(CD4.data$CD4))
#data_summary(CD4.data, varname="CD4", 
  #           groupnames=c("SampleID", "Celltype"))
CD4.data$sd <- sd(unlist(CD4.data["CD4"]), na.rm=TRUE)
#sd(CD4.data["CD4"], na.rm=TRUE)

unique(plot.data$Celltype)


ggplot(CD4.data, aes(Celltype, CD4)) + 
  geom_line(aes(group = SampleID), stat = "summary", fun = median) + 
  #geom_point() +
  geom_errorbar(aes(ymin=CD4-sd, ymax=CD4+sd), width=.2,
                position=position_dodge(0.05))
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
    #            position=position_dodge(0.05))

stage <- c("DN", "gdTcell", "ISP", "DP", "CD4", "CD8", "NK-cell", "pDC", "cDC", "Unknown", "Monocyte")
first_minmax <- as.data.frame(sapply(all.df[c("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d")], function(x) ( (x-min(x))/(max(x)-min(x))) ))
first_minmax$Celltype <- factor(all.df$Celltype, levels = stage)
median.data_minmax <- first_minmax %>% group_by(Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), .funs = lst(median, sd), )
median.data_minmax <- median.data_minmax[! median.data_minmax$Celltype %in% c("B-cell", "Monocyte"),]
ggplot(median.data_minmax, aes(x=Celltype), color = colors) + 
  geom_line(aes(y = CD4_median, group = 1, color = "CD4")) +
  geom_errorbar(aes(ymin=CD4_median-CD4_sd, ymax=CD4_median+CD4_sd), width=.2,
                position=position_dodge(0.05))


median.data_minmax$CD4_median







median.data_new <- all.df %>% group_by(Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), .funs = lst(median, sd), )
colnames(median.data_new)
median.data_new <- median.data_new[! median.data_new$Celltype %in% c("B-cell", "Monocyte"),]
plot.data_new <- as.data.frame(sapply(median.data_new[3:18], function(x) ( (x-min(x))/(max(x)-min(x))) ))

#colnames(median.data_new[19:25])
#colnames(median.data_new[18:25])
plot.data_new <- cbind(plot.data_new, median.data_new[19:25])

plot.data_new <- as.data.frame(sapply(median.data_new[3:35], function(x) ( (x-min(x))/(max(x)-min(x))) ))
plot.data_new$Celltype <- factor(median.data_new$Celltype, levels = stage)
ggplot(plot.data_new, aes(x=Celltype), color = colors) + 
  geom_line(aes(y = CD4_median, group = 1, color = "CD4")) +
  geom_errorbar(aes(ymin=CD4_median-CD4_sd, ymax=CD4_median+CD4_sd), width=.2,
                position=position_dodge(0.05))







######################################################################## test old fashioned summarized plots 
#median.data_old 
stage <- c("DN", "gdTcell", "ISP", "DP", "CD4", "CD8", "NK-cell", "pDC", "cDC", "Unknown")
median.data_old <- all.df %>% group_by(Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), median)
median.data_old <- median.data_old[! median.data_old$Celltype %in% c("B-cell", "Monocyte"),]
plot.data_old <- as.data.frame(sapply(median.data_old[3:18], function(x) ( (x-min(x))/(max(x)-min(x))) ))
plot.data_old$Celltype <- factor(median.data_old$Celltype, levels = stage)
#plot.data_old$SampleID <- median.data_old$SampleID
plot.data_old

#plot.data_old <- plot.data_old[match(stage, plot.data_old$Celltype), ]
#plot.data_old$Celltype <- as.factor(plot.data_old$Celltype)


colors = c("CD45"="Green", 
           "CD7"="Yellow", 
           "CD5"="Red", 
           "CD1a"="Orange", 
           "CD3"="Pink",
           "TCRg_d"="Gold",
           "CD4"="Black",
           "CD8a"="steelblue",
           "TCRa_b"="Purple")
#pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/CombinedThymus_line_simple.pdf")
ggplot(plot.data_old, aes(x=Celltype), color = colors) + 
  geom_line(aes(y = CD4, group = 1, color = "CD4")) +
  geom_line(aes(y = CD8a, group = 1, color = "CD8a")) +
  geom_line(aes(y = CD3, group = 1, color = "CD3")) +
  geom_line(aes(y = CD5, group = 1, color = "CD5")) +
  geom_line(aes(y = CD45, group = 1, color = "CD45")) +
  geom_line(aes(y = CD7, group = 1, color = "CD7")) +
  geom_line(aes(y = CD1a, group = 1, color = "CD1a")) +
  geom_line(aes(y = TCRa_b, group = 1, color = "TCRa_b")) +
  geom_line(aes(y = TCRg_d, group = 1, color = "TCRg_d")) + 
  #geom_point() +
  labs(x = "Cell type", y = "Relative expression", color = "Legend") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
#dev.off()