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


# Exporting FCS files from R - without the issues with truncation of values
# Code by Dmitry Tebaykin
write_modified_FCS <- function(fcs_exprs, fcs_name, channel_descriptions = NULL, 
                               reference_description = NULL, old_description = NULL) {
  if (requireNamespace('flowCore', quietly = TRUE)) {
    if (is.matrix(fcs_exprs)) {
      
      # Build metadata for FCS file
      pd <- c()  # 'params' phenoData
      dl <- list()  # 'description' list
      
      dl[['$DATATYPE']] <- 'F'
      
      if (!is.null(reference_description)) {
        if (!is.null(reference_description[['$DATE']])){
          dl[['$DATE']] <- reference_description[['$DATE']]
        }
        if (!is.null(reference_description[['$BTIM']])){
          dl[['$BTIM']] <- reference_description[['$BTIM']]
        }
        if (!is.null(reference_description[['$ETIM']])){
          dl[['$ETIM']] <- reference_description[['$ETIM']]
        }
        if (!is.null(reference_description[['$CYT']])){
          dl[['$CYT']] <- reference_description[['$CYT']]
        }
        if (!is.null(reference_description[['$CYTSN']])){
          dl[['$CYTSN']] <- reference_description[['$CYTSN']]
        }      
      }
      
      for (c in 1:ncol(fcs_exprs)) {
        c_name <- colnames(fcs_exprs)[c]
        c_desc <- colnames(fcs_exprs)[c]
        if (!is.null(channel_descriptions)){
          if (!is.na(channel_descriptions[c])){
            c_desc <- channel_descriptions[c]  
          }
        }
        
        c_min <- floor(min(0, min(fcs_exprs[, c])) )  # Prevents flowCore from shifting range
        c_max <- ceiling(max(fcs_exprs[, c]))
        c_rng <- c_max - c_min + 1
        
        pl <- matrix(c(c_name, c_desc, c_rng, c_min, c_max), nrow = 1)
        colnames(pl) <- c('name', 'desc', 'range', 'minRange', 'maxRange')
        rownames(pl) <- paste('$P', c, sep = '') 
        pd <- rbind(pd, pl)
        
        if (!is.null(reference_description[[paste0('P', c, 'DISPLAY')]])){
          dl[[paste0('P', c, 'DISPLAY')]] <- reference_description[[paste0('P', c, 'DISPLAY')]]
        } 
        
        if (!is.null(reference_description[[paste0('$P', c, 'G')]])){
          dl[[paste0('$P', c, 'G')]] <- reference_description[[paste0('$P', c, 'G')]]
        } 
        
        if (!is.null(reference_description[[paste0('$P', c, 'R')]])){
          dl[[paste0('$P', c, 'R')]] <- reference_description[[paste0('$P', c, 'R')]]
        } else {
          dl[[paste('$P', c, 'R',sep = '')]] <- toString(c_rng); # Range
        }
        
        if (!is.null(reference_description[[paste0('$P', c, 'B')]])){
          dl[[paste0('$P', c, 'B')]] <- reference_description[[paste0('$P', c, 'B')]]
        } else {
          dl[[paste('$P', c, 'B', sep = '')]] <- '32';      # Number of bits
        }
        
        if (!is.null(reference_description[[paste0('$P', c, 'E')]])){
          dl[[paste0('$P', c, 'E')]] <- reference_description[[paste0('$P', c, 'E')]]
        } else { 
          dl[[paste('$P', c, 'E', sep = '')]] <- '0,0';      # Exponent
        }
        
        dl[[paste('$P', c, 'N', sep = '')]] <- c_name;	    # Name
        dl[[paste('$P', c, 'S', sep = '')]] <- c_desc;	    # Desc	
      }	
      
      if (!is.null(old_description)) {
        if (!is.null(old_description[['$CYT']])){
          dl[['$CYT']] <- old_description[['$CYT']]
        }		  
        if (!is.null(old_description[['$DATE']])){
          dl[['$DATE']] <- old_description[['$DATE']]
        }
        if (!is.null(old_description[['$BTIM']])){
          dl[['$BTIM']] <- old_description[['$BTIM']]
        }
        if (!is.null(old_description[['$ETIM']])){
          dl[['$ETIM']] <- old_description[['$ETIM']]
        }
      }
      
      fcs_exprs <- flowCore::flowFrame(fcs_exprs, 
                                       as(data.frame(pd), 'AnnotatedDataFrame'), 
                                       description = dl) 
    }
    
    suppressWarnings(flowCore::write.FCS(fcs_exprs, fcs_name, what = 'numeric'))
  }
}


### create flow Frame from Data frame 
# from: https://rdrr.io/github/ssmpsn2/flowAssist/src/R/DFtoFF.R
DFtoFF<-function(DF){
  if(class(DF) == "data.frame"){
    return(flowFrame(as.matrix(DF)))
  }
  if(class(DF) == "list"){
    frameList<-as.list(NULL)
    length(frameList)<-length(DF)
    for(i in 1:length(DF)){
      if(class(DF[[i]]) == "data.frame"){
        frameList[[i]]<-flowFrame(as.matrix(DF[[i]]))
        names(frameList)[[i]]<-names(DF)[[i]]
      }
      else{
        warning(paste("Object at index",i,"not of type data.frame"))
      }
    }
    return(frameList)
  }
  else {
    stop(" Object is not of type data.frame")
  }
}


### Load data 
all.df <- readRDS("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/ThymusPatients/ThyPtsCombined.rds")
all.df$TissueType <- ifelse(grepl("Thy",all.df$SampleIDshort),'Thymus','Patient')
thy.df <- all.df[all.df$TissueType == "Thymus", ]
pts.df <- all.df[all.df$TissueType == "Patient", ]



pts.df.sub <- pts.df[pts.df$SampleIDshort == "pt335", ]
fsom.df <- rbind(thy.df, pts.df.sub)
fsom.df.sub <- fsom.df[,1:17] %>% select(- contains("CD45"))
fSOM.input <- DFtoFF(fsom.df.sub)

fSOM.10 <- FlowSOM(fSOM.input, colsToUse = colnames(fsom.df.sub), nClus = 10)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FlowSOM_CellType_Meta_pt335.pdf")
PlotPies(fSOM.10, cellTypes = fsom.df$Celltype, backgroundValues = fSOM.10$metaclustering)
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FlowSOM_SampleName_Meta_pt335.pdf")
PlotPies(fSOM.10, cellTypes = fsom.df$SampleIDshort, backgroundValues = fSOM.10$metaclustering)
dev.off()

#percentages <- GetFeatures(fsom = flowSOM.res, files = paste0("ff_tmp",1:10,".fcs"), type = "percentages")
percentages <- GetFeatures(fsom = fSOM.10, type = "percentages", files = fsom.df$SampleIDshort)
counts <- GetFeatures(fsom = fSOM.10, type = "percentages", files = c(fSOM.input[1:length(rownames(thy.df)), ], fSOM.input[length(rownames(thy.df)):length(rownames(thy.df))+10000, ]) )
#c(ff[1001:2000, ], ff[2001:3000, ]


groups <- list("Group 1" = 1,"Group 2" = 2)

MC_stats <- GroupStats(counts[["metacluster_percentages"]], groups )
C_stats <- GroupStats(counts[["cluster_percentages"]], groups)

fold_changes <- C_stats["fold changes", ]


gr_1 <- PlotStars(fSOM.10, title = "Group 1",
                  nodeSizes = C_stats["medians Group 1", ],
                  refNodeSize = max(C_stats[c("medians Group 1", "medians Group 2"),]),
                  backgroundValues = fold_changes,
                  backgroundColors = c("white", "red", "blue"),
                  list_insteadof_ggarrange = TRUE)
gr_2 <- PlotStars(fSOM.10, title = "Group 2",
                    nodeSizes = C_stats["medians Group 2", ],
                    refNodeSize = max(C_stats[c("medians Group 1", "medians Group 2"),]),
                    backgroundValues = fold_changes,
                    backgroundColors = c("white", "red", "blue"),
                    list_insteadof_ggarrange = TRUE)


p <- ggpubr::ggarrange(plotlist = list(gr_1$tree, gr_2$tree,
                                       gr_2$starLegend, gr_2$backgroundLegend),
                       ncol = 2, nrow = 2,
                       heights = c(4, 1))
print(p, newpage = FALSE)


fSOM.10$data
fSOM.10$metaData
fSOM.10$map
fSOM.10$MST
fSOM.10$map





########################################################################### Test new set-up
### Input data 
thy.df <- all.df[all.df$TissueType == "Thymus", ]
thy.df <- thy.df[thy.df$SampleIDshort != "ptThy073", ]
thy.df.sub <- all.df[all.df$SampleID == "A4_ptThy073_70Mdepleted", ]
pts.df.sub <- all.df[all.df$SampleIDshort == "pt335", ]


### Convert to flowFrame
thy.df.sub.ff <- thy.df.sub[,1:17] %>% select(- contains("CD45"))
thy.df.sub.ff <- DFtoFF(thy.df.sub.ff)
write_modified_FCS(thy.df.sub.ff, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/Thy073.fcs")
pts.df.sub.ff <- pts.df.sub[,1:17] %>% select(- contains("CD45"))
pts.df.sub.ff <- DFtoFF(pts.df.sub.ff)
write_modified_FCS(pts.df.sub.ff, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/pt335.fcs")


### Create FlowSOM
thy.df.ff <- thy.df[,1:17] %>% select(- contains("CD45"))
thy.df.ff <- DFtoFF(thy.df.ff)
fSOM.10 <- FlowSOM(thy.df.ff, colsToUse = colnames(thy.df.ff), nClus = 10)


# Get the count matrix
percentages <- GetFeatures(fsom = fSOM.10, files = c("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/Thy073.fcs",
                                                     "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/pt335.fcs"), type = "percentages")
# Perform the statistics
groups <- list("Group 1" = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/Thy073.fcs",
               "Group 2" = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/pt335.fcs")
MC_stats <- GroupStats(percentages[["metacluster_percentages"]], groups)
C_stats <- GroupStats(percentages[["cluster_percentages"]], groups)
# Process the fold changes vector
fold_changes <- C_stats["fold changes", ]
fold_changes <- factor(ifelse(fold_changes < -3, "Underrepresented compared to Group 1",
                              ifelse(fold_changes > 3, "Overrepresented compared to Group 1", 
                                     "--")), levels = c("--", "Underrepresented compared to Group 1",
                                                        "Overrepresented compared to Group 1"))
fold_changes[is.na(fold_changes)] <- "--"
# Show in figure
## Fold change
gr_1 <- PlotStars(fSOM.10, title = "Group 1",
                  nodeSizes = C_stats["medians Group 1", ],
                  refNodeSize = max(C_stats[c("medians Group 1", "medians Group 2"),]),
                  backgroundValues = fold_changes,
                  backgroundColors = c("white", "red", "blue"),
                  list_insteadof_ggarrange = TRUE)
gr_2 <- PlotStars(fSOM.10, title = "Group 2",
                  nodeSizes = C_stats["medians Group 2", ],
                  refNodeSize = max(C_stats[c("medians Group 1", "medians Group 2"),]),
                  backgroundValues = fold_changes,
                  backgroundColors = c("white", "red", "blue"),
                  list_insteadof_ggarrange = TRUE)
p <- ggpubr::ggarrange(plotlist = list(gr_1$tree, gr_2$tree,
                                       gr_2$starLegend, gr_2$backgroundLegend),
                       ncol = 2, nrow = 2,
                       heights = c(4, 1))
print(p, newpage = FALSE)

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/SinglePtsVsThymus073/pt335_Thy073_diff.pdf")
print(p, newpage = TRUE)
dev.off()









### Create FlowSOM
all.df.ff <- all.df[,1:17] %>% select(- contains("CD45"))
all.df.ff <- DFtoFF(all.df.ff)
fSOM.all.10 <- FlowSOM(all.df.ff, colsToUse = colnames(all.df.ff), nClus = 10)

# Get the count matrix
percentages.all <- GetFeatures(fsom = fSOM.all.10, files = c("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/Thy073.fcs",
                                                     "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/pt335.fcs"), type = "percentages")
# Perform the statistics
groups <- list("Group 1" = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/Thy073.fcs",
               "Group 2" = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/FcsFiles/pt335.fcs")
MC_stats.all <- GroupStats(percentages[["metacluster_percentages"]], groups)
C_stats.all <- GroupStats(percentages[["cluster_percentages"]], groups)
# Process the fold changes vector
fold_changes.all <- C_stats.all["fold changes", ]
fold_changes.all <- factor(ifelse(fold_changes.all < -3, "Underrepresented compared to Group 1",
                              ifelse(fold_changes.all > 3, "Overrepresented compared to Group 1", 
                                     "--")), levels = c("--", "Underrepresented compared to Group 1",
                                                        "Overrepresented compared to Group 1"))
fold_changes.all[is.na(fold_changes.all)] <- "--"
# Show in figure
## Fold change
gr_1.all <- PlotStars(fSOM.all.10, title = "Group 1",
                  nodeSizes = C_stats.all["medians Group 1", ],
                  refNodeSize = max(C_stats.all[c("medians Group 1", "medians Group 2"),]),
                  backgroundValues = fold_changes.all,
                  backgroundColors = c("white", "red", "blue"),
                  list_insteadof_ggarrange = TRUE)
gr_2.all <- PlotStars(fSOM.all.10, title = "Group 2",
                  nodeSizes = C_stats.all["medians Group 2", ],
                  refNodeSize = max(C_stats.all[c("medians Group 1", "medians Group 2"),]),
                  backgroundValues = fold_changes.all,
                  backgroundColors = c("white", "red", "blue"),
                  list_insteadof_ggarrange = TRUE)
p.all <- ggpubr::ggarrange(plotlist = list(gr_1.all$tree, gr_2.all$tree,
                                       gr_2.all$starLegend, gr_2.all$backgroundLegend),
                       ncol = 2, nrow = 2,
                       heights = c(4, 1))


print(p.all, newpage = TRUE)


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/SinglePtsVsThymus073/pt335_Thy073_diff_all.pdf")
print(p.all, newpage = TRUE)
dev.off()





PlotPies(fSOM.all.10, cellTypes = all.df$SampleIDshort, backgroundValues = fSOM.all.10$metaclustering)
PlotPies(fSOM.all.10, cellTypes = ifelse(all.df$SampleIDshort %in% c("ptThy073", "pt335"), all.df$SampleIDshort, "Other"), backgroundValues = fSOM.all.10$metaclustering)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/SinglePtsVsThymus073/fSOM_SampleContribution.pdf")
PlotPies(fSOM.all.10, cellTypes = ifelse(all.df$SampleIDshort %in% c("ptThy073", "pt335"), all.df$SampleIDshort, "Other"))
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/SinglePtsVsThymus073/fSOM_all.pdf")
PlotStars(fSOM.all.10 )
dev.off()
PlotPies(fSOM.10, cellTypes = thy.df$SampleIDshort, backgroundValues = fSOM.10$metaclustering)




#thy.df <- thy.df[thy.df$SampleIDshort != "ptThy073", ]
#pts.df.sub <- all.df[all.df$SampleIDshort == "pt335", ]







########################################################################### Test adding data to existing tree
thy.dat <- all.df[all.df$TissueType == "Thymus", ]
thy.dat.ff <- thy.dat[,1:17] %>% select(- contains("CD45"))
thy.dat.ff <- DFtoFF(thy.dat.ff)
pts.dat <- all.df[all.df$TissueType == "Patient", ]
pts.dat.ff <- pts.dat[,1:17] %>% select(- contains("CD45"))
pts.dat.ff <- DFtoFF(pts.dat.ff)


fSOM.thy.10 <- FlowSOM(thy.dat.ff, colsToUse = colnames(thy.dat.ff), nClus = 10, seed = 42)


fSOM.merged.10 <- NewData(fsom = fSOM.thy.10, 
                          input = pts.dat.ff)
#Warning message:
#  In NewData(fsom = fSOM.thy.10, input = pts.dat.ff) :
#  27860 cells (10.03%) seem far from their cluster centers.

f1 <- PlotStars(fSOM.thy.10)
f2 <- PlotStars(fSOM.merged.10)

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/fSOM_Thyall.pdf")
f1
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/fSOM_Thyall_merged.pdf")
f2
dev.off()

fSOM.pt335.10 <- NewData(fsom = fSOM.thy.10, 
                          input = pts.df.sub.ff)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/fSOM_Thyall_pt335.pdf")
PlotStars(fSOM.pt335.10)
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/FlowSOM/fSOM_Thyall_celltypes.pdf")
PlotPies(fSOM.thy.10, cellTypes = thy.dat$Celltype, backgroundValues = fSOM.thy.10$metaclustering)
dev.off()









