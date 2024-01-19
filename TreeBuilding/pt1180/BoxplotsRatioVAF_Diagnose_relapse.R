#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(rstatix)

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt1180/branches_vcfs/overlapBulkSeq/")

### Colour selection
pal = c("#FFA500","#FFA500","#33A02C","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4", "#A4A4A4", "#385725", "#385725","#E7298A")
names(pal) <- c("Dx_ALLBULK", "Rx_ALLBULK","DN", "Dx_DNCD1aNeg","Rx_DNCD1aNeg","Dx_SPCD4", "DP", "Dx_MSCBULK", "MONOBulk","Dx_DNCD1aPos","Rx_DNCD1aPos", "Dx_iSPCD4")
#SampleOrder <- c("MONOBULK", "ALLBULK", "DN", "SPCD4","DP", "iSPCD4")
SampleOrder <- c("Dx_ALLBULK", "Dx_MSCBULK","Dx_DNCD1aNeg", "Dx_DNCD1aPos", "Dx_iSPCD4", "Dx_SPCD4", 
                 "Rx_ALLBULK", "Rx_DNCD1aNeg", "Rx_DNCD1aPos")

## load df molten file for relapse as list files
#diagnosis_files <- list.files(path = './plots/', pattern= '.csv', full.names = T)
#relapse_files <- list.files(path = './../overlapRelapse/plots', pattern = ".csv", full.names = T)


for (csv_file in list.files(path='./plots',pattern='.csv',full.names = F)) {
  if (file.exists(paste0('../overlapRelapse/plots/', csv_file))) {
    mydf2 <- read.csv(paste0('../overlapRelapse/plots/', csv_file))
  } else { next
  }
  mydf <- read.csv(paste0('./plots/', csv_file))
  
  pdf_fname <- gsub("df.csv","VAFplot.pdf",csv_file)
  pdf_fname <- gsub("./","",pdf_fname)
  plot_title <- gsub("VAFplot.pdf", "", pdf_fname)
  
  mydf$bulk <- paste0("Dx_", mydf$bulk)
  mydf2$bulk <- paste0("Rx_", mydf2$bulk)
  mydf_comb <- rbind(mydf, mydf2)
  
  # T-test
  n1 <- ggqqplot(mydf_comb, "value", facet.by = "bulk")
  
  mydf_comb.filter <- mydf_comb %>% filter(bulk %in% c("Dx_ALLBULK", "Dx_DNCD1aNeg", "Dx_DNCD1aPos", "Dx_iSPCD4", "Dx_SPCD4", 
                                                       "Rx_ALLBULK", "Rx_DNCD1aNeg", "Rx_DNCD1aPos"))
  res.aov <- anova_test(data = mydf_comb, dv = value, wid = variant, within = bulk) 
  get_anova_table(res.aov)

  pwc <- mydf_comb %>% pairwise_t_test(value ~ bulk, paired = FALSE,
                                       p.adjust.method = "bonferroni") # paired was not possible as in relapse other variants are round to overlap than in diagnosis, and they are not in the same VCF file. 
  
  # Visualization: box plots with p-values
  pwc <- pwc %>% add_xy_position(x = "value")
  pwc1 <- pwc %>% filter(group1 %in% c("Dx_ALLBULK"))
  
  p <- ggplot(mydf_comb, aes(bulk, value, color = bulk, fill = bulk)) +
    geom_boxplot(alpha = 0.4) +
    theme_classic() +
    geom_hline(yintercept=1, linetype="dotted")+
    scale_color_manual(values=pal)+ 
    scale_fill_manual(values=pal)+
    scale_x_discrete(limits=SampleOrder)+
    ggtitle(plot_title)+
    ylab("observed/expected VAF")+
    labs(subtitle = get_test_label(res.aov, detailed = TRUE),
         caption = get_pwc_label(pwc))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p1 <- p +
    stat_pvalue_manual(pwc1, xmin = "group1", xmax = "group2", inherit.aes = F) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc1)
    )
  
  p2 <- p + 
    stat_pvalue_manual(pwc, xmin = "group1", xmax = "group2", inherit.aes = F) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc)
    )
  
  pdf(paste0("../overlapDiagnosisRelapse/",pdf_fname))
  print(n1)
  print(p)
  print(p1)
  #print(p2)
  dev.off()
  
} 

#getwd()








library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(rstatix)


setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt1180/branches_vcfs/overlapBulkSeq/")
#dir.create("OERatio_var", showWarnings = TRUE, recursive = FALSE, mode = "0777")
#colors
### Colour selection
pal = c("#FFA500","#FFA500","#33A02C","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4", "#A4A4A4", "#385725", "#385725","#E7298A")
names(pal) <- c("Dx_ALLBULK", "Rx_ALLBULK","DN", "Dx_DNCD1aNeg","Rx_DNCD1aNeg","Dx_SPCD4", "DP", "Dx_MSCBULK", "MONOBulk","Dx_DNCD1aPos","Rx_DNCD1aPos", "Dx_iSPCD4")
#SampleOrder <- c("MONOBULK", "ALLBULK", "DN", "SPCD4","DP", "iSPCD4")
SampleOrder <- c("Dx_ALLBULK", "Dx_MSCBULK","Dx_DNCD1aNeg", "Dx_DNCD1aPos", "Dx_iSPCD4", "Dx_SPCD4",
                 "Rx_ALLBULK","Rx_DNCD1aNeg", "Rx_DNCD1aPos")
#list1<-list.files(path='.',pattern='_branch.vcf',full.names = T)

#test <- readVcf('../overlapRelapse/', vcf_fname)
#paste0('../overlap')
for (vcf_fname in list.files(path='.',pattern='_BulkSeqSMuRF.vcf',full.names = F)) {
  myvcf <- readVcf(vcf_fname)
  myvcf2 <- readVcf(paste0('../overlapRelapse/', vcf_fname))
  pdf_fname <- gsub(".vcf","_vafRatio.pdf",vcf_fname)
  pdf_fname <- gsub("./","",pdf_fname)
  plot_title <- gsub("_BulkSeqSMuRF_vafRatio.pdf", "", pdf_fname)
  if (nrow(myvcf) == 0) {
    next
  }
  mydf1 <- data.frame()
  for (sname in samples(header(myvcf))) {
    ad <- sapply(geno(myvcf)$AD[,sname],"[[",2)
    dp <- geno(myvcf)$DP[,sname]
    VAF <- ad/dp
    sname <- gsub(".+-(.+)","\\1",sname)
    mydf1 <- rbind( mydf1, cbind( sname, as.numeric(VAF), names(ad)) )
  }
  # do this loop on the myvcf2 adapt de sname to _Relapse and then add to the mydf1. This way the 
  #mydf2 <- data.frame()
  for (sname in samples(header(myvcf2))) {
    ad <- sapply(geno(myvcf2)$AD[,sname],"[[",2)
    dp <- geno(myvcf2)$DP[,sname]
    VAF <- ad/dp
    sname <- gsub(".+-(.+)","\\1",sname)
    sname <- paste0("Rx_", sname)
    mydf1 <- rbind( mydf1, cbind(sname, as.numeric(VAF), names(ad)))
  }
  colnames(mydf1) <- c("bulk","VAF", "variant")
  #colnames(mydf2) <- c("bulk","VAF", "variant")
  #mydf1$bulk <- paste("Dx", mydf1$bulk, sep="_")
  #mydf2$bulk <- paste("Rx", mydf2$bulk, sep="_")

  #mydf <- rbind(mydf1, mydf2)
  
  # compare the overlapping ad names between diagnosis and relapse, en haal de teveel aan relapse er dan uit. 
  
  mydf1$VAF <- as.numeric(mydf1$VAF)
  
  splitdf<- split(mydf1, mydf1$bulk)
  # add variants to column, with NA for each variant missing. merge lists to dataframe based on this column
  #splitdf2 <- merge(mydf1, mydf2, by= "variant")
  #splitlist <- splitlist[splitlist["Dx_ALLBULK"]$Dx_ALLBULK$VAF > 0,]
  
  ratiodf <- data.frame(DNCD1aNeg = splitdf$DNCD1aNeg$VAF / splitdf$ALLBULK$VAF)
  ratiodf$DNCD1aPos <- splitdf$DNCD1aPos$VAF / splitdf$ALLBULK$VAF
  ratiodf$iSPCD4 <- splitdf$iSPCD4$VAF / splitdf$ALLBULK$VAF
  ratiodf$SPCD4 <- splitdf$SPCD4$VAF / splitdf$ALLBULK$VAF
  ratiodf$MSCBULK <- splitdf$MSCBULK$VAF / splitdf$ALLBULK$VAF
  ratiodf$ALLBULK <- splitdf$ALLBULK$VAF / median(splitdf$ALLBULK$VAF)
  # does not work since the data frame is generated based on the size of the diagnosis which has fewer variants. 
  ratiodf$Rx_ALLBULK <- splitdf$Rx_ALLBULK$VAF / median(splitdf$ALLBULK$VAF)
  ratiodf$Rx_DNCD1aPos <- splitdf$Rx_DNCD1aPos$VAF / splitdf$ALLBULK$VAF
  ratiodf$Rx_DNCD1aNeg <- splitdf$Rx_DNCD1aNeg$VAF / splitdf$ALLBULK$VAF
  ratiodf$RowNames <- rownames(ratiodf)
  if (nrow(ratiodf) == 1 ) {
    next
  }
  df.molten <- melt( ratiodf, id.vars="RowNames")
  
  if (mean(df.molten$value) == 1 ) {
    next
  }
  colnames(df.molten) <- c("variant", "bulk", "value")
  #df.molten wegschrijven en combineren met diagnose + bulk en dan plotten
  
  # assumption tests for one way Annova
  # normal distribution
  n1 <- ggqqplot(df.molten, "value", facet.by = "bulk")
  
  df.molten.filter <- df.molten %>% filter(!bulk %in% c("Dx_ALLBULK", "DNCD1aPos"))
  res.aov <- anova_test(data = df.molten.filter, dv = value, wid = variant, within = bulk) 
  get_anova_table(res.aov)
  
  # pairwise comparisons
  pwc <- df.molten %>% pairwise_t_test(value ~ bulk, paired = TRUE,
                                       p.adjust.method = "bonferroni")
  
  # Visualization: box plots with p-values
  pwc <- pwc %>% add_xy_position(x = "value")
  pwc1 <- pwc %>% filter(group2 %in% c("ALLBULK"))
  
  p <- ggplot(df.molten, aes(bulk, value, color = bulk, fill = bulk)) +
    geom_boxplot(alpha = 0.4) +
    theme_classic() +
    geom_hline(yintercept=1, linetype="dotted")+
    scale_color_manual(values=pal)+ 
    scale_fill_manual(values=pal)+
    scale_x_discrete(limits=SampleOrder)+
    ggtitle(plot_title)+
    ylab("observed/expected VAF")+
    labs(subtitle = get_test_label(res.aov, detailed = TRUE),
         caption = get_pwc_label(pwc))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p1 <- p +
    stat_pvalue_manual(pwc1, xmin = "group1", xmax = "group2", inherit.aes = F) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc1)
    )
  
  p2 <- p + 
    stat_pvalue_manual(pwc, xmin = "group1", xmax = "group2", inherit.aes = F) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc)
    )
  
  pdf(paste0("plots/ANNOVA_",pdf_fname))
  print(n1)
  print(p)
  print(p1)
  print(p2)
  dev.off()
}
#getwd()

#unique(mydf2$bulk)
#head(mydf1)
#tail(mydf1)
sname
samples(header(myvcf2))
listvcf1 <- sapply(geno(myvcf)$AD[,"pt1180-DX1PB-ALLBULK"],"[[",2)
listvcf2 <- sapply(geno(myvcf2)$AD[,"pt1909-DX1BM-ALLBULK"],"[[",2)

length(listvcf1)
length(listvcf2)

intersect(listvcf1, listvcf2)
