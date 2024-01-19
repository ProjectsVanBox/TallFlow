#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(rstatix)


setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt1909/branches_vcfs/overlapDiagnosis/")
#dir.create("OERatio_var", showWarnings = TRUE, recursive = FALSE, mode = "0777")
#colors
### Colour selection
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4", "#A4A4A4", "#385725", "#E7298A")
names(pal) <- c("ALLBULK", "DN", "DNCD1aNeg","SPCD4", "DP", "MSCBULK", "MONOBULK","DNCD1aPos", "iSPCD4")
#SampleOrder <- c("MONOBULK", "ALLBULK", "DN", "SPCD4","DP", "iSPCD4")
SampleOrder <- c("ALLBULK", "MSCBULK","DNCD1aNeg", "DNCD1aPos", "iSPCD4", "SPCD4")
#list1<-list.files(path='.',pattern='_branch.vcf',full.names = T)

path <- "~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt1180/branches_vcfs/overlapRelapse/"
vcf_fname = "TALL9_TALL7_branch_BulkSeqSMuRF.vcf"
for (vcf_fname in vcf_fname) {
  myvcf <- readVcf(vcf_fname)
  pdf_fname <- gsub(".vcf","_vafRatio.pdf",vcf_fname)
  pdf_fname <- gsub("./","",pdf_fname)
  plot_title <- gsub("_BulkSeqSMuRF_vafRatio.pdf", "", pdf_fname)
 if (nrow(myvcf) == 0 ) {
   next
   }
  mydf <- data.frame()
  for (sname in samples(header(myvcf))) {
    ad <- sapply(geno(myvcf)$AD[,sname],"[[",2)
    dp <- geno(myvcf)$DP[,sname]
    VAF <- ad/dp
    sname <- gsub(".+-(.+)","\\1",sname)
    mydf <- rbind( mydf, cbind( sname, as.numeric(VAF)) )
  }
  colnames(mydf) <- c("bulk","VAF")
  mydf$VAF <- as.numeric(mydf$VAF)
  
  splitdf <- as.data.frame(split(mydf, mydf$bulk))
  #splitdf <- splitdf[splitdf$ALLBULK.VAF > 0,]
  write.csv(splitdf, "Tall9_Tall9_iSPCD4overlap.csv")
  
  ratiodf <- data.frame(DNCD1aNeg = splitdf$DNCD1aNeg.VAF / splitdf$ALLBULK.VAF)
  ratiodf$DNCD1aPos <- splitdf$DNCD1aPos.VAF / splitdf$ALLBULK.VAF
  ratiodf$iSPCD4 <- splitdf$iSPCD4.VAF / splitdf$ALLBULK.VAF
  ratiodf$SPCD4 <- splitdf$SPCD4.VAF / splitdf$ALLBULK.VAF
  ratiodf$MSCBULK <- splitdf$MSCBULK.VAF / splitdf$ALLBULK.VAF
  ratiodf$ALLBULK <- splitdf$ALLBULK.VAF / median(splitdf$ALLBULK.VAF)
  ratiodf$RowNames <- rownames(ratiodf)
  if (nrow(ratiodf) == 1 ) {
    next
  }
  df.molten <- melt( ratiodf, id.vars="RowNames")
  
  if (mean(df.molten$value) == 1 ) {
    next
  }
  colnames(df.molten) <- c("variant", "bulk", "value")
  
  # assumption tests for one way Annova
  # normal distribution
  n1 <- ggqqplot(df.molten, "value", facet.by = "bulk")
  
  df.molten.filter <- df.molten %>% filter(!bulk %in% c("ALLBULK", "MSCBULK"))
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
getwd()

myvcf <- readVcf("../../../pt1180/branches_vcfs/overlapRelapse/TALL9_TALL7_branch_BulkSeqSMuRF.vcf")

myvcf <- readVcf("/Users/verapoort/hpc/projects/TallFlow/3_Output/TreeBuilding/pt1180/branches_vcfs/overlapBulkSeq/TALL9_TALL7_branch_BulkSeqSMuRF.vcf")
