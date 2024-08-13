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


setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt3291/branches_vcfs/overlapBulkSeq/")
#dir.create("OERatio_var", showWarnings = TRUE, recursive = FALSE, mode = "0777")
#colors
### Colour selection
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4", "#A4A4A4", "#385725", "#E7298A", "#E28028")
names(pal) <- c("ALLBULK", "DN", "DNCD1aNeg","SPCD4", "DP", "MSCBULK", "MONOBulk","DNCD1aPos", "iSPCD4", "SPCD8")
#SampleOrder <- c("MONOBULK", "ALLBULK", "DN", "SPCD4","DP", "iSPCD4")
SampleOrder <- c("ALLBULK", "MONOBULK", "DNCD1aNeg", "DNCD1aPos", "DP", "SPCD8")
#list1<-list.files(path='.',pattern='_branch.vcf',full.names = T)


for (vcf_fname in list.files(path='.',pattern='_BulkSeqSMuRF.vcf',full.names = T)) {
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
  splitdf <- splitdf[splitdf$ALLBULK.VAF > 0,]
  
  ratiodf <- data.frame(DNCD1aPos = splitdf$plus.VAF / splitdf$ALLBULK.VAF)
  ratiodf$DNCD1aNeg <- splitdf$min.VAF / splitdf$ALLBULK.VAF
  ratiodf$DP <- splitdf$DP.VAF / splitdf$ALLBULK.VAF
  ratiodf$SPCD8 <- splitdf$SPCD8.VAF / splitdf$ALLBULK.VAF
  ratiodf$MONOBULK <- splitdf$MONOBULK.VAF / splitdf$ALLBULK.VAF
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
  
  df.molten.filter <- df.molten %>% filter(bulk %in% c("DNCD1aNeg", "DNCD1aPos", "DP", "SPCD8"))
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
  
  write.csv(df.molten, paste0("plots/", plot_title, "_df.csv"))
}


