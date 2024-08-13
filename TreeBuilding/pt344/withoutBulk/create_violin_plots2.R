#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("VariantAnnotation")

library(VariantAnnotation)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(reshape2)

setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt344/branches_vcfs/overlapBulkSeq/")
dir.create("../../Violin", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#colors
### Colour selection 
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A", "#002eb3")
names(pal) <- c("ALLBULK", "DN", "DNCD1aNeg","SPCD4", "DP", "MSCBULK", "MONOBULK","DNCD1aPos", "iSPCD4", "SPCD8")
#SampleOrder <- c("MONOBULK", "ALLBULK",  "DN", "SPCD4","DP", "iSPCD4")
#SampleOrder <- c("MSCBULK", "ALLBULK", "DNCD1aNeg", "DNCD1aPos", "iSPCD4", "DP")
#SampleOrder <- c("MONOBULK", "ALLBULK",  "DN", "iSPCD4","DP", "SPCD4")
SampleOrder <- c("MONOBULK", "ALLBULK",  "DNCD1aNeg", "DNCD1aPos", "DP", "SPCD8")



#list1<-list.files(path='.',pattern='_branch.vcf',full.names = T)

#myvcf <- readVcf("root_branch_BulkSeqSMuRF.vcf")

for (vcf_fname in list.files(path='.',pattern='_BulkSeqSMuRF.vcf',full.names = T)) {
  pdf_fname <- gsub(".vcf","_violin_colors_ANNOVA.pdf",vcf_fname)
  pdf_fname <- gsub("./","",pdf_fname)
  plot_title <- gsub("_branch_BulkSeqSMuRF_violin_colors_ANNOVA.pdf", "", pdf_fname)
  myvcf <- readVcf(vcf_fname)
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
  splitdf$RowNames <- rownames(splitdf)
  splitdf <- splitdf %>% select(!matches("bulk"))
  df.molten <- melt(splitdf, id.vars="RowNames")
  colnames(df.molten) <- c("variant", "bulk", "value")
  df.molten$bulk <- gsub(".VAF", "", df.molten$bulk)
  
  if (nrow(df.molten) < 6 ) {
    print("only 1 variant")
    print(pdf_fname)
    next
  }

  # statistics
  n1 <- ggqqplot(df.molten, "value", facet.by = "bulk")
  
  df.molten.filter <- df.molten %>% filter(bulk %in% c("DNCD1aNeg", "DNCD1aPos", "DP", "SPCD8"))
  res.aov <- anova_test(data = df.molten.filter, dv = value, wid = variant, within = bulk) 
  get_anova_table(res.aov)
  
  # pairwise comparisons
  pwc <- df.molten %>% pairwise_t_test(value ~ bulk, paired = TRUE,
                                       p.adjust.method = "bonferroni")
  
  # Visualization: box plots with p-values
  pwc <- pwc %>% add_xy_position(x = "value")
  
  p <- ggplot( mydf, aes(x=bulk, y=VAF, fill = bulk)) +
    geom_violin(show.legend = FALSE, alpha = 0.4) +
    geom_jitter(shape=16, position=position_jitter(0.2), show.legend = FALSE) +
    theme_classic() +
    scale_fill_manual(values = pal) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))+ 
    scale_x_discrete(limits=SampleOrder) +
    scale_y_continuous(breaks=c(seq(0,1,0.5)), limits = c(-0.1, 1.4))+
    ggtitle(plot_title)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
    #theme(text = element_text(size = 25))

  p1 <- p + 
    stat_pvalue_manual(pwc, xmin = "group1", xmax = "group2", inherit.aes = F) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc)
    )
  
  pdf(paste0("../../Violin/", pdf_fname))
  print(n1)
  print(p)
  print(p1)
  dev.off()
}

#getwd()

