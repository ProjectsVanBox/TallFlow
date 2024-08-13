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


setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt2322/withoutBulk/branches_vcfs/overlapBulkSeq/")
dir.create("OERatio_var", showWarnings = TRUE, recursive = FALSE, mode = "0777")
#colors
### Colour selection
pal = c("#FFA500","#33A02C","#33A02C", "#E7298A", "#B17BA6", "#A4A4A4", "#A4A4A4", "#385725", "#1965B0")
names(pal) <- c("ALLBULK", "DN", "DNCD1aNeg","SPCD4", "DP", "MSCBULK", "MONOBulk","DNCD1aPos", "iSPCD4")
SampleOrder <- c("ALLBULK", "MONOBULK","DN", "SPCD4","DP", "iSPCD4")
#SampleOrder <- c("ALLBULK", "MSCBULK","DNCD1aNeg", "DNCD1aPos", "iSPCD4", "DP")

# group comparisons
my_comparisons <- list( c("ALLBULK", "MONOBULK"), c("ALLBULK", "DN"), c("ALLBULK", "SPCD4"), 
                        c("ALLBULK", "DP"), c("ALLBULK", "iSPCD4"))

#samples(header(myvcf))
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
  
  ratiodf <- data.frame(DN = splitdf$DN.VAF / splitdf$ALLBULK.VAF))
  ratiodf$SPCD4 <- splitdf$SPCD4.VAF / splitdf$ALLBULK.VAF
  ratiodf$DP <- splitdf$DP.VAF / splitdf$ALLBULK.VAF
  ratiodf$iSPCD4 <- splitdf$iSPCD4.VAF / splitdf$ALLBULK.VAF
  ratiodf$MONOBULK <- splitdf$MONOBULK.VAF / splitdf$ALLBULK.VAF
  ratiodf$ALLBULK <- splitdf$ALLBULK.VAF / median(splitdf$ALLBULK.VAF)
  ratiodf$RowNames <- rownames(ratiodf)
  
  df.molten <- melt( ratiodf, id.vars="RowNames")
  df.molten <- melt( ratiodf, id.vars="RowNames")
  p <- ggplot(df.molten, aes(variable, value, color = variable, fill = variable)) +
    geom_boxplot(alpha = 0.4) +
    stat_summary(fun = "median", fun.min = "median", fun.max= "median", size= 0.3, geom = "crossbar") +
    theme_classic() +
    geom_hline(yintercept=1, linetype="dotted")+
    scale_color_manual(values=pal)+ 
    scale_fill_manual(values=pal)+
    scale_x_discrete(limits=SampleOrder)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 3)+
    ggtitle(plot_title)+
    ylab("observed/expected VAF")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  pdf(paste0("OERatio_var/",pdf_fname))
  print(p)
  dev.off()
}

getwd()

myvcf <- readVcf("root_branch_BulkSeqSMuRF.vcf")

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
  
ratiodf <- data.frame(DN = splitdf$DN.VAF / splitdf$ALLBULK.VAF)
ratiodf$SPCD4 <- splitdf$SPCD4.VAF / splitdf$ALLBULK.VAF
ratiodf$DP <- splitdf$DP.VAF / splitdf$ALLBULK.VAF
ratiodf$iSPCD4 <- splitdf$iSPCD4.VAF / splitdf$ALLBULK.VAF
ratiodf$MONOBULK <- splitdf$MONOBULK.VAF / splitdf$ALLBULK.VAF
ratiodf$ALLBULK <- splitdf$ALLBULK.VAF / median(splitdf$ALLBULK.VAF)
ratiodf$RowNames <- rownames(ratiodf)

df.molten <- melt( ratiodf, id.vars="RowNames")
df.molten <- melt( ratiodf, id.vars="RowNames")
p <- ggplot(df.molten, aes(variable, value, color = variable, fill = variable)) +
  geom_boxplot(alpha = 0.4) +
  stat_summary(fun = "median", fun.min = "median", fun.max= "median", linewidth= 0.3, geom = "crossbar") +
  theme_classic() +
  geom_hline(yintercept=1, linetype="dotted")+
  scale_color_manual(values=pal)+ 
  scale_fill_manual(values=pal)+
  scale_x_discrete(limits=SampleOrder)+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 3)+
  ggtitle("plot_title")+
  ylab("observed/expected VAF")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(paste0("OERatio/",pdf_fname))
print(p)
dev.off()  
  
  


head(splitdf)
mydf
ratiodf
dev.off()

print("hello")
