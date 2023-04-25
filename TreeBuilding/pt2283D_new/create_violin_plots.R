### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

library(VariantAnnotation)
library(ggplot2)


setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt2283D_new/manual_filter/branches_vcfs/overlapBulkSeq/")


#-----colors
### Colour selection 
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A")
names(pal) <- c("ALLBULK", "DN", "DNCD1aNeg","SPCD4", "DP", "MSCBULK", "MONOBulk","DNCD1aPos", "iSPCD4")
SampleOrder <- c("MSCBULK", "ALLBULK",  "DN", "iSPCD4","DP", "SPCD4")


#vcf_fname
#pdf_fname
## -- create and save plots
for (vcf_fname in list.files(path='.',pattern='.BulkSeqSMuRF.vcf',full.names = T)) {
  pdf_fname <- gsub(".vcf","_violin.pdf",vcf_fname)
  myvcf <- readVcf(vcf_fname)
  if (nrow(myvcf) == 0 ) {
    next
  }
  mydf <- data.frame()
  for (sname in samples(header(myvcf))) {
    ad <- sapply(geno(myvcf)$AD[,sname],"[[",2)
    dp <- geno(myvcf)$DP[,sname]
    vaf <- ad/dp
    sname <- gsub(".+-(.+)","\\1",sname)
    mydf <- rbind( mydf, cbind( sname, as.numeric(vaf)) )
  }
  colnames(mydf) <- c("bulk","vaf")
  
  mydf$vaf <- as.numeric(mydf$vaf)
  
  p <- ggplot( mydf, aes(x=bulk, y=vaf, fill = bulk)) +
    geom_violin(show.legend = FALSE) +
    geom_jitter(shape=16, position=position_jitter(0.2), show.legend = FALSE) +
    theme_classic() +
    scale_fill_manual(values = pal) + 
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))+ 
    scale_x_discrete(limits=SampleOrder) +
    scale_y_continuous(breaks=c(seq(0,1,0.5)), limits = c(-0.1, 1.1))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(text = element_text(size = 25))
  
  pdf(paste0("../../ViolinPlots/", pdf_fname))
  print(p)
  dev.off()
}
