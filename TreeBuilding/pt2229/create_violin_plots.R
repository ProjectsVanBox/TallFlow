
library(VariantAnnotation)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/branches_vcfs/overlapBulkSeq/")

for (vcf_fname in list.files(path='.',pattern='.BulkSeqSMuRF.vcf',full.names = T)) {
  pdf_fname <- gsub(".vcf","_violin2.pdf",vcf_fname)
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
  
  p <- ggplot( mydf, aes(x=bulk, y=vaf)) +
    geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2))
  
  pdf(pdf_fname)
  print(p)
  dev.off()
}