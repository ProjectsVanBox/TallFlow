
library(ggplot2)
library(reshape2)
library(dplyr)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FreecBAF/new")

tcr_genes <- read.table("../../../1_Input/FreecBAF/TCR_genes.bed",sep="\t")
colnames(tcr_genes) <- c("CHROM","START","END","GENE")
tcr_genes$TCR <- substr(tcr_genes$GENE,1,3)
tcr_genes <- tcr_genes[tcr_genes$TCR != "",]
tcr_genes$VDJ <- substr(tcr_genes$GENE,4,4)

MSC_BAF <- read.table("pt2229-DX1BM-MSCBULK_dedup.mpileup_BAF.txt",sep="\t",header=T)
MSC_ratio <- read.table("pt2229-DX1BM-MSCBULK_dedup.mpileup_ratio_subset.bed",sep="\t")
colnames(MSC_ratio) <- c("Chromosome","Start","Position","Ratio","MedianRatio","CopyNumber",	"BAFratio","estimatedBAF","Genotype","UncertaintyOfGT")
MSC_BAF_melted <- melt(MSC_BAF[,1:3],id.vars = c("Chromosome","Position"), )
MSC_ratio_melted <- melt(MSC_ratio[,c(1,3,4,6)],id.vars = c("Chromosome","Position"), )
MSC_BAF_ratio <- rbind(MSC_BAF_melted, MSC_ratio_melted)

for ( myfile in list.files(".","_BAF.txt")) {
  my_baf_file <- myfile
  my_ratio_file <- gsub("_BAF.txt","_ratio_subset.bed",myfile)
  
  if ( !file.exists(my_baf_file) || !file.exists(my_ratio_file) ) {
    next
  }
  
  sname <- gsub(".+DX1BM-(.+)_dedup.+","\\1",my_baf_file)
  if ( grepl("MSC",sname) ) {
    next
  }

  sample_BAF <- read.table(my_baf_file,sep="\t",header=T)
  sample_ratio <- read.table(my_ratio_file,sep="\t")
  colnames(sample_ratio) <- c("Chromosome","Start","Position","Ratio","MedianRatio","CopyNumber",	"BAFratio","estimatedBAF","Genotype","UncertaintyOfGT")
  sample_BAF_melted <- melt(sample_BAF[,1:3],id.vars = c("Chromosome","Position"), )
  sample_ratio_melted <- melt(sample_ratio[,c(1,3,4,6)],id.vars = c("Chromosome","Position"), )
  sample_BAF_ratio <- rbind(sample_BAF_melted, sample_ratio_melted)
  
  for (myTCR in unique(tcr_genes$TCR)) {
    tcr_gene_tmp <- tcr_genes[tcr_genes$TCR == myTCR,]
    tcr_gene_tmp_gene_order <- tcr_gene_tmp[order(tcr_gene_tmp$START),]$GENE
    tcr_gene_tmp$GENE <- factor(tcr_gene_tmp$GENE, levels=tcr_gene_tmp_gene_order)
    
    MSC_BAF_ratio_tcr <- MSC_BAF_ratio[ MSC_BAF_ratio$Chromosome == unique(tcr_gene_tmp$CHROM) & MSC_BAF_ratio$Position >= min(tcr_gene_tmp$START) & MSC_BAF_ratio$Position <= max(tcr_gene_tmp$END),]
    MSC_BAF_ratio_tcr_info <- merge( MSC_BAF_ratio_tcr, tcr_genes, all = T) %>% filter( Chromosome == CHROM, Position >= START, Position <= END)
    MSC_BAF_ratio_tcr_info$GENE <- factor(MSC_BAF_ratio_tcr_info$GENE, levels=tcr_gene_tmp_gene_order)
    
    sample_BAF_ratio_tcr <- sample_BAF_ratio[ sample_BAF_ratio$Chromosome == unique(tcr_gene_tmp$CHROM) & sample_BAF_ratio$Position >= min(tcr_gene_tmp$START) & sample_BAF_ratio$Position <= max(tcr_gene_tmp$END),]
    sample_BAF_ratio_tcr_info <- merge( sample_BAF_ratio_tcr, tcr_genes, all = T) %>% filter( Chromosome == CHROM, Position >= START, Position <= END)
    sample_BAF_ratio_tcr_info$GENE <- factor(sample_BAF_ratio_tcr_info$GENE, levels=tcr_gene_tmp_gene_order)
    
    dummy <- as.data.frame(MSC_BAF_ratio_tcr_info) %>% group_by(GENE) %>%
        summarise_at(vars(Position),list(min = min)) %>%
        slice(rep(1:n(), each = 6))
    
    colnames(dummy) <- c("GENE","Position")
    dummy$variable <- rep(c("BAF","Ratio","CopyNumber"),nrow(dummy)/3)
    dummy$value <- rep(c(0,0,0,1,2,4),nrow(dummy)/6)
  
    
    p1 <- ggplot(data=MSC_BAF_ratio_tcr_info,aes(x=Position,y=value)) +
      geom_point(alpha=0.4) +
      geom_blank(data=dummy) +
      geom_point(data=sample_BAF_ratio_tcr_info,aes(x=Position,y=value),col="red",alpha=0.4) +
      facet_grid( variable ~ GENE, scales = "free")
    
    pdf(paste0(sname,"_",myTCR,".pdf"),width=30,height=15)
    print(p1)
    dev.off()

  }
}
