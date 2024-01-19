library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FreecBAF/new")

tcr_genes <- read.table("../../../1_Input/FreecBAF/TCR_genes.bed",sep="\t")
colnames(tcr_genes) <- c("CHROM","START","END","GENE")
tcr_genes$TCR <- substr(tcr_genes$GENE,1,3)
tcr_genes <- tcr_genes[tcr_genes$TCR != "",]
tcr_genes$VDJ <- substr(tcr_genes$GENE,4,4)

MSC_ratio <- read.table("pt2229-DX1BM-MSCBULK_dedup.mpileup_ratio_awk.txt",header=T)
sample_ratio <- read.table("pt2229-DX1BM-ALLBULK_dedup.mpileup_ratio_awk.txt",header=T)

for (myTCR in "TRG") {
  tcr_gene_tmp <- tcr_genes[tcr_genes$TCR == myTCR,]
  tcr_gene_tmp_gene_order <- tcr_gene_tmp[order(tcr_gene_tmp$START),]$GENE
  tcr_gene_tmp$GENE <- factor(tcr_gene_tmp$GENE, levels=tcr_gene_tmp_gene_order)
  
  MSC_ratio_tmp <- MSC_ratio[ MSC_ratio$Chromosome == unique(tcr_gene_tmp$CHROM) & MSC_ratio$Start >= min(tcr_gene_tmp$START) & MSC_ratio$Start <= max(tcr_gene_tmp$END),]
  MSC_ratio_tmp_tcr <- merge( MSC_ratio_tmp, tcr_genes, all = T) %>% filter( Chromosome == CHROM, Start >= START, Start <= END)
  
  MSC_ratio_tmp_tcr <- MSC_ratio_tmp_tcr[MSC_ratio_tmp_tcr$TCR == myTCR,]
  MSC_ratio_tmp_tcr$GENE <- factor(MSC_ratio_tmp_tcr$GENE, levels=tcr_gene_tmp_gene_order)
  
  sample_ratio_tmp <- sample_ratio[ sample_ratio$Chromosome == unique(tcr_gene_tmp$CHROM) & sample_ratio$Start >= min(tcr_gene_tmp$START) & sample_ratio$Start <= max(tcr_gene_tmp$END),]
  sample_ratio_tmp_tcr <- merge( sample_ratio_tmp, tcr_genes, all = T) %>% filter( Chromosome == CHROM, Start >= START, Start <= END)
  
  sample_ratio_tmp_tcr <- sample_ratio_tmp_tcr[sample_ratio_tmp_tcr$TCR == myTCR,]
  sample_ratio_tmp_tcr$GENE <- factor(sample_ratio_tmp_tcr$GENE, levels=tcr_gene_tmp_gene_order)
  
  p1 <- ggplot(data=tcr_gene_tmp) +
    geom_rect(data=tcr_gene_tmp,mapping = aes(xmin=as.numeric(START),xmax=as.numeric(END),ymin=0.5,ymax=2, fill=VDJ)) +
    facet_grid(.~GENE, scales="free_x") +
    theme(legend.position = "none") +
    ylab("x") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()) +
    ylim(c(0.5,2))
  
  p2 <- ggplot(data=MSC_ratio_tmp_tcr) +
    geom_point(aes(x=Start,y=Ratio)) +
    geom_point(data=sample_ratio_tmp_tcr, aes(x=Start,y=Ratio), col="red") +
    facet_grid(.~GENE, scales="free_x") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()) +
    theme(strip.background = element_blank(),strip.text.x = element_blank())
  
  p3 <- ggplot(data=MSC_ratio_tmp_tcr) +
    geom_point(aes(x=Start,y=CopyNumber)) +
    geom_point(data=sample_ratio_tmp_tcr, aes(x=Start,y=CopyNumber), col="red") +
    facet_grid(.~GENE, scales="free_x") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()) +
    theme(strip.background = element_blank(),strip.text.x = element_blank()) +
    scale_y_continuous(limits=c(0,4), breaks=c(0,1,2,3,4),labels=function(x) sprintf("%.1f", x)) 
    
  p4 <- ggplot(data=MSC_ratio_tmp_tcr) +
    geom_point(aes(x=Start,y=estimatedBAF)) +
    geom_point(data=sample_ratio_tmp_tcr, aes(x=Start,y=estimatedBAF), col="red") +
    facet_grid(.~GENE, scales="free_x") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()) +
    theme(strip.background = element_blank(),strip.text.x = element_blank()) +
    scale_y_continuous(limits = c(0, 1), breaks=c(0,0.5,1)) 
    
  ggarrange(p1, p2, p3, p4, ncol = 1, nrow = 4)
}
  
  # baf_tcr <- merge( baf, tcr_genes, all = T) %>% filter( Chromosome == CHROM, Position >= START, Position <= END)
  # test2 <- baf_tcr[baf_tcr$TCR == "TRD",]
  # 
  # ggplot(data=test2) +
  #   geom_point(aes(x=Position,y=BAF)) +
  #   facet_grid(.~GENE, scales="free_x")
  # 
  # ratio <- read.table("pt2229-DX1BM-MSCBULK_dedup.mpileup_ratio_awk.txt",header=T)
  # 
  # ratio_tcr <- merge( ratio, tcr_genes, all = T) %>% filter( Chromosome == CHROM, Start >= START, Start <= END)
  #   