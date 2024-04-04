
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(ggrepel)

setwd("~/hpc/projects/TallFlow/3_Output/FreecBAF/pt2229/output")

tcr_regions <- read.table("../../../../1_Input/FreecBAF/TCR_subset.bed")
tcr_regions <- tcr_regions[-c(1:3),]
tcr_regions[tcr_regions$V4 == "TRA",c(2:3)] <- c(21621904,22552132)
tcr_regions[tcr_regions$V4 == "TRD",c(2:3)] <- c(22891537,22935569)
tcr_regions[tcr_regions$V4 == "TRG",c(2:3)] <- c(38240024,38368055)
tcr_regions[tcr_regions$V4 == "TRB",c(2:3)] <- c(142299011,142813287)

tcr_genes <- read.table("../../../../1_Input/FreecBAF/TCR_genes.bed",sep="\t")
colnames(tcr_genes) <- c("CHROM","START","END","GENE")
tcr_genes$TCR <- substr(tcr_genes$GENE,1,3)
tcr_genes <- tcr_genes[tcr_genes$TCR != "",]
tcr_genes$VDJ <- substr(tcr_genes$GENE,4,4)

MSC_BAF <- read.table("../output/pt2229-DX1BM-MSCBULK_dedup.mpileup_BAF.txt",sep="\t",header=T)
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
  
  for (myTCR in unique(tcr_regions$V4)) {
    tcr_region_tmp <- tcr_regions[tcr_regions$V4 == myTCR,]
    
    MSC_BAF_ratio_tcr <- MSC_BAF_ratio[ MSC_BAF_ratio$Chromosome == unique(tcr_region_tmp$V1) & MSC_BAF_ratio$Position >= min(tcr_region_tmp$V2) & MSC_BAF_ratio$Position <= max(tcr_region_tmp$V3),]
    MSC_BAF_ratio_tcr$TCR <- myTCR 
    
    sample_BAF_ratio_tcr <- sample_BAF_ratio[ sample_BAF_ratio$Chromosome == unique(tcr_region_tmp$V1) & sample_BAF_ratio$Position >= min(tcr_region_tmp$V2) & sample_BAF_ratio$Position <= max(tcr_region_tmp$V3),]
    sample_BAF_ratio_tcr$TCR <- myTCR
    
    dummy <- as.data.frame(MSC_BAF_ratio_tcr) %>% group_by(TCR) %>%
        summarise_at(vars(Position),list(min = min)) %>%
        slice(rep(1:n(), each = 6))
    
    colnames(dummy) <- c("GENE","Position")
    dummy$variable <- rep(c("BAF","Ratio","CopyNumber"),nrow(dummy)/3)
    dummy$value <- rep(c(0,0,0,1,2,4),nrow(dummy)/6)
  
    
    p1 <- ggplot(data=MSC_BAF_ratio_tcr,aes(x=Position,y=value)) +
      stat_binhex(bins = nrow (MSC_BAF_ratio_tcr) / 250, col = "black", show.legend = FALSE)+
      geom_blank(data=dummy) +
      stat_binhex(data=sample_BAF_ratio_tcr,aes(x=Position,y=value),col="red", bins = nrow (sample_BAF_ratio_tcr) / 250, show.legend = FALSE) +
      facet_grid( variable ~ ., scales = "free") 
    
    #pdf(paste0(sname,"_",myTCR,".pdf"),width=30,height=15)
    print(p1)
    #dev.off()
    
    tcr_region2 <- tcr_genes[tcr_genes$TCR == myTCR,] %>% group_by( VDJ ) %>% summarise_at(vars(START, END),list(min = min, max = max))
    tcr_region2$TCR <- myTCR
    
    p2 <- ggplot() +
      geom_rect(data=tcr_region2,mapping = aes(xmin=as.numeric(START_min),xmax=as.numeric(END_max),ymin=0.5,ymax=2, fill=VDJ)) +
      ylab("x") +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()) +
      theme(legend.position = "none") +
      facet_grid( TCR ~ . ) +
      scale_y_continuous(limits=c(0,2), labels=function(x) sprintf("%.2f", x)) 
    
    tcr_region3 <- tcr_genes[tcr_genes$TCR == myTCR,]
      
    p3 <- ggplot() +
      #geom_rect(data=tcr_region3,mapping = aes(xmin=as.numeric(START),xmax=as.numeric(END),ymin=0.5,ymax=2, fill=VDJ)) +
      ylab("x") +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()) +
      theme(legend.position = "top") +
      facet_grid( TCR ~ . ) +
      scale_y_continuous(limits=c(0,2), labels=function(x) sprintf("%.2f", x)) +
      geom_text_repel(data=tcr_region3, mapping=aes(x=((START+END)/2), y=0, label=GENE, col=VDJ), 
                size=2, force_pull   = 0, max.overlaps = Inf,
                nudge_y      = 2,
                direction    = "x",
                angle        = 90,
                hjust        = 0,
                segment.size = 0.2,  )
    #,position=position_jitter(width=0,height=2))
    
    pdf(paste0(sname,"_",myTCR,"_3.pdf"))
    grid.arrange(p3, p2, p1, ncol= 1, nrow = 3, heights= c(3,1,4))
    dev.off()
    
  }
}
