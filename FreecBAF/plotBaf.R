
library(ggplot2)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FreecBAF/new")
tcr_regions <- read.table("../../1_Input/FreecBAF/TCR_subset.bed")
tcr_regions <- tcr_regions[-c(1:3),]
tcr_regions[tcr_regions$V4 == "TRA",c(2:3)] <- c(21621904,22552132)
tcr_regions[tcr_regions$V4 == "TRD",c(2:3)] <- c(22891537,22935569)
tcr_regions[tcr_regions$V4 == "TRG",c(2:3)] <- c(38240024,38368055)
tcr_regions[tcr_regions$V4 == "TRB",c(2:3)] <- c(142299011,142813287)

baf <- read.table("pt2229-DX1BM-MSCBULK_dedup.mpileup_BAF.txt",header=T)

tcr_baf <- matrix(ncol=10)[-1,]

for ( i in c(1:nrow(tcr_regions)) ) {
  tcr <- tcr_regions[i,]
  print(tcr)
  tcr_baf_tmp <- baf[ baf$Chromosome == tcr$V1 & baf$Position >= tcr$V2 & baf$Position <= tcr$V3, ]
  tcr_baf_tmp$TCR <- tcr$V4
  tcr_baf_tmp$Sample <- ".MSC"
  tcr_baf <- rbind(tcr_baf, tcr_baf_tmp)
}

baf_files <- list.files(path=".",pattern="*BAF.txt")

for ( baf_file in baf_files ) {
  sname <- gsub(".+DX1BM-(.+)_dedup.+","\\1",baf_file)
  if ( grepl("MSC",sname) ) {
    next
  }
  
  baf <- read.table(baf_file,header=T)
  
  for ( i in c(1:nrow(tcr_regions)) ) {
    tcr <- tcr_regions[i,]
    print(tcr)
    tcr_baf_tmp <- baf[ baf$Chromosome == tcr$V1 & baf$Position >= tcr$V2 & baf$Position <= tcr$V3, ]
    if (nrow(tcr_baf_tmp) > 0) {
      tcr_baf_tmp$TCR <- tcr$V4
      tcr_baf_tmp$Sample <- sname
      tcr_baf <- rbind(tcr_baf, tcr_baf_tmp)
    }
  }
  
  pdf(paste0(sname,"_BAF.pdf"),width=10,height=5)
  p <- ggplot(data=tcr_baf, aes(x=Position,y=BAF, col=Sample)) +
    geom_point(size=1) +
    facet_grid( . ~ TCR, scales = "free_x")
  print( p )
  dev.off()
  
  tcr_baf <- tcr_baf[ tcr_baf$Sample == ".MSC",]
}

# pdf("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FreecBAF/bafplots.pdf",width=10,height=20)
# p <- ggplot(data=tcr_baf, aes(x=Position,y=BAF,col=Sample )) +
#   geom_point(size=0.5) +
#   facet_grid( Y ~ TCR, scales = "free_x")
# print( p )
# dev.off()
