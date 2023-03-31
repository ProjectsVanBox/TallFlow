library(VariantAnnotation)
library(UpSetR)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2229/")

myvcf <- readVcf("SMuRF.filtered.coords.vcf")
rownames(myvcf) <- paste(as.character(seqnames(myvcf)),start(myvcf),sep=":")

myvcf_gt_df <- as.data.frame(geno(myvcf)$GT)  
myvcf_gt_df[is.na(myvcf_gt_df)] <- 0
myvcf_gt_df[myvcf_gt_df == '0/0'] <- 0
myvcf_gt_df[myvcf_gt_df == '0|0'] <- 0
myvcf_gt_df[myvcf_gt_df == './.'] <- 0
myvcf_gt_df[myvcf_gt_df == '.|.'] <- 0
myvcf_gt_df[myvcf_gt_df == '0/1'] <- 1
myvcf_gt_df[myvcf_gt_df == '0|1'] <- 1
myvcf_gt_df[myvcf_gt_df == '1/1'] <- 1
myvcf_gt_df[myvcf_gt_df == '1|1'] <- 1

upset_matrix_raw <- matrix(as.numeric(unlist(myvcf_gt_df)),ncol=ncol(myvcf_gt_df))
rownames(upset_matrix_raw) <- rownames(myvcf_gt_df)
colnames(upset_matrix_raw) <- colnames(myvcf_gt_df)
upset_matrix <- upset_matrix_raw


upset_df <- as.data.frame(upset_matrix)
branch_6_1_2 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL6` == 1 & 
  upset_df$`pt2229-DX1BM-TALL2` == 1 & 
  upset_df$`pt2229-DX1BM-TALL7` == 0 & 
  upset_df$`pt2229-DX1BM-TALL10` == 0 & 
  upset_df$`pt2229-DX1BM-TALL8` == 0 & 
  upset_df$`pt2229-DX1BM-TALL9` == 0 & 
  upset_df$`pt2229-DX1BM-TALL5` == 0 & 
  upset_df$`pt2229-DX1BM-TALL3` == 0 & 
  upset_df$`pt2229-DX1BM-TALL4` == 0
  , ]

branch_1_2 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL1` == 1 & 
  upset_df$`pt2229-DX1BM-TALL2` == 1 & 
  upset_df$`pt2229-DX1BM-TALL6` == 0 & 
  upset_df$`pt2229-DX1BM-TALL7` == 0 & 
  upset_df$`pt2229-DX1BM-TALL10` == 0 & 
  upset_df$`pt2229-DX1BM-TALL8` == 0 & 
  upset_df$`pt2229-DX1BM-TALL9` == 0 & 
  upset_df$`pt2229-DX1BM-TALL5` == 0 & 
  upset_df$`pt2229-DX1BM-TALL3` == 0 & 
  upset_df$`pt2229-DX1BM-TALL4` == 0
  , ]

branch_7_10_8_9_5_3_4 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL1` == 0 & 
    upset_df$`pt2229-DX1BM-TALL2` == 0 & 
    upset_df$`pt2229-DX1BM-TALL6` == 0 & 
    upset_df$`pt2229-DX1BM-TALL10` == 1 & 
    upset_df$`pt2229-DX1BM-TALL9` == 1 & 
    upset_df$`pt2229-DX1BM-TALL5` == 1 & 
    upset_df$`pt2229-DX1BM-TALL4` == 1
  , ]

branch_10_8_9 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL1` == 0 & 
    upset_df$`pt2229-DX1BM-TALL2` == 0 & 
    upset_df$`pt2229-DX1BM-TALL6` == 0 & 
    upset_df$`pt2229-DX1BM-TALL7` == 0 & 
    upset_df$`pt2229-DX1BM-TALL10` == 1 & 
    upset_df$`pt2229-DX1BM-TALL9` == 1 & 
    upset_df$`pt2229-DX1BM-TALL5` == 0 & 
    upset_df$`pt2229-DX1BM-TALL3` == 0 & 
    upset_df$`pt2229-DX1BM-TALL4` == 0
  , ]

branch_8_9 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL1` == 0 & 
    upset_df$`pt2229-DX1BM-TALL2` == 0 & 
    upset_df$`pt2229-DX1BM-TALL6` == 0 & 
    upset_df$`pt2229-DX1BM-TALL7` == 0 & 
    upset_df$`pt2229-DX1BM-TALL10` == 0 & 
    upset_df$`pt2229-DX1BM-TALL9` == 1 & 
    upset_df$`pt2229-DX1BM-TALL8` == 1 & 
    upset_df$`pt2229-DX1BM-TALL5` == 0 & 
    upset_df$`pt2229-DX1BM-TALL3` == 0 & 
    upset_df$`pt2229-DX1BM-TALL4` == 0
  , ]


branch_5_3_4 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL1` == 0 & 
    upset_df$`pt2229-DX1BM-TALL2` == 0 & 
    upset_df$`pt2229-DX1BM-TALL6` == 0 & 
    upset_df$`pt2229-DX1BM-TALL7` == 0 & 
    upset_df$`pt2229-DX1BM-TALL10` == 0 & 
    upset_df$`pt2229-DX1BM-TALL9` == 0 & 
    upset_df$`pt2229-DX1BM-TALL8` == 0 & 
    upset_df$`pt2229-DX1BM-TALL5` == 1 & 
    upset_df$`pt2229-DX1BM-TALL3` == 0 & 
    upset_df$`pt2229-DX1BM-TALL4` == 1
  , ]


branch_3_4 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL1` == 0 & 
    upset_df$`pt2229-DX1BM-TALL2` == 0 & 
    upset_df$`pt2229-DX1BM-TALL6` == 0 & 
    upset_df$`pt2229-DX1BM-TALL7` == 0 & 
    upset_df$`pt2229-DX1BM-TALL10` == 0 & 
    upset_df$`pt2229-DX1BM-TALL9` == 0 & 
    upset_df$`pt2229-DX1BM-TALL8` == 0 & 
    upset_df$`pt2229-DX1BM-TALL5` == 0 & 
    upset_df$`pt2229-DX1BM-TALL3` == 1 & 
    upset_df$`pt2229-DX1BM-TALL4` == 1
  , ]

shared_1 <- upset_df[ 
    upset_df$`pt2229-DX1BM-TALL1` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 5
  , ]
pdf("shared_TALL1_upset.pdf")
upset(shared_1,nsets=ncol(shared_1), order.by = "freq")
dev.off()

shared_2 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL2` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 4
  , ]
pdf("shared_TALL2_upset.pdf")
upset(shared_2,nsets=ncol(shared_2), order.by = "freq")
dev.off()

shared_3 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL3` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 5
  , ]
pdf("shared_TALL3_upset.pdf")
upset(shared_3,nsets=ncol(shared_3), order.by = "freq")
dev.off()

shared_4 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL4` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 4
  , ]
pdf("shared_TALL4_upset.pdf")
upset(shared_4,nsets=ncol(shared_4), order.by = "freq")
dev.off()

shared_5 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL5` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 4
  , ]
pdf("shared_TALL5_upset.pdf")
upset(shared_5,nsets=ncol(shared_5), order.by = "freq")
dev.off()

shared_6 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL6` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 4
  , ]
pdf("shared_TALL6_upset.pdf")
upset(shared_6,nsets=ncol(shared_6), order.by = "freq")
dev.off()

shared_7 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL7` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 5
  , ]
pdf("shared_TALL7_upset.pdf")
upset(shared_7,nsets=ncol(shared_7), order.by = "freq")
dev.off()

shared_8 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL8` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9","pt2229-DX1BM-TALL10")]  ) < 5
  , ]
pdf("shared_TALL8_upset.pdf")
upset(shared_8,nsets=ncol(shared_8), order.by = "freq")
dev.off()

shared_9 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL9` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL10")]  ) < 4
  , ]
pdf("shared_TALL9_upset.pdf")
upset(shared_9,nsets=ncol(shared_9), order.by = "freq")
dev.off()

shared_10 <- upset_df[ 
  upset_df$`pt2229-DX1BM-TALL10` == 1 &
    rowSums( upset_df[,c("pt2229-DX1BM-TALL2","pt2229-DX1BM-TALL4","pt2229-DX1BM-TALL5","pt2229-DX1BM-TALL6","pt2229-DX1BM-TALL9")]  ) < 4
  , ]
pdf("shared_TALL10_upset.pdf")
upset(shared_10,nsets=ncol(shared_10), order.by = "freq")
dev.off()
