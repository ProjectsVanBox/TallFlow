library(VariantAnnotation)
library(UpSetR)

setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding_indel/pt2283D/")

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
nrow(myvcf_gt_df)
upset_matrix_raw <- matrix(as.numeric(unlist(myvcf_gt_df)),ncol=ncol(myvcf_gt_df))
rownames(upset_matrix_raw) <- rownames(myvcf_gt_df)
colnames(upset_matrix_raw) <- colnames(myvcf_gt_df)
upset_matrix <- upset_matrix_raw
upset_matrix <- upset_matrix[which(rowSums(upset_matrix) > 1),]
upset_matrix <- upset_matrix[which(rowSums(upset_matrix) < 12),]
pdf("UpsetPlot.pdf")
upset(as.data.frame(upset_matrix),nsets=ncol(upset_matrix), order.by = "freq")
dev.off()
upset_matrix2 <- cbind(apply(upset_matrix,1,paste, collapse=""), upset_matrix)
upset_matrix2_table <- table(upset_matrix2[,1])

names_to_check <- names(upset_matrix2_table[which(upset_matrix2_table < 50 & upset_matrix2_table > 5)])
for ( x in names_to_check ) {
  upset_matrix_to_check <- upset_matrix2[upset_matrix2[,1] %in% x,]
  label <- colnames(upset_matrix_to_check)[upset_matrix_to_check[1,] == "1"]
  label <- gsub("pt2283D-DX1BM-","",label)
  label <- paste(label,collapse = "_")
  write.table(upset_matrix_to_check[,-1],paste(label,"_upset_matrix_to_check.txt",sep=""),quote=F,sep=" ")
  
  myvcf_to_check <- myvcf[rownames(myvcf) %in% rownames(upset_matrix_to_check),]
  writeVcf( myvcf_to_check, paste(label,"_to_check.vcf",sep=""))
}

upset_df <- as.data.frame(upset_matrix)





shared_Bulk <- upset_df[ 
  upset_df$`pt2283D-DX1BM-ALLBULK` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2", "pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10",
                         "pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 5
  , ]
pdf("shared_BULK_upset.pdf")
upset(shared_Bulk,nsets=ncol(shared_Bulk), order.by = "freq")
dev.off()



shared_1 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL1` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10",
                         "pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 5
  , ]
pdf("shared_TALL1_upset.pdf")
upset(shared_1,nsets=ncol(shared_1), order.by = "freq")
dev.off()

shared_2 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL2` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10","pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 4
  , ]
pdf("shared_TALL2_upset.pdf")
upset(shared_2,nsets=ncol(shared_2), order.by = "freq")
dev.off()

shared_3 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL3` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10","pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 5
  , ]
pdf("shared_TALL3_upset.pdf")
upset(shared_3,nsets=ncol(shared_3), order.by = "freq")
dev.off()

shared_4 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL4` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10","pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 4
  , ]
pdf("shared_TALL4_upset.pdf")
upset(shared_4,nsets=ncol(shared_4), order.by = "freq")
dev.off()

shared_5 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL5` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10","pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 4
  , ]
pdf("shared_TALL5_upset.pdf")
upset(shared_5,nsets=ncol(shared_5), order.by = "freq")
dev.off()

shared_6 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL6` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10", "pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 4
  , ]
pdf("shared_TALL6_upset.pdf")
upset(shared_6,nsets=ncol(shared_6), order.by = "freq")
dev.off()

shared_7 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL7` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10"),"pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12"]  ) < 5
  , ]
pdf("shared_TALL7_upset.pdf")
upset(shared_7,nsets=ncol(shared_7), order.by = "freq")
dev.off()

shared_8 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL8` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL10","pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 5
  , ]
pdf("shared_TALL8_upset.pdf")
upset(shared_8,nsets=ncol(shared_8), order.by = "freq")
dev.off()

shared_9 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL9` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL10","pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 4
  , ]
pdf("shared_TALL9_upset.pdf")
upset(shared_9,nsets=ncol(shared_9), order.by = "freq")
dev.off()

shared_10 <- upset_df[ 
  upset_df$`pt2283D-DX1BM-TALL10` == 1 &
    rowSums( upset_df[,c("pt2283D-DX1BM-TALL2","pt2283D-DX1BM-TALL4","pt2283D-DX1BM-TALL5","pt2283D-DX1BM-TALL6","pt2283D-DX1BM-TALL9","pt2283D-DX1BM-TALL11", "pt2283D-DX1BM-TALL12")]  ) < 4
  , ]
pdf("shared_TALL10_upset.pdf")
upset(shared_10,nsets=ncol(shared_10), order.by = "freq")
dev.off()

