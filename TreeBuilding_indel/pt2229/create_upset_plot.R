library(VariantAnnotation)
library(UpSetR)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding_indel/pt2322/")

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
upset_matrix <- upset_matrix[which(rowSums(upset_matrix) > 1),]
upset_matrix <- upset_matrix[which(rowSums(upset_matrix) < 11),]
pdf("UpsetPlot.pdf")
upset(as.data.frame(upset_matrix),nsets=ncol(upset_matrix), order.by = "freq")
dev.off()
upset_matrix2 <- cbind(apply(upset_matrix,1,paste, collapse=""), upset_matrix)
upset_matrix2_table <- table(upset_matrix2[,1])

names_to_check <- names(upset_matrix2_table[which(upset_matrix2_table < 50 & upset_matrix2_table > 5)])
for ( x in names_to_check ) {
  upset_matrix_to_check <- upset_matrix2[upset_matrix2[,1] %in% x,]
  label <- colnames(upset_matrix_to_check)[upset_matrix_to_check[1,] == "1"]
  label <- gsub("pt2322-DX1BM-","",label)
  label <- paste(label,collapse = "_")
  write.table(upset_matrix_to_check[,-1],paste(label,"_upset_matrix_to_check.txt",sep=""),quote=F,sep=" ")
  
  myvcf_to_check <- myvcf[rownames(myvcf) %in% rownames(upset_matrix_to_check),]
  writeVcf( myvcf_to_check, paste(label,"_to_check.vcf",sep=""))
}