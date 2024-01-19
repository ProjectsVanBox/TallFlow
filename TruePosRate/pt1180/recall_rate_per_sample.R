
library(VariantAnnotation)
library(ggplot2)
library(reshape)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TruePosRate/pt1180/")

myBulkVcf_fname <- "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt1180/231017_HMFreg2090_pt1180_Bulk.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_pt1180-DX1PB-ALLBULK.SMuRF.filtered.vcf"

myBulkVcf <- readVcf(myBulkVcf_fname)
myBulkVcfgr <- granges(myBulkVcf)

mySingleSampleSMuRF_fname <- "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt1180/intermediate/short_variants/SMuRF/pt1180/231016_HMFreg2090_pt1180.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted.SMuRF.filtered.vcf"

mySingleSampleSMuRF <- readVcf(mySingleSampleSMuRF_fname)
mySingleSampleSMuRFgr <- granges(mySingleSampleSMuRF)

subset <- subsetByOverlaps( mySingleSampleSMuRFgr, myBulkVcfgr )

mySingleSampleSMuRF_subset <- mySingleSampleSMuRF[rownames(mySingleSampleSMuRF) %in% names(subset),]

clonal_table <- table(info(mySingleSampleSMuRF_subset)$CLONAL_SAMPLE_NAMES)
colnames(clonal_table) <- gsub(".+PB-(.+)","\\1",colnames(clonal_table))

clonal_table_2 <- as.data.frame(apply(clonal_table,2,paste))
clonal_table_df <- clonal_table_2
w <- which(clonal_table_2=="1",arr.ind=TRUE)
clonal_table_2[w] <- names(clonal_table_2)[w[,"col"]]
clonal_table_2[clonal_table_2 == 0] <- NA

TALL1_missed <- clonal_table[clonal_table[,"TALL1"] == 0 & rowSums(clonal_table) >= 10,]
TALL2_missed <- clonal_table[clonal_table[,"TALL2"] == 0 & rowSums(clonal_table) >= 10,]
TALL3_missed <- clonal_table[clonal_table[,"TALL3"] == 0 & rowSums(clonal_table) >= 10,]

TALL4_missed <- clonal_table[clonal_table[,"TALL4"] == 0 & rowSums(clonal_table) >= 10,]
TALL5_missed <- clonal_table[clonal_table[,"TALL5"] == 0 & rowSums(clonal_table) >= 10,]
TALL6_missed <- clonal_table[clonal_table[,"TALL6"] == 0 & rowSums(clonal_table) >= 10,]

TALL7_missed <- clonal_table[clonal_table[,"TALL7"] == 0 & rowSums(clonal_table) >= 10,]
TALL8_missed <- clonal_table[clonal_table[,"TALL8"] == 0 & rowSums(clonal_table) >= 10,]
TALL9_missed <- clonal_table[clonal_table[,"TALL9"] == 0 & rowSums(clonal_table) >= 10,]

TALL10_missed <- clonal_table[clonal_table[,"TALL10"] == 0 & rowSums(clonal_table) >= 10,]
TALL11_missed <- clonal_table[clonal_table[,"TALL11"] == 0 & rowSums(clonal_table) >= 10,]
TALL12_missed <- clonal_table[clonal_table[,"TALL12"] == 0 & rowSums(clonal_table) >= 10,]

TALL1_found <- clonal_table[clonal_table[,"TALL1"] == 1,]
TALL2_found <- clonal_table[clonal_table[,"TALL2"] == 1,]
TALL3_found <- clonal_table[clonal_table[,"TALL3"] == 1,]

TALL4_found <- clonal_table[clonal_table[,"TALL4"] == 1,]
TALL5_found <- clonal_table[clonal_table[,"TALL5"] == 1,]
TALL6_found <- clonal_table[clonal_table[,"TALL6"] == 1,]

TALL7_found <- clonal_table[clonal_table[,"TALL7"] == 1,]
TALL8_found <- clonal_table[clonal_table[,"TALL8"] == 1,]
TALL9_found <- clonal_table[clonal_table[,"TALL9"] == 1,]

TALL10_found <- clonal_table[clonal_table[,"TALL10"] == 1,]
TALL11_found <- clonal_table[clonal_table[,"TALL11"] == 1,]
TALL12_found <- clonal_table[clonal_table[,"TALL12"] == 1,]

mydf <- data.frame(
  "Blast" = rep(colnames(clonal_table),each=1),
  "Recall_Rate" = c(
    1-nrow(TALL1_missed)/(nrow(TALL1_found)+nrow(TALL1_missed)),
    1-nrow(TALL10_missed)/(nrow(TALL10_found)+nrow(TALL10_missed)),
    1-nrow(TALL11_missed)/(nrow(TALL11_found)+nrow(TALL11_missed)),
    1-nrow(TALL12_missed)/(nrow(TALL12_found)+nrow(TALL12_missed)),
    1-nrow(TALL2_missed)/(nrow(TALL2_found)+nrow(TALL2_missed)),
    1-nrow(TALL3_missed)/(nrow(TALL3_found)+nrow(TALL3_missed)),
    1-nrow(TALL4_missed)/(nrow(TALL4_found)+nrow(TALL4_missed)),
    1-nrow(TALL5_missed)/(nrow(TALL5_found)+nrow(TALL5_missed)),
    1-nrow(TALL6_missed)/(nrow(TALL6_found)+nrow(TALL6_missed)),
    1-nrow(TALL7_missed)/(nrow(TALL7_found)+nrow(TALL7_missed)),
    1-nrow(TALL8_missed)/(nrow(TALL8_found)+nrow(TALL8_missed)),
    1-nrow(TALL9_missed)/(nrow(TALL9_found)+nrow(TALL9_missed))
  )
)

pdf("recall_rate.pdf")
ggplot(mydf,aes(x=Blast,y=Recall_Rate))+
  geom_bar(stat='identity')
dev.off()

write.table(mydf,"recall_rate.txt",sep="",quote=F)

