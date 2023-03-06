
library(VariantAnnotation)
library(ggplot2)
library(reshape)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TruePosRate/pt2322/")

myBulkVcf_fname <- "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/intermediate/short_variants/somatic_vcfs/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-ALLBULK.SMuRF.filtered.vcf.gz"

myBulkVcf <- readVcf(myBulkVcf_fname)

clonal_table <- table(info(myBulkVcf)$CLONAL_SAMPLE_NAMES)
colnames(clonal_table) <- gsub(".+BM-(.+)","\\1",colnames(clonal_table))

clonal_table_2 <- as.data.frame(apply(clonal_table,2,paste))
clonal_table_df <- clonal_table_2
w <- which(clonal_table_2=="1",arr.ind=TRUE)
clonal_table_2[w] <- names(clonal_table_2)[w[,"col"]]
clonal_table_2[clonal_table_2 == 0] <- NA

TALL1_missed <- clonal_table[clonal_table[,"TALL1"] == 0 & (clonal_table[,"TALL3"] == 1 | rowSums(clonal_table) >= 11),]
TALL2_missed <- clonal_table[clonal_table[,"TALL2"] == 0 & ( (( clonal_table[,"TALL1"] == 1 | clonal_table[,"TALL3"] == 1  ) & (clonal_table[,"TALL7"] == 1 | clonal_table[,"TALL12"] == 1) ) | rowSums(clonal_table) >= 11),] 
TALL3_missed <- clonal_table[clonal_table[,"TALL3"] == 0 & (clonal_table[,"TALL1"] == 1 | rowSums(clonal_table) >= 11),]

TALL4_missed <- clonal_table[clonal_table[,"TALL4"] == 0 & (clonal_table[,"TALL5"] == 1 | rowSums(clonal_table) >= 11),]
TALL5_missed <- clonal_table[clonal_table[,"TALL5"] == 0 & (clonal_table[,"TALL4"] == 1 | rowSums(clonal_table) >= 11),]
TALL6_missed <- clonal_table[clonal_table[,"TALL6"] == 0 & ( ((clonal_table[,"TALL4"] == 1 | clonal_table[,"TALL5"] == 1  ) & clonal_table[,"TALL8"] == 1) | rowSums(clonal_table) >= 11),] 

TALL7_missed <- clonal_table[clonal_table[,"TALL7"] == 0 & (clonal_table[,"TALL12"] == 1 | rowSums(clonal_table) >= 11),]
TALL8_missed <- clonal_table[clonal_table[,"TALL8"] == 0 & rowSums(clonal_table) >= 11,]
TALL9_missed <- clonal_table[clonal_table[,"TALL9"] == 0 & (clonal_table[,"TALL11"] == 1 | rowSums(clonal_table) >= 11),]

TALL10_missed <- clonal_table[clonal_table[,"TALL10"] == 0 & ( ((clonal_table[,"TALL9"] == 1 | clonal_table[,"TALL11"] == 1  ) & (clonal_table[,"TALL1"] == 1 | clonal_table[,"TALL3"] == 1 | clonal_table[,"TALL2"] | clonal_table[,"TALL7"] | clonal_table[,"TALL12"]  )) | rowSums(clonal_table) >= 11),] 
TALL11_missed <- clonal_table[clonal_table[,"TALL11"] == 0 & (clonal_table[,"TALL9"] == 1 | rowSums(clonal_table) >= 11),]
TALL12_missed <- clonal_table[clonal_table[,"TALL12"] == 0 & (clonal_table[,"TALL7"] == 1 | rowSums(clonal_table) >= 11),]

TALL1_missed <- clonal_table[clonal_table[,"TALL1"] == 0 & rowSums(clonal_table) >= 11,]
TALL2_missed <- clonal_table[clonal_table[,"TALL2"] == 0 & rowSums(clonal_table) >= 11,]
TALL3_missed <- clonal_table[clonal_table[,"TALL3"] == 0 & rowSums(clonal_table) >= 11,]

TALL4_missed <- clonal_table[clonal_table[,"TALL4"] == 0 & rowSums(clonal_table) >= 11,]
TALL5_missed <- clonal_table[clonal_table[,"TALL5"] == 0 & rowSums(clonal_table) >= 11,]
TALL6_missed <- clonal_table[clonal_table[,"TALL6"] == 0 & rowSums(clonal_table) >= 11,]

TALL7_missed <- clonal_table[clonal_table[,"TALL7"] == 0 & rowSums(clonal_table) >= 11,]
TALL8_missed <- clonal_table[clonal_table[,"TALL8"] == 0 & rowSums(clonal_table) >= 11,]
TALL9_missed <- clonal_table[clonal_table[,"TALL9"] == 0 & rowSums(clonal_table) >= 11,]

TALL10_missed <- clonal_table[clonal_table[,"TALL10"] == 0 & rowSums(clonal_table) >= 11,]
TALL11_missed <- clonal_table[clonal_table[,"TALL11"] == 0 & rowSums(clonal_table) >= 11,]
TALL12_missed <- clonal_table[clonal_table[,"TALL12"] == 0 & rowSums(clonal_table) >= 11,]

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
  "Blast" = rep(colnames(clonal_table)[-1],each=1),
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

# mydf <- data.frame(
#   "cell" = rep(colnames(clonal_table)[-1],each=3),
#   "label" = rep(c("total","missed","recall"),12),
#   "x" = c(
#     nrow(TALL1_found)+nrow(TALL1_missed),nrow(TALL1_missed),1-nrow(TALL1_missed)/(nrow(TALL1_found)+nrow(TALL1_missed)),
#     nrow(TALL10_found)+nrow(TALL10_missed),nrow(TALL10_missed),1-nrow(TALL10_missed)/(nrow(TALL10_found)+nrow(TALL10_missed)),
#     nrow(TALL11_found)+nrow(TALL11_missed),nrow(TALL11_missed),1-nrow(TALL11_missed)/(nrow(TALL11_found)+nrow(TALL11_missed)),
#     nrow(TALL12_found)+nrow(TALL12_missed),nrow(TALL12_missed),1-nrow(TALL12_missed)/(nrow(TALL12_found)+nrow(TALL12_missed)),
#     nrow(TALL2_found)+nrow(TALL2_missed),nrow(TALL2_missed),1-nrow(TALL2_missed)/(nrow(TALL2_found)+nrow(TALL2_missed)),
#     nrow(TALL3_found)+nrow(TALL3_missed),nrow(TALL3_missed),1-nrow(TALL3_missed)/(nrow(TALL3_found)+nrow(TALL3_missed)),
#     nrow(TALL4_found)+nrow(TALL4_missed),nrow(TALL4_missed),1-nrow(TALL4_missed)/(nrow(TALL4_found)+nrow(TALL4_missed)),
#     nrow(TALL5_found)+nrow(TALL5_missed),nrow(TALL5_missed),1-nrow(TALL5_missed)/(nrow(TALL5_found)+nrow(TALL5_missed)),
#     nrow(TALL6_found)+nrow(TALL6_missed),nrow(TALL6_missed),1-nrow(TALL6_missed)/(nrow(TALL6_found)+nrow(TALL6_missed)),
#     nrow(TALL7_found)+nrow(TALL7_missed),nrow(TALL7_missed),1-nrow(TALL7_missed)/(nrow(TALL7_found)+nrow(TALL7_missed)),
#     nrow(TALL8_found)+nrow(TALL8_missed),nrow(TALL8_missed),1-nrow(TALL8_missed)/(nrow(TALL8_found)+nrow(TALL8_missed)),
#     nrow(TALL9_found)+nrow(TALL9_missed),nrow(TALL9_missed),1-nrow(TALL9_missed)/(nrow(TALL9_found)+nrow(TALL9_missed))
#   )
#   )
# 
# ggplot(mydf[mydf$label =='recall',],aes(x=cell,y=x)) +
#   geom_bar(stat="identity")

