#!/usr/bin/env Rscript


### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries 
library(ggplot2)
library(VariantAnnotation)
library(pheatmap)
library(stringr)


### Read in data
pt2322_vcf <- readVcf("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2322/DriverSelection/pt2322_DriverVafs.vcf")
#pt2229_vcf <- readVcf("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2229/DriverSelection/pt2229_DriverVafs.vcf")
pt2229_vcf <- readVcf("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2229inclDP_RH/DriverSelection/pt2229_DriverVafs.vcf")


### Parse variables for pt2322
pt2322_df <- data.frame(geno(pt2322_vcf)$VAF)
pt2322_df <- pt2322_df[c("pt2322.DX1BM.MONOBULK", "pt2322.DX1BM.ALLBULK", "pt2322.DX1BM.ALLBULK.DN", "pt2322.DX1BM.ALLBULK.SPCD4", "pt2322.DX1BM.ALLBULK.DP", "pt2322.DX1BM.ALLBULK.iSPCD4")]
pt2322_list <- lapply(info(pt2322_vcf)$ANN,function(x){strsplit(x,"\\|")})
pt2322_genes <- c(sapply(sapply(pt2322_list, "[", 1),"[",4))
pt2322_muts <- c(sapply(sapply(pt2322_list, "[", 1),"[",10))
pt2322_rownames <- paste0(pt2322_genes, "_", pt2322_muts)
pt2322_HmInput <- data.matrix(pt2322_df)
rownames(pt2322_HmInput) <- pt2322_rownames


### Parse variables for pt2229
pt2229_df <- data.frame(geno(pt2229_vcf)$VAF)
pt2229_df <- pt2229_df[c("pt2229.DX1BM.MSCBULK", "pt2229.DX1BM.ALLBULK", "pt2229.DX1BM.ALLBULK.DNCD1aNeg", "pt2229.DX1BM.ALLBULK.DNCD1aPos", "pt2229.DX1BM.ALLBULK.iSPCD4", "pt2229.DX1BM.ALLBULK.DP")]
pt2229_list <- lapply(info(pt2229_vcf)$ANN,function(x){strsplit(x,"\\|")})
pt2229_genes <- c(sapply(sapply(pt2229_list, "[", 1),"[",4))
pt2229_muts <- c(sapply(sapply(pt2229_list, "[", 1),"[",10))
pt2229_rownames <- paste0(pt2229_genes, "_", pt2229_muts)
pt2229_HmInput <- data.matrix(pt2229_df)
rownames(pt2229_HmInput) <- pt2229_rownames

            
### Create heatmap 
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/DriverVaf/pt2322_Hm.pdf")
pheatmap(pt2322_HmInput, display_numbers = T, cluster_cols = F)
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/DriverVaf/pt2229_Hm.pdf")
pheatmap(pt2229_HmInput, display_numbers = T, cluster_cols = F)
dev.off()














