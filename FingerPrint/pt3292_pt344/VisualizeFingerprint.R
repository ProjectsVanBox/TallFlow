#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)
#setwd("~/hpc/pmc_vanboxtel/projects/_external_projects/GenomeTox/3_Output/FingerPrint/")
setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FingerPrint/pt3292_pt344/")

library(ggplot2)
library(gplots)
library(reshape2)
library(VariantAnnotation)
library(UpSetR)
library(tidyverse)
library(ggpubr)
library(BSgenome.Hsapiens.NCBI.GRCh38)
ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38'
library(MutationalPatterns)


# Read in data
fingerprint_files <- list.files(path=".",pattern=".vcf$", full.names = T, recursive = T)
vcfs = lapply(fingerprint_files, function(file) readVcf(file))


# Function to convert the genotyping in solid values
convert_GT <- function( vcf ) {
  x = geno(vcf)$GT
  x[,1] <- c(as.numeric(substr(x,1,1))+as.numeric(substr(x,3,3)))
  return(x)
}


# Visualize the fingerprint data 
df <- as.data.frame(lapply(vcfs, function(x) convert_GT(x)))
df <- as.matrix(sapply(df, as.numeric))
pdf("GenomeToxFingerPrint.pdf")
pheatmap::pheatmap(na.omit(df), fontsize_col = 6)
dev.off()









