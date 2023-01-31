#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(ggplot2)
library(MutationalPatterns)
library(BSgenome)

### Read in genomes
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)


### Read in data 
campbell_data <- read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", sep = "\t", header = T)
colnames(campbell_data)
campbell_data$DataType <- campbell_data$Cell.type2



### General plot 
ggplot(campbell_data, aes(Age, Nmut)) + 
  geom_point() + 
  theme_bw()


### Select T-naive cells 
Tnaive <- campbell_data[campbell_data$Cell.type2 == "Naive T", ]

Tnaive_ped <- Tnaive[Tnaive$Age < 20, ]


ggplot(Tnaive, aes(Age, Nmut)) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  theme_bw()
ggplot(Tnaive_ped, aes(Age, Nmut)) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  theme_bw()
#table(campbell_data$CellType)
#table(campbell_data$Cell.type2)




### Read in our data 
pt2322_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt2322_sample_names <- c("pt2322_ALLBULK", "pt2322_TALL1", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12", "pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4",
                         "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", "pt2322_TALL9") 
pt2322_grl <- read_vcfs_as_granges(pt2322_vcf_files, pt2322_sample_names, ref_genome)
pt2322_type_occurrences <- mut_type_occurrences(pt2322_grl, ref_genome)
pt2322_type_occurrences$Nmut <- rowSums(pt2322_type_occurrences[,1:6])
pt2322_type_occurrences$Age <- rep.int(9, nrow(pt2322_type_occurrences))
pt2322_type_occurrences$DataType <- rep.int("pt2322", nrow(pt2322_type_occurrences))

pt2229_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2229/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt2229_sample_names <- c("pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL10", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", 
                         "pt2229_TALL7", "pt2229_TALL8", "pt2229_TALL9")
pt2229_grl <- read_vcfs_as_granges(pt2229_vcf_files, pt2229_sample_names, ref_genome)
pt2229_type_occurrences <- mut_type_occurrences(pt2229_grl, ref_genome)
pt2229_type_occurrences$Nmut <- rowSums(pt2229_type_occurrences[,1:6])
pt2229_type_occurrences$Age <- rep.int(5, nrow(pt2229_type_occurrences))
pt2229_type_occurrences$DataType <- rep.int("pt2229", nrow(pt2229_type_occurrences))



pt_tnaive <- rbind(pt2322_type_occurrences[,c("Age", "Nmut", "DataType")], 
                       pt2229_type_occurrences[,c("Age", "Nmut", "DataType")], 
                       Tnaive[,c("Age", "Nmut", "DataType")])
HSC <- campbell_data[campbell_data$Cell.type2 == "HSC", ]
Tmem <- campbell_data[campbell_data$Cell.type2 == "Memory T", ]
pt_HSC <- rbind(pt2322_type_occurrences[,c("Age", "Nmut", "DataType")], 
                pt2229_type_occurrences[,c("Age", "Nmut", "DataType")], 
                HSC[,c("Age", "Nmut", "DataType")])
pt_Tmem <- rbind(pt2322_type_occurrences[,c("Age", "Nmut", "DataType")], 
                 pt2229_type_occurrences[,c("Age", "Nmut", "DataType")], 
                 Tmem[,c("Age", "Nmut", "DataType")])





norm_tnaive_val <- 9*22+164
pt2322_norm <- pt2322_type_occurrences[,c("Age", "Nmut", "DataType")]
pt2322_norm$Nmut <- pt2322_norm$Nmut/norm_tnaive_val
norm_tnaive_val <- 5*22+164
pt2229_norm <- pt2229_type_occurrences[,c("Age", "Nmut", "DataType")]
pt2229_norm$Nmut <- pt2229_norm$Nmut/norm_tnaive_val
pt_tnaive_norm <- rbind(pt2322_norm, pt2229_norm)

pt_tnaive_norm$DataType <- factor(pt_tnaive_norm$DataType, levels = c("pt2322", "pt2229"))


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/AgeLineBoxplot_test.pdf")
ggplot(data = pt_tnaive_norm) + 
  geom_boxplot(aes(x = DataType, y = Nmut)) +
  geom_point(aes(x = DataType, y = Nmut))
dev.off()

ggplot(pt_tnaive_norm, aes(Age, Nmut, color = DataType)) + 
  geom_point() + 
  geom_abline(intercept = 164, slope = 22) +
  theme_bw()



#pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Ageline_V1.pdf")
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Ageline_Tnaive.pdf")
ggplot(pt_tnaive, aes(Age, Nmut, color = DataType)) + 
  geom_point() + 
  geom_abline(intercept = 164, slope = 22) +
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Ageline_V2.pdf")
ggplot(pt_tnaive, aes(Age, Nmut, color = DataType)) + 
  geom_point() + 
  geom_abline(intercept = 164, slope = 22, color = "NavyBlue") + 
  geom_abline(intercept = 382, slope = 25, color = "Green") + 
  theme_bw()
dev.off()






pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Ageline_HSC.pdf")
ggplot(pt_HSC, aes(Age, Nmut, color = DataType)) + 
  geom_point() + 
  geom_abline(intercept = 105, slope = 16) + 
  theme_bw()
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Ageline_Tmem.pdf")
ggplot(pt_Tmem, aes(Age, Nmut, color = DataType)) + 
  geom_point() + 
  geom_abline(intercept = 382, slope = 25) + 
  theme_bw()
dev.off()

HSC_T <- campbell_data[campbell_data$Cell.type2 %in% c("HSC", "Naive T", "Memory T"),]
pt_all <- rbind(pt2322_type_occurrences[,c("Age", "Nmut", "DataType")], 
                   pt2229_type_occurrences[,c("Age", "Nmut", "DataType")], 
                HSC_T[,c("Age", "Nmut", "DataType")])
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Ageline_HSC_T.pdf")
ggplot(pt_all, aes(Age, Nmut, color = DataType)) + 
  geom_point() + 
  geom_abline(intercept = 164, slope = 22, color = "NavyBlue") + 
  geom_abline(intercept = 382, slope = 25, color = "Orange") +
  geom_abline(intercept = 105, slope = 16, color = "Green") + 
  theme_bw()
dev.off()
# + geom_abline(intercept = 37, slope = -5)

COLORS6 <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
)
COLORS7 <- c(
  "#2EBAED", "#000000", 
  "#E98C7B", "#DE1C14","#D4D2D2", "#ADCC54",
  "#F0D0CE"
)
plot_spectrum_stacked(pt2322_type_occurrences , max_value_to_show = 10)