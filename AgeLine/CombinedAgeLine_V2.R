#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries 
library(ggplot2)
library(nlme)
library(ggeffects)
library(lme4)
library(dplyr)
library(ggpubr)
library(RColorBrewer)


### Read in van boxtel ageline data 
setwd("~/hpc/projects/Ageline/3_Output/ageline")
input_ageline_table <- read.table("input_ageline_table.txt",header=T)
ages <- read.table("ages.txt",header=T)
hg38_autosomal_nonN_genome_size <- 2745186691
input_df <- merge(input_ageline_table, ages, by="DONOR")
input_df$SNV_LOAD_NORM <- input_df$SNV_LOAD/input_df$CALLABLE*hg38_autosomal_nonN_genome_size
input_df$INDEL_LOAD_NORM <- input_df$INDEL_LOAD/input_df$CALLABLE*hg38_autosomal_nonN_genome_size

# Get the linear model
lme_age = lme(SNV_LOAD_NORM ~  AGE, random = ~ - 1 + AGE | DONOR, data = input_df)
healthy_intcpt <- lme_age$coefficients$fixed[1]
healthy_slp <- lme_age$coefficients$fixed[2]
ci_interval <- ggpredict(lme_age, ci.lvl = 0.95)$AGE
age_pval = as.data.frame(summary(lme_age)$tTable["AGE","p-value"])
colnames(age_pval) <- "pval"
age_confint = intervals(lme_age)$fixed["AGE",]
age_confint <- as.data.frame(t(age_confint))
age_confint$Tissue = factor("Blood")

# create data.frame with linear fits of fixed effect
mut_expected = expand.grid(Tissue = "Blood", AGE = c(min(input_df$AGE), max(input_df$AGE)))
mut_expected$fit = predict(lme_age, level=0, newdata=mut_expected)



### Read in campbell data 
campbell_data <- read.table("~/hpc/projects/TallFlow/1_Input/CampbellAgeline/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", sep = "\t", header = T)
colnames(campbell_data)
campbell_data$DataType <- campbell_data$Cell.type2
campbell_HSC <- campbell_data[campbell_data$Cell.type2 == "HSC", ]
campbell_naive <- campbell_data[campbell_data$Cell.type2 == "Naive T", ]
campbell_mem <- campbell_data[campbell_data$Cell.type2 == "Memory T", ]

### Read in TallFlow data 
# pt2322
pt2322_data <- read.table("~/hpc/projects/TallFlow/2_Code/AgeLine/pt2322_ageline_table.txt", sep = " ", header = T)
pt2322_data$age <- rep.int(x = 9.4, times = length(rownames(pt2322_data)))
pt2322_data$nmuts_norm <- (pt2322_data$SNV_LOAD/pt2322_data$CALLABLE)*hg38_autosomal_nonN_genome_size
# pt2229
pt2229_data <- read.table("~/hpc/projects/TallFlow/2_Code/AgeLine/pt2229_ageline_table.txt", sep = " ", header = T)
pt2229_data$age <- rep.int(x = 5.5, times = length(rownames(pt2229_data)))
pt2229_data$nmuts_norm <- (pt2229_data$SNV_LOAD/pt2229_data$CALLABLE)*hg38_autosomal_nonN_genome_size
# pt2283D
pt2283D_data <- read.table("~/hpc/projects/TallFlow/2_Code/AgeLine/pt2283D_ageline_table.txt", sep = " ", header = T)
pt2283D_data$age <- rep.int(x = 10, times = length(rownames(pt2283D_data)))
pt2283D_data$nmuts_norm <- (pt2283D_data$SNV_LOAD/pt2283D_data$CALLABLE)*hg38_autosomal_nonN_genome_size




### Add source data 
input_df$Source <- rep.int("vanBoxtel_HSCT", length(rownames(input_df)))
campbell_HSC$Source <- rep.int("Campbell_HSC", length(rownames(campbell_HSC)))
campbell_naive$Source <- rep.int("Campbell_naive", length(rownames(campbell_naive)))
campbell_mem$Source <- rep.int("Campbell_mem", length(rownames(campbell_mem)))
pt2322_data$Source <- pt2322_data$DONOR#rep.int("pt2322", length(rownames(pt2322_data)))
pt2229_data$Source <- pt2229_data$DONOR#rep.int("pt2229", length(rownames(pt2229_data)))
pt2283D_data$Source <- pt2283D_data$DONOR#rep.int("pt2229", length(rownames(pt2229_data)))


### Subset in order to combine 
input_df_sub <- input_df[ colnames(input_df) %in% c("DONOR", "SAMPLE", "SNV_LOAD_NORM", "AGE", "Source")]
campbell_HSC_sub <- campbell_HSC[ colnames(campbell_HSC) %in% c("Donor", "colony", "Nmut_hsc_as", "Age", "Source") ]
campbell_naive_sub <- campbell_naive[ colnames(campbell_naive) %in% c("Donor", "colony", "Nmut_hsc_as", "Age", "Source") ]
campbell_mem_sub <- campbell_mem[ colnames(campbell_mem) %in% c("Donor", "colony", "Nmut_hsc_as", "Age", "Source") ]
pt2322_data_sub <- pt2322_data[colnames(pt2322_data) %in% c("DONOR", "SAMPLE", "nmuts_norm", "age", "Source")]
pt2229_data_sub <- pt2229_data[colnames(pt2229_data) %in% c("DONOR", "SAMPLE", "nmuts_norm", "age", "Source")]
pt2283D_data_sub <- pt2283D_data[colnames(pt2283D_data) %in% c("DONOR", "SAMPLE", "nmuts_norm", "age", "Source")]


# Rename colnames 
colnames(input_df_sub) <- c("DONOR", "SAMPLE", "AGE", "SNV_LOAD_NORM", "SOURCE")
colnames(campbell_HSC_sub) <- c("SAMPLE", "DONOR", "AGE", "SNV_LOAD_NORM", "SOURCE")
colnames(campbell_naive_sub) <- c("SAMPLE", "DONOR", "AGE", "SNV_LOAD_NORM", "SOURCE")
colnames(campbell_mem_sub) <- c("SAMPLE", "DONOR", "AGE", "SNV_LOAD_NORM", "SOURCE")
colnames(pt2322_data_sub) <- c("DONOR", "SAMPLE", "AGE", "SNV_LOAD_NORM", "SOURCE")
colnames(pt2229_data_sub) <- c("DONOR", "SAMPLE", "AGE", "SNV_LOAD_NORM", "SOURCE")
colnames(pt2283D_data_sub) <- c("DONOR", "SAMPLE", "AGE", "SNV_LOAD_NORM", "SOURCE")


### Combine the data 
combined_df <- rbind(input_df_sub, campbell_HSC_sub, campbell_naive_sub, campbell_mem_sub, pt2322_data_sub, pt2229_data_sub, pt2283D_data_sub)


### Values from Campbell lines 
HSC_slope <- 16
HSC_intercept <- 105
naive_slope <- 22
naive_intercept <- 164
mem_slope <- 25
mem_intercept <- 382


### Calculate predicted muts and ratios 
# HSC
combined_df$HSC_PRED_MUT = HSC_intercept + HSC_slope * combined_df$AGE
combined_df$HSC_RATIO = combined_df$SNV_LOAD_NORM / combined_df$HSC_PRED_MUT
# naive
combined_df$naive_PRED_MUT = naive_intercept + naive_slope * combined_df$AGE
combined_df$naive_RATIO = combined_df$SNV_LOAD_NORM / combined_df$naive_PRED_MUT
# mem
combined_df$mem_PRED_MUT = mem_intercept + mem_slope * combined_df$AGE
combined_df$mem_RATIO = combined_df$SNV_LOAD_NORM / combined_df$mem_PRED_MUT
# Get the correct levels for the plot 
combined_df$SOURCE <- factor(combined_df$SOURCE, levels = c("Campbell_HSC", "Campbell_naive","Campbell_mem", "vanBoxtel_HSCT", "pt2322", "pt2229", "pt2283D"))


### plot all figures 
# HSC
HSC_comp <- list( c("Campbell_HSC", "vanBoxtel_HSCT"), c("Campbell_HSC", "pt2322"), c("Campbell_HSC", "pt2229"), c("Campbell_HSC", "pt2283D"), c("pt2322", "pt2229"), c("pt2322", "pt2283D"), c("pt2229", "pt2283D"))
p1 <- ggplot(data = combined_df[combined_df$SOURCE %in% c("vanBoxtel_HSCT", "Campbell_HSC", "pt2229", "pt2322", "pt2283D"),],
       aes(x = SOURCE, y = HSC_RATIO, color = SOURCE, fill = SOURCE)) + 
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  #geom_point() +
  geom_jitter(position = position_jitter(0.15) ) +
  scale_color_manual(values = c("grey", "black", "red", "orange", "darkgreen")) +
  scale_fill_manual(values = c("grey", "black", "red", "orange", "darkgreen")) + 
  theme_classic() + 
  stat_compare_means(comparisons = HSC_comp, label = "p.signif") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# Naive T-cell
naive_comp <- list( c("Campbell_naive", "vanBoxtel_HSCT"), c("Campbell_naive", "pt2322"), c("Campbell_naive", "pt2229"), c("Campbell_naive", "pt2283D"), c("pt2322", "pt2229"), c("pt2229", "pt2283D") )
p2 <- ggplot(data = combined_df[combined_df$SOURCE %in% c("vanBoxtel_HSCT", "Campbell_naive",  "pt2229", "pt2322", "pt2283D"),],
       aes(x = SOURCE, y = naive_RATIO, color = SOURCE, fill = SOURCE)) + 
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  #geom_point() +
  geom_jitter(position = position_jitter(0.15) ) +
  scale_color_manual(values = c("grey", "black", "red", "orange", "darkgreen")) +
  scale_fill_manual(values = c("grey", "black", "red", "orange", "darkgreen")) + 
  theme_classic() + 
  stat_compare_means(comparisons = naive_comp, label = "p.signif") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# Memory T-cell 
mem_comp <- list( c("Campbell_mem", "vanBoxtel_HSCT"), c("Campbell_mem", "pt2322"), c("Campbell_mem", "pt2229"),  c("Campbell_mem", "pt2283D"), c("pt2322", "pt2229") , c("pt2229", "pt2283D") )
p3 <- ggplot(data = combined_df[combined_df$SOURCE %in% c("vanBoxtel_HSCT", "Campbell_mem",  "pt2229","pt2322", "pt2283D"),],
       aes(x = SOURCE, y = mem_RATIO, color = SOURCE, fill = SOURCE)) + 
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  #geom_point() +
  geom_jitter(position = position_jitter(0.15) ) +
  scale_color_manual(values = c("grey", "black", "red", "orange", "darkgreen")) +
  scale_fill_manual(values = c("grey", "black", "red", "orange", "darkgreen")) + 
  theme_classic() + 
  stat_compare_means(comparisons = mem_comp, label = "p.signif") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


#
p1
p2
p3


### Save figures 
ggsave("/Users/verapoort/hpc/projects/TallFlow/3_Output/Ageline/Boxplot_HSC.pdf", p1)
ggsave("/Users/verapoort/hpc/projects/TallFlow/3_Output/Ageline/Boxplot_naive.pdf", p2)
ggsave("/Users/verapoort/hpc/projects/TallFlow/3_Output/Ageline/Boxplot_mem.pdf", p3)


write.table(combined_df, "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/AgeLine/Ratios/Boxplot_matrix.txt", sep = "\t", col.names = NA)


mean(pt2229_data$nmuts_norm)
sd(pt2229_data$nmuts_norm)
mean(pt2322_data$nmuts_norm)
sd(pt2322_data$nmuts_norm)
mean(pt2283D_data$nmuts_norm)
sd(pt2283D_data$nmuts_norm)













### Check the mutational load from our PMC patients
pmc_pts <- read.csv("/Users/verapoort/hpc/projects/TallFlow/1_Input/CampbellAgeline/tall_age.csv")
pmc_pts_2 <- read.csv("/Users/verapoort/hpc/projects/TallFlow/1_Input/CampbellAgeline/tall_age_V2.csv")
colnames(pmc_pts)



pdf("Added_PMC_PTS_2.pdf")
ggplot() +
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='grey', color='grey', shape=21, size=1) + 
  geom_abline(intercept = 105, slope = 16, color='grey') + 
  geom_point(campbell_naive, mapping = aes(Age, Nmut_hsc_as), fill='chocolate4', color='chocolate4', shape=21, size=1) + 
  geom_abline(intercept = 164, slope = 22, color='chocolate4') + 
  geom_point(campbell_mem, mapping = aes(Age, Nmut_hsc_as), fill='darkorchid4', color='darkorchid4', shape=21, size=1) + 
  geom_abline(intercept = 382, slope = 25, color='darkorchid4') + 
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=1) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  geom_point(data = pt2322_data, mapping = aes( y = nmuts_norm, x = age), fill='orange', color='orange', shape=17, size=1) + 
  geom_point(data = pt2229_data, mapping = aes( y = nmuts_norm, x = age), fill='red', color='red', shape=17, size=1) +
  geom_point(data = pt2283D_data, mapping = aes( y = nmuts_norm, x = age), fill='darkgreen', color='darkgreen', shape=17, size=1) + 
  geom_point(pmc_pts_2, mapping = aes(age_at_diagnosis, mutations_after_filters), fill='gold1', color='gold1', shape=17, size=1) +
  #geom_point(pmc_pts, mapping = aes(age_at_diagnosis, mutations_after_filters), fill='darkgrey', color='darkgrey', shape=21, size=1) +
  theme_classic() + 
  ylab("Normalised mutations")
dev.off()





### Include Target 
SuppTable <- read.table("/Users/verapoort/hpc/projects/TallFlow/1_Input/CampbellAgeline/SuppTable1.txt", sep = "\t", header = T)
TargetBurden <- read.table("/Users/verapoort/hpc/projects/TallFlow/1_Input/CampbellAgeline/MutBurdens.txt", sep = "\t", header = T, dec = ".", quote = '"')



SuppTable.sub <- SuppTable[(SuppTable$Lineage == "T") & (SuppTable$WGS == "Yes" ), ]
TargetBurden.sub <- TargetBurden[TargetBurden$PID %in% SuppTable.sub$PatientID, ]
TargetBurden.sub2 <- TargetBurden.sub[TargetBurden.sub$CANCER_TYPE == "TALL", ]
TargetBurden.sub2$AGE <- gsub(pattern = "age=", "", TargetBurden.sub2$AGE)
#float(TargetBurden.sub2$AGE)
TargetBurden.sub2$AGE <- as.integer(TargetBurden.sub2$AGE)
ggplot() +
  #geom_point(pmc_pts, mapping = aes(age_at_diagnosis, mutations_after_filters), fill='darkgrey', color='darkgrey', shape=21, size=1) +
  geom_point(TargetBurden.sub2, mapping = aes(AGE, SBS), fill='SALMON', color='SALMON', shape=21, size=1) +
  theme_bw() + 
  ylab("Normalised mutations")

getwd()

pdf("Added_PMC_PTS_2_withTarget.pdf")
ggplot() +
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='grey', color='grey', shape=21, size=3) + 
  geom_abline(intercept = 105, slope = 16, color='grey') + 
  geom_point(campbell_naive, mapping = aes(Age, Nmut_hsc_as), fill='grey45', color='grey45', shape=17, size=3) + 
  geom_abline(intercept = 164, slope = 22, color='grey45') + 
  geom_point(campbell_mem, mapping = aes(Age, Nmut_hsc_as), fill='grey32', color='grey32', shape=23, size=3) + 
  geom_abline(intercept = 382, slope = 25, color='grey32') + 
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=3) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  geom_point(data = pt2322_data, mapping = aes( y = nmuts_norm, x = age), fill='orange', color='orange', shape=21, size=3) + 
  geom_point(data = pt2229_data, mapping = aes( y = nmuts_norm, x = age), fill='red', color='red', shape=21, size=3) +
  geom_point(data = pt2283D_data, mapping = aes( y = nmuts_norm, x = age), fill='darkgreen', color='darkgreen', shape=21, size=3) + 
  geom_point(pmc_pts_2, mapping = aes(age_at_diagnosis, mutations_after_filters), fill='purple', color='purple', shape=22, size=3) +
  #geom_point(pmc_pts, mapping = aes(age_at_diagnosis, mutations_after_filters), fill='darkgrey', color='darkgrey', shape=21, size=1) +
  geom_point(TargetBurden.sub2, mapping = aes(AGE, SBS), fill='cornflowerblue', color='cornflowerblue', shape=22, size=3) +
  ylab("Normalised mutations") +
  theme_classic() +
  ylim(0, 4000)+ 
  theme(axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text = element_text(size=30),
        legend.text = element_text(size=30))
dev.off()

getwd()
# add the bulk sample to the Agelines, mutations are also normalized to callable loci in script for mutational pattern analysis
# see supplement figure 3C
Nmut_bulk <- c(667,317,97)
pt <- c("pt2229", "pt2322", "pt2283D")
Age <- c(5.5, 9.4, 10)
BulkMuts <- data.frame(Nmut_bulk, pt, Age)


pdf("/Users/verapoort/hpc/projects/TallFlow/3_Output/Ageline/Naive_Ageline_TALL.pdf")
ggplot() +
  geom_point(BulkMuts, mapping = aes(Age, Nmut_bulk), color = "black", fill = "blue", shape = 25, size = 3)+
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='grey', color='grey', shape=21, size=3) + 
  geom_abline(intercept = 105, slope = 16, color='grey') + 
  geom_point(campbell_naive, mapping = aes(Age, Nmut_hsc_as), fill='grey45', color='grey45', shape=17, size=3) + 
  geom_abline(intercept = 164, slope = 22, color='grey45') + 
  #geom_point(campbell_mem, mapping = aes(Age, Nmut_hsc_as), fill='grey32', color='grey32', shape=23, size=3) + 
  #geom_abline(intercept = 382, slope = 25, color='grey32') + 
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=3) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  geom_point(data = pt2322_data, mapping = aes( y = nmuts_norm, x = age), fill='orange', color='orange', shape=21, size=3) + 
  geom_point(data = pt2229_data, mapping = aes( y = nmuts_norm, x = age), fill='red', color='red', shape=21, size=3) + 
  geom_point(data = pt2283D_data, mapping = aes( y = nmuts_norm, x = age), fill='darkgreen', color='darkgreen', shape=21, size=3) +
  theme_classic() + 
  ylab("Normalised mutations") + 
  ylim(0, 2500)+ 
  theme(axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text = element_text(size=30),
        legend.text = element_text(size=30))
dev.off()




geom_point(color='black', shape=21, size=4, aes(fill=factor(DataType))) + 
  scale_fill_manual(values=c('black','Orange', 'red'))+
  #geom_abline(intercept = 164, slope = 22, color = "black") + #Tnaive
  geom_abline(intercept = 382, slope = 25, color = "black") + # Tmem
  #geom_abline(intercept = 105, slope = 16, color = "Black") + # HSC
  theme_classic() +
  ylim(0,4000)+
  ylab("Base substitutions per genome")+
  theme(axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text = element_text(size=30),
        legend.text = element_text(size=30))












