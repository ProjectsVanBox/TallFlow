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
library(MutationalPatterns)
library(BSgenome)
library(ggpubr)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

setwd("~/hpc/pmc_vanboxtel/projects/Ageline/3_Output/ageline")

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


# Plot
age_plot = ggplot() +
  geom_ribbon(data=pi, aes(x=x, ymin=conf.low, ymax=conf.high),fill = "#EEEEEE") +
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high),fill = "#CCCCCC") +
  geom_line(data=mut_expected, aes(y = fit, x = AGE), size=1.5, col="black") +
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black", linetype = "dotted") +
  geom_point(data = input_df, aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=5) +
  ylab("Base substitutions per genome") +
  geom_text(data = age_pval, aes(x = 35, y = 1400,label=paste("P =",format(pval,scientific = F,digits = 2))), size=3.5, colour="black") +
  scale_y_continuous(expand = c(0,0), limits = c(-65,1500),breaks = c(0,500,1000,1500)) +
  scale_x_continuous(expand = c(0,0),limits = c(-2,70), breaks = c(seq(0, 70, 10))) +
  geom_segment(aes(x = 0, xend = 70, y = -Inf, yend = -Inf), col = "black", size = 0.5) +
  geom_segment(aes(y = 0, yend = 1500, x = -Inf, xend = -Inf), col = "black") +
  theme_classic() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=10), strip.text = element_text(size=15), legend.text = element_text(size=13)) +
  theme(legend.position="right",axis.line = element_blank()) +
  xlab("Age (years)")
age_plot



### Read in campbell data 
campbell_data <- read.table("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", sep = "\t", header = T)
colnames(campbell_data)
campbell_data$DataType <- campbell_data$Cell.type2
campbell_HSC <- campbell_data[campbell_data$Cell.type2 == "HSC", ]


campbell_lme_age = lme(Nmut_hsc_as ~  Age, random = ~ - 1 + Age | Donor, data = campbell_HSC)
campbell_intcpt <- campbell_lme_age$coefficients$fixed[1]
campbell_slp <- campbell_lme_age$coefficients$fixed[2]
campbell_ci_interval <- ggpredict(campbell_lme_age, ci.lvl = 0.95)$Age

#campbell_data$Cell.type2 = factor(campbell_data$Cell.type2, levels=c("Naive B", "HSC", "Memory B", "Naive T",  "Memory T"))
m3Allsub = subset(campbell_data, Cell.type2 != "Treg")
full.lymph.Age_cell.interaction.lme <- lme(fixed = Nmut_hsc_as ~ Age * Cell.type2, 
                                           random = ~ 1 | Donor / Cell.type2, 
                                           weights = varIdent(form= ~ 1 | Cell.type2),
                                           data = m3Allsub, method="ML")



ggplot(campbell_HSC, aes(Age, Nmut_hsc_as, color = DataType)) + 
  geom_point() + 
  geom_abline(intercept = 105, slope = 16) + 
  theme_bw()



pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/AgeLine/Age_points_combined.pdf")
ggplot() +
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='red', color='red', shape=21, size=1) + 
  geom_abline(intercept = 105, slope = 16, color='red') + 
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=1) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  theme_bw()
dev.off()




ggplot() +
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='red', color='red', shape=21, size=1) + 
  geom_abline(intercept = 105, slope = 16, color='red', linetype = "dotted") + 
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=1) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  theme_bw()



### Read in our data 
# pt2322
pt2322_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/pt2322_ageline_table.txt", sep = " ", header = T)
pt2322_data$age <- rep.int(x = 9.4, times = length(rownames(pt2322_data)))
pt2322_data$nmuts_norm <- (pt2322_data$SNV_LOAD/pt2322_data$CALLABLE)*hg38_autosomal_nonN_genome_size
pt2322_data$predicted_muts = healthy_intcpt + healthy_slp * pt2322_data$age
pt2322_data$ratio = pt2322_data$nmuts_norm / pt2322_data$predicted_muts
# pt2229
pt2229_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/pt2229_ageline_table.txt", sep = " ", header = T)
pt2229_data$age <- rep.int(x = 5.5, times = length(rownames(pt2229_data)))
pt2229_data$nmuts_norm <- (pt2229_data$SNV_LOAD/pt2229_data$CALLABLE)*hg38_autosomal_nonN_genome_size
pt2229_data$predicted_muts = healthy_intcpt + healthy_slp * pt2229_data$age
pt2229_data$ratio = pt2229_data$nmuts_norm / pt2229_data$predicted_muts


### Add the predicted and ratios to the healthy controls
input_df$PREDICTED_MUTS <- healthy_intcpt + healthy_slp * input_df$AGE
input_df$RATIO = input_df$SNV_LOAD_NORM / input_df$PREDICTED_MUTS
campbell_HSC$PREDICTED_MUTS <- healthy_intcpt + healthy_slp * campbell_HSC$Age
campbell_HSC$RATIO = campbell_HSC$Nmut_hsc_as / campbell_HSC$PREDICTED_MUTS
input_df$Source <- rep.int("Box", length(rownames(input_df)))
campbell_HSC$Source <- rep.int("campbell", length(rownames(campbell_HSC)))


pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/AgeLine/Age_points_HSC.pdf")
ggplot() +
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='red', color='red', shape=21, size=1) + 
  geom_abline(intercept = 105, slope = 16, color='red', linetype = "dotted") + 
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=1) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  geom_point(data = pt2322_data, mapping = aes( y = nmuts_norm, x = age), fill='blue', color='blue', shape=21, size=1) + 
  geom_point(data = pt2229_data, mapping = aes( y = nmuts_norm, x = age), fill='darkgreen', color='darkgreen', shape=21, size=1) + 
  theme_bw() + 
  ylab("Normalised mutations")
dev.off()


p1 <- ggplot() +
  geom_boxplot(data = input_df, mapping = aes(x = DONOR, y = RATIO), alpha = 0.4, fill = "#CCCCCC", color = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes(x = DONOR, y = RATIO), alpha = 0.4, fill = "#CCCCCC", color = "#CCCCCC") +
  geom_boxplot(data = campbell_HSC, mapping = aes(x = Donor, y = RATIO), fill='red', color='red', alpha = 0.4) +
  geom_point(data = campbell_HSC, mapping = aes(x = Donor, y = RATIO), fill='red', color='red', alpha = 0.4) +
  geom_boxplot(data = pt2322_data, mapping = aes(x = DONOR, y = ratio), fill='blue', color='blue', alpha = 0.4) +
  geom_point(data = pt2322_data, mapping = aes(x = DONOR, y = ratio), fill='blue', color='blue') +
  geom_boxplot(data = pt2229_data, mapping = aes(x = DONOR, y = ratio), fill='darkgreen', color='darkgreen', alpha = 0.4) +
  geom_point(data = pt2229_data, mapping = aes(x = DONOR, y = ratio), fill='darkgreen', color='darkgreen') +
  theme_bw()
p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/AgeLine/Boxplot_HSC.pdf")
p1
dev.off()


p2 <- ggplot() +
  geom_boxplot(data = input_df, mapping = aes(x = Source, y = RATIO), alpha = 0.4, fill = "#CCCCCC", color = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes(x = Source, y = RATIO), alpha = 0.4, fill = "#CCCCCC", color = "#CCCCCC") +
  geom_boxplot(data = campbell_HSC, mapping = aes(x = Source, y = RATIO), fill='red', color='red', alpha = 0.4) +
  geom_point(data = campbell_HSC, mapping = aes(x = Source, y = RATIO), fill='red', color='red', alpha = 0.4) +
  geom_boxplot(data = pt2322_data, mapping = aes(x = DONOR, y = ratio), fill='blue', color='blue', alpha = 0.4) +
  geom_point(data = pt2322_data, mapping = aes(x = DONOR, y = ratio), fill='blue', color='blue') +
  geom_boxplot(data = pt2229_data, mapping = aes(x = DONOR, y = ratio), fill='darkgreen', color='darkgreen', alpha = 0.4) +
  geom_point(data = pt2229_data, mapping = aes(x = DONOR, y = ratio), fill='darkgreen', color='darkgreen') +
  theme_bw()
p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/AgeLine/Boxplot_HSC_consise.pdf")
p2
dev.off()



p2 + stat_compare_means(method = "kruskal.test")
p2 + stat_compare_means(method = "kruskal.test", label = "p.signif", label.x = 1.5, label.y = 40)





rownames(pt2322_data)
mut_counts =  data.frame(age = MATCHED_AGE_PER_SAMPLE,
                         cl = CALLABLE_AREA,
                         nmuts = N_SNV_PER_SAMPLE %>%
                           mutate(nmuts_norm = (nmuts/cl)*autosomal))
mut_counts$patient = rownames(mut_counts)
mut_counts$predicted_muts = intercept + slope * mut_counts$age
mut_counts$ratio = mut_counts$nmuts_norm / mut_counts$predicted_muts






# , linetype = "dotted"
ggplot() +
  geom_ribbon(data=campbell_ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "red") +
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='red', color='red', shape=21, size=1) + 
  geom_abline(intercept = 105, slope = 16, color='red', linetype = "dotted") + 
  geom_abline(intercept = campbell_intcpt, slope = campbell_slp, color='red') + 
  geom_point(campbell_HSC, mapping = aes(Age, Nmut), fill='green', color='green', shape=21, size=1) + 
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=1) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  theme_bw()



ggplot() +
  #geom_ribbon(data=campbell_ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "red") +
  #geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='red', color='red', shape=21, size=1) + 
  #geom_abline(intercept = 105, slope = 16, color='red', linetype = "dotted") + 
  #geom_abline(intercept = campbell_intcpt, slope = campbell_slp, color='red') + 
  geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=1) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  theme_bw()

ggplot() +
  #geom_ribbon(data=campbell_ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "red") +
  geom_point(campbell_HSC, mapping = aes(Age, Nmut_hsc_as), fill='red', color='red', shape=21, size=1) + 
  geom_abline(intercept = 105, slope = 16, color='red', linetype = "dotted") + 
  geom_abline(intercept = campbell_intcpt, slope = campbell_slp, color='red') + 
  #geom_ribbon(data=ci_interval, aes(x=x, ymin=conf.low, ymax=conf.high), alpha = 0.4, fill = "#CCCCCC") +
  geom_point(data = input_df, mapping = aes( y = SNV_LOAD_NORM, x = AGE), fill='black', color='black', shape=21, size=1) + 
  geom_abline(intercept = healthy_intcpt, slope = healthy_slp, color="black") +
  theme_bw()











### Old 
pt2322_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt2322_sample_names <- c("pt2322_ALLBULK", "pt2322_TALL1", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12", "pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4",
                         "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", "pt2322_TALL9") 
pt2322_grl <- read_vcfs_as_granges(pt2322_vcf_files, pt2322_sample_names, ref_genome)
pt2322_type_occurrences <- mut_type_occurrences(pt2322_grl, ref_genome)
pt2322_type_occurrences$Nmut <- rowSums(pt2322_type_occurrences[,1:6])
pt2322_type_occurrences$Age <- rep.int(9.4, nrow(pt2322_type_occurrences))
pt2322_type_occurrences$DataType <- rep.int("pt2322", nrow(pt2322_type_occurrences))


