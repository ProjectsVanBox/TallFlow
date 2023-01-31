library(ggplot2)
library(nlme)
library(ggeffects)
library(lme4)

setwd("~/hpc/pmc_vanboxtel/projects/Ageline/3_Output/ageline")

input_ageline_table <- read.table("input_ageline_table.txt",header=T)
ages <- read.table("ages.txt",header=T)
hg38_autosomal_nonN_genome_size <- 2745186691

input_df <- merge(input_ageline_table, ages, by="DONOR")
input_df$SNV_LOAD_NORM <- input_df$SNV_LOAD/input_df$CALLABLE*hg38_autosomal_nonN_genome_size
input_df$INDEL_LOAD_NORM <- input_df$INDEL_LOAD/input_df$CALLABLE*hg38_autosomal_nonN_genome_size

lme_age = lme(SNV_LOAD_NORM ~  AGE, random = ~ - 1 + AGE | DONOR, data = input_df)

lme_age2 = lmer(SNV_LOAD_NORM ~ AGE +( -1 + AGE | DONOR), data = input_df)
pi = ggpredict(lme_age2,"AGE",type="re")

healthy_intcpt <- lme_age$coefficients$fixed[1]
healthy_slp <- lme_age$coefficients$fixed[2]

ci_interval <- ggpredict(lme_age, ci.lvl = 0.95)$AGE
head(ci_interval)

age_pval = as.data.frame(summary(lme_age)$tTable["AGE","p-value"])
colnames(age_pval) <- "pval"
age_confint = intervals(lme_age)$fixed["AGE",]
age_confint <- as.data.frame(t(age_confint))
age_confint$Tissue = factor("Blood")

# create data.frame with linear fits of fixed effect
mut_expected = expand.grid(Tissue = "Blood", AGE = c(min(input_df$AGE), max(input_df$AGE)))
mut_expected$fit = predict(lme_age, level=0, newdata=mut_expected)

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

pdf("ageline.pdf")
print(age_plot)
dev.off()


