## compare populations found with flow cytometry with CITE seqeuncing 
### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

## libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(reshape2)
library(tidyverse)

## set dir
setwd("/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallFlow/3_Output/")

# color scheme
pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A", "#E28028", "#963C25", "grey")
names(pal) <- c("ALLBULK", "DN-like", "DNCD1aNeg","CD4-like", "DP-like", "MSCBULK", "MONOBulk","DN3-like", "ISP-like", "CD8-like", "gdTcell-like", "Unknown")

## plotting order
plot_order <- c("DN-like", "DN3-like", "ISP-like", "gdTcell-like", "DP-like", "CD4-like", "CD8-like")

## load data 
# population proportions in the flow cytometry
#flow_pop <- read.table('/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallFlow/3_Output/Barcharts/PatientStrict/OnlyT_Populations/PercentagesLabel_Patients_TPop.txt', header = TRUE)
flow_pop <- read.table('Barcharts/PatientStrict/OnlyT_Populations/LabelCountsTable_Patients_TPop_V3.txt', header = TRUE)
#flow_thy <- read.csv('/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallFlow/3_Output/Barcharts/')
#thy_flow_files <- list.files(path='PythonGating/ThyCITE/',pattern='.csv',full.names = T)
  
thy_pop_list <- lapply(thy_flow_files, read.csv)
names(thy_pop_list) <- c("thy01_undepl", "thy02_undepl", "thy03_undepl", "Thy01", "Thy02", "Thy03")
thy_pop <- bind_rows(thy_pop_list, .id = "id")

thy_pop_count <- table(thy_pop$id, thy_pop$Celltype)
thy_pop_m <- melt(thy_pop_count)
colnames(thy_pop_m) <- c("SampleID", "Celltype", "n")
thy_pop_m$platform <- "flow"

thy_pop_m_mean <- thy_pop_m %>% 
  group_by(SampleID, Celltype, platform) %>% 
  summarise(across(where(is.numeric), mean)) %>% na.omit()
colnames(thy_pop_m_mean) <- c("SampleID", "Celltype", "platform", "n")

## prep flow pop file
flow_pop_t <- t(flow_pop)
#rownames(flow_pop_t) <- colnames(flow_pop)
#colnames(flow_pop_t) <- rownames(flow_pop)
colnames(flow_pop_t) <- flow_pop_t[1,]
flow_pop_t <- flow_pop_t[-1,]
flow_pop_m <- melt(flow_pop_t)
colnames(flow_pop_m) <- c("SampleID", "Celltype", "n")
flow_pop_m$platform <- "flow"
flow_pop_m <- flow_pop_m %>% mutate(SampleID2 = gsub(".*_", "", SampleID))
flow_pop_m$n <- as.numeric(flow_pop_m$n)
head(flow_pop_m)


flow_pop_m_mean <- flow_pop_m %>% 
  group_by(SampleID2, Celltype, platform) %>% 
  summarise(across(where(is.numeric), mean)) %>% na.omit()
colnames(flow_pop_m_mean) <- c("SampleID", "Celltype", "platform", "n")

# population proportions with CITE
# -> this still needs to be exported
CITE_pop <- read.table('~/hpc/projects/TallClonal/3_Output/10x_ThyPtsComb/cellgating.csv', sep = ",", header =TRUE)
CITE_pop$SampleID <- sub("_[^_]+$*", "", CITE_pop$X)

head(CITE_pop)
CITE_pop <- CITE_pop %>% mutate(Celltype = case_when(grepl("SPCD8", x) ~ "CD8-like", 
                                        grepl("DNearly",x) ~ "DN-like", 
                                        grepl("SPCD4", x) ~ "CD4-like", 
                                        grepl("DN3", x) ~ "DN3-like",
                                        grepl("Bcell", x) ~ "B-cell-like", 
                                        grepl("Monocyte", x) ~ "Monocyte-like", 
                                        grepl("NK", x) ~ "NK-cell-like", 
                                        grepl("DP", x) ~ "DP-like", 
                                        grepl("Dendritic", x) ~ "pDC-like", 
                                        grepl("unknown", x) ~ "Unknown", 
                                        grepl("gdTcell", x) ~ "gdTcell-like", 
                                        grepl("iSPCD4", x) ~ "ISP-like"))


#colnames(CITE_pop) <- c("X", "CITE", "index", "V1", "Celltype")
CITE_pop1 <- as.data.frame(table(CITE_pop$SampleID, CITE_pop$Celltype))
colnames(CITE_pop1) <- c("SampleID", "Celltype", "n")
CITE_pop1$platform <- "CITE"

# get naming the same -> of patients but also of the populations
#write.csv(colnames(flow_pop), 'namingFlowCiteComp.csv')
#name_conf <- read.table('hpc/projects/TallClonal/3_Output/10x_ThyPtsComb/namingFlowCiteComp.csv', sep = ",", header = TRUE)

#colnames(name_conf) <- c("SampleID", "Flow", "CITE")
#head(name_conf)

## remove non T populations + Unknown
thy_pop_m_mean <- thy_pop_m_mean[thy_pop_m_mean$Celltype %in% plot_order,]
CITE_pop1 <- CITE_pop1[CITE_pop1$Celltype %in% plot_order,]

## add proportions to whole data set > group.by summarize 

flow_pop_m_mean <- flow_pop_m_mean %>% 
  group_by(SampleID) %>%
  mutate(ncells = sum(n)) %>% 
  mutate(freq = n/ncells)

thy_pop_m_mean <- thy_pop_m_mean %>% 
  group_by(SampleID) %>%
  mutate(ncells = sum(n)) %>% 
  mutate(freq = n/ncells)

CITE_pop1_mean <- CITE_pop1 %>% 
  group_by(SampleID) %>%
  mutate(ncells = sum(n)) %>% 
  mutate(freq = n/ncells)

# filter to keep only Samples that we processed on both platforms
flow_pop_filt <- flow_pop_m_mean %>% filter(SampleID %in% CITE_pop1$SampleID)
thy_pop_filt <- thy_pop_m_mean %>% filter(SampleID %in% CITE_pop1$SampleID)
CITE_pop1_filt <- CITE_pop1_mean %>% filter(SampleID %in% thy_pop_m_mean$SampleID) # filtered against thymus
CITE_pop1_filt2 <- CITE_pop1_mean %>% filter(SampleID %in% flow_pop_m_mean$SampleID) # filtered against pts
#unique(CITE_pop1_filt$SampleID)
#unique(CITE_pop1_filt$ncells)


combined <- bind_rows(thy_pop_filt, CITE_pop1_filt) #%>% full_join(CITE_pop, by="CITE")
combined_pts <- bind_rows(flow_pop_filt, CITE_pop1_filt2) #%>% full_join(CITE_pop, by="CITE")

head(combined_pts)
tail(combined_pts)

# plot paired barplots per patient, what is found in CITE and what in Flow.
# using both relative and absolute values.
# absolute makes no sense as we checked 500.000 cells for flow and ~ 5000 for CITE

# statistics required? 
#https://radiant-rstats.github.io/docs/basics/compare_props.html 

p1 <- ggplot(combined, aes(fill = factor(combined$Celltype, levels = rev(plot_order)), y = freq, x = platform)) +
  geom_bar(position = "stack", stat = "identity") + labs(fill = "Celltype") +
  facet_grid(~SampleID) +  scale_fill_manual(values = pal) + theme_classic()

p1
pdf("/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallFlow/3_Output/FlowVSCiteGatingThymus.pdf")
p1
dev.off()


## plotting for patients
p2 <- ggplot(combined_pts, aes(fill = factor(combined_pts$Celltype, levels = rev(plot_order)), y = freq, x = platform)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(~SampleID) +  scale_fill_manual(values = pal) + theme_classic() + labs(fill = "Celltype")

p2
pdf("/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallFlow/3_Output/FlowVSCiteGatingPatients_V3.pdf")
p2
dev.off()

