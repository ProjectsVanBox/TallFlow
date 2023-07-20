#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape)
library(RColorBrewer)
library(gridExtra)
library(egg)
library(randomcoloR)

#working directory
setwd("~/hpc/projects/TallFlow/3_Output/MIXCR4.3.2/Output/")

#sample name
sample_name = "pt2283D"
dir <- getwd()
Outputdir = paste0(dir, "TCR_Output_", sample_name, "/")
if (!dir.exists(Outputdir)) {dir.create(Outputdir)}   # if permission denied ask Rico to give you permission

# alternative output directory
Outputdir <- "/Users/verapoort/Desktop/output_TCR/" # change if you cannot wait for Rico ;)

# color pallet
n <- 100
palette <- distinctColorPalette(n)


#get dir/file names
nonempty <- list.dirs(recursive=F)[which(lengths(lapply(list.dirs(recursive=F), list.files, pattern = "\\.tsv$")) > 1)]
subdir <- nonempty[grepl(sample_name, nonempty) ]
# exclude directories with no .tsv output files


ff<-do.call(rbind, lapply(subdir, function(x) {
  ff<-list.files(x, "\\.tsv$", include.dirs = FALSE, full.names = TRUE )
  data.frame(dir=basename(x), file=basename(ff), 
             fullpath=ff, stringsAsFactors=F)
}))

tsv_files_list <- lapply(ff$fullpath, function(x) {read.delim(x, header = T, sep ="\t")})
names(tsv_files_list) <- gsub(".*M-","",ff$dir)

# Combine them
combined_df <- do.call("rbind", lapply(tsv_files_list, as.data.frame)) 
combined_df$SampleId <- rep(names(tsv_files_list), sapply(tsv_files_list, nrow))

# remove columns that are not of interest <- not needed but makes checking easier
combined_df <- combined_df[,c("readCount","readFraction","allVHitsWithScore",
               "allDHitsWithScore","allJHitsWithScore",
               "allCHitsWithScore","nSeqCDR3", "aaSeqCDR3","SampleId")]
#remove scores after TR genes
combined_df[,c("allVHitsWithScore",
               "allDHitsWithScore","allJHitsWithScore",
               "allCHitsWithScore")] <- apply(combined_df[,c("allVHitsWithScore",
                     "allDHitsWithScore","allJHitsWithScore",
                     "allCHitsWithScore")], 2, function(y) gsub( x = gsub(pattern = "\\*00\\(.*[0-9]+\\),",replacement =  ",", x = y),
                                                                 pattern = "\\*00.*", replacement = ""))

#merge columns to TCR recombination names, there are never C hits in my data so I do not merge this column
combined_df <- mutate(combined_df, TCR = paste(allVHitsWithScore,
                                                     allDHitsWithScore,allJHitsWithScore, sep = "_"))

# mutate column for gene segment name
combined_df <- mutate(combined_df, segment = case_when(grepl("TRGV", allVHitsWithScore) ~ "Gamma", 
                      grepl("TRBV", allVHitsWithScore) ~ "Beta", 
                      grepl("TRAV", allVHitsWithScore) ~ "Alpha",
                      grepl("TRGV", allVHitsWithScore) ~ "Gamma",
                      grepl("TRDV", allVHitsWithScore) ~ "Delta",
                      grepl("IGKV", allVHitsWithScore) ~ "IGK",
                      grepl("IGLV", allVHitsWithScore) ~ "IGL",
                      grepl("IGHV", allVHitsWithScore) ~ "IGH"))

# Save tables
write.table(combined_df, file = paste0(Outputdir, sample_name, "_combinedTCR.txt"), sep = "\t",
            row.names = TRUE, col.names = TRUE)


# if want to start from saved table: 
combined_df <- read.csv2("Desktop/output_TCR/pt2229_combinedTCR.txt", sep = "\t")
combined_df <- transform(combined_df, cloneFraction = as.numeric(readFraction))
combined_df <- transform(combined_df, cloneCount = as.numeric(readCount))

# remove rearangements that are noise -> in this data set cloneCount of 2 or lower
#filtered_df <- combined_df[combined_df$cloneFraction > 0.05 && combined_df$segment != "Alpha",]

#filtered_df <- filter(combined_df, segment != "Alpha") %>% combined_df[combined_df$cloneFraction > 0.05, ]

# filter segments that you don't want:
#for now only the light chains are filtered out. 
filtered_df <- combined_df %>% filter(!segment %in% c("IGL")) %>% 
  filter(!segment %in% c("IGK")) %>% 
  filter(!segment %in% c("IGH"))

# if no filter used 
filtered_df <- combined_df

# create data set per recombination
Grecomb <- filtered_df %>% filter(str_detect(allVHitsWithScore, "TRGV"))
Arecomb <- filtered_df %>% filter(str_detect(allVHitsWithScore, "TRAV"))
Brecomb <- filtered_df %>% filter(str_detect(allVHitsWithScore, "TRBV"))
Drecomb <- filtered_df %>% filter(str_detect(allVHitsWithScore, "TRDV"))

# the TCR alpha are diluted over all the cells in the DP populations so this should not be filtered
#Arecomb <- combined_df %>% filter(str_detect(allVHitsWithScore, "TRAV"))

#unique(filtered_df$segment)
# layouting for plots
#unique(filtered_df$SampleId)
#pt2229  <- Diagnosis + single cells
SampleOrder <- c("ALLBULK",  "ALLBULK-DNCD1aNeg","ALLBULK-DNCD1aPos", "ALLBULK-iSPCD4", "ALLBULK-DP",
                 "TALL1", "TALL2", "TALL3" ,"TALL4", "TALL5" , "TALL6", "TALL7", "TALL8", "TALL9", "TALL10")
#pt2322  <- Diagnosis + single cells
SampleOrder <- c("ALLBULK",  "ALLBULK-DN","ALLBULK-SPCD4", "ALLBULK-DP", "ALLBULK-iSPCD4",
                "TALL1", "TALL2", "TALL3" ,"TALL4", "TALL5" , "TALL6", "TALL7", 
                "TALL8", "TALL9", "TALL10", "TALL11", "TALL12")

#pt2283D <- Diagnosis + single cells
SampleOrder <- c("ALLBULK",  "ALLBULK-DN","ALLBULK-iSPCD4", "ALLBULK-DP", "ALLBULK-SPCD4",
                 "TALL1", "TALL2", "TALL3" ,"TALL4", "TALL5" , "TALL6", "TALL7", 
                 "TALL8", "TALL9", "TALL10", "TALL11", "TALL12")
#pt2283R <- Relapse only
SampleOrder <- c("ALLBULK",  "ALLBULK-DN","ALLBULK-iSPCD4", "ALLBULK-DP", "ALLBULK-SPCD4")

            
# plot  
## plotting cloneFraction or cloneCount makes no difference
TRG_plot <-ggplot(Grecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_brewer(palette = "YlOrBr", direction = -1) + 
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

TRG_plot_stacked <-ggplot(Grecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_brewer(palette = "YlOrBr", direction = -1) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

TRA_plot <- ggplot(Arecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = palette) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

TRA_plot_stacked <- ggplot(Arecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = palette) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

TRB_plot <- ggplot(Brecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_brewer(palette = "YlGnBu") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

TRB_plot_stacked <- ggplot(Brecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_brewer(palette = "YlGnBu") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

TRD_plot <- ggplot(Drecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_brewer(palette = "YlGn") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

TRD_plot_stacked <- ggplot(Drecomb, aes(fill=TCR, y=readFraction, x=SampleId)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_brewer(palette = "YlGn") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#print PDFs with plots
pdf(paste0(Outputdir, sample_name, "_rel_fixed.pdf"))
grid.arrange(grobs=lapply(list(
  TRG_plot, TRA_plot),
  set_panel_size), ncol = 1)
grid.arrange(grobs=lapply(list(
  TRB_plot, TRD_plot),
  set_panel_size), ncol = 1)
dev.off()


pdf(paste0(Outputdir, sample_name,"_Abs_fixed.pdf"))
grid.arrange(grobs=lapply(list(
  TRG_plot_stacked, TRA_plot_stacked),
  set_panel_size), ncol = 1)
grid.arrange(grobs=lapply(list(
  TRB_plot_stacked, TRD_plot_stacked),
  set_panel_size), ncol = 1)
dev.off()

pdf(paste0(Outputdir, sample_name, "_Alpha.pdf"))
TRA_plot
dev.off()

pdf(paste0(Outputdir, sample_name, "_Alpha_abs.pdf"))
grid.arrange(grobs=lapply(list(TRA_plot_stacked),
                          set_panel_size), ncol = 1)
dev.off()

pdf(paste0(Outputdir, sample_name,"_TRB_abs.pdf"))
grid.arrange(grobs=lapply(list(TRB_plot_stacked),
  set_panel_size), ncol = 1)
dev.off()

pdf(paste0(Outputdir, sample_name, "_TRD_abs.pdf"))
grid.arrange(grobs=lapply(list(TRD_plot_stacked),
                          set_panel_size), ncol = 1)
dev.off()



# ---- colors for next plots
colors <- c( Delta = "#ADD18D", Gamma = "#FCC34F", Beta = "cornflowerblue",
             Alpha = "#C181B6", IGK = "lightgrey", IGL = "darkgrey", IGH = "grey")


# --------------- Plot which gene segments are used

segement_contr <- ggplot(combined_df, aes(fill=segment, y=readFraction, x=SampleId)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = colors) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 2, ncol = 1)))+ 
  scale_x_discrete(limits=SampleOrder) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste0(Outputdir, "segmentcount_", sample_name, ".pdf"))
segement_contr
dev.off()
    

## plot donut
# more color tweaking: https://stackoverflow.com/questions/68095243/piedonut-how-to-change-color-of-pie-and-donut 
library(webr)

for (sample in unique(filtered_df$SampleId)){
  donut_df <- filtered_df[filtered_df$SampleId == sample,] %>% group_by(TCR, segment) %>% summarise(n = n())
  pdf_fname <- paste0(sample,"_donut_filtered.pdf")
  #pdf(pdf_fname)
  print(PieDonut(donut_df, aes(segment, TCR, count=n)))
  #dev.off()
}



### alternative option:
#install.packages("ggpie")
library(ggpie)
ggnestedpie(data = combined_df[combined_df$SampleId == sample,], group_key = c("TCR", "segment"), count_type = "full",
            inner_label_info = "all", inner_label_split = NULL,inner_label_size = 2,
            outer_label_type = "circle", outer_label_pos = "in", outer_label_info = "all")

