#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(MutationalPatterns)
library(BSgenome)
library(ChIPpeakAnno)
library(ggplot2)
library(NMF)
library(RColorBrewer)
library(tibble)
library(reshape2)
library(grid)
mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2 <- brewer.pal(8, "Dark2")
### Read in genomes
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
hg38_autosomal_nonN_genome_size <- 2745186691


## Plot mutation spectra as stacked bar
plot_spectrum_stacked <- function(type_occurrences, mode = "absolute", max_value_to_show = 30, facet = "", total_value = TRUE, round_number = 0, title = ""){
  type_occurrences_stack <- type_occurrences[,c(1,2, 7, 8,4,5,6)]
  ylabel <- "Absolute contribution"
  y_extension <- 0.1
  if(mode == "relative"){
    type_occurrences_stack <- type_occurrences_stack / rowSums(type_occurrences_stack)
    round_number <- 2
    ylabel <- "Relative contribution"
    y_extension <- 0
  }
  
  type_occurrences_stack$Sample <- row.names(type_occurrences_stack)
  type_occurrences_stack$Sample <- factor(type_occurrences_stack$Sample, levels = row.names(type_occurrences_stack))
  if(length(facet > 1)){
    type_occurrences_stack$Facet <- facet
    type_occurrences_stack_m <- melt(type_occurrences_stack, id.vars = c("Sample", "Facet"))
    
  } else {
    type_occurrences_stack_m <- melt(type_occurrences_stack, id.vars = "Sample")
    
  }
  type_occurrences_stack_m$Label <- type_occurrences_stack_m$value
  # The values below max_value_to_show are not plotted
  type_occurrences_stack_m$Label[type_occurrences_stack_m$Label < max_value_to_show] <- ""
  # type_occurrences_total <- aggregate(type_occurrences_stack_m$value ~ type_occurrences_stack_m$Sample, FUN = "sum")
  # colnames(type_occurrences_total) <- c("Sample", "Count")
  if(length(facet) < 2){
    
    plot <- ggplot(data = type_occurrences_stack_m, aes(x = Sample, 
                                                        y = value, 
                                                        label = as.character(round(as.numeric(Label),round_number )))) + 
      geom_bar(stat = "identity", aes(fill = variable), width = 0.8) +
      geom_text(aes(fill = variable), size = 2, position = position_stack(vjust = 0.5)) +
      #facet_grid(.~Type, space = "free",scales = "free", switch = "both") +
      scale_fill_manual(values = COLORS7) +
      scale_y_continuous(expand = expansion(mult = c(0,y_extension))) +
      theme_classic() +
      ggtitle(title) + 
      labs(x = "Sample", y = ylabel, fill = "Point mutation type")  +
      theme(axis.text.x= element_text(angle = 45, hjust = 1), 
            strip.placement = "outside", 
            strip.background = element_rect(fill = "lightgray"), 
            plot.title = element_text(hjust = 0.5))
    if(total_value == TRUE){
      plot +
        geom_text(aes(label = round(stat(y), round_number), group = Sample), 
                  stat = 'summary', fun = sum, vjust = -1, size = 3) 
    } else {
      plot
    }
  } else {
    ggplot(data = type_occurrences_stack_m, aes(x = Sample, y = value, label = as.character(round(as.numeric(Label),round_number )))) + 
      geom_bar(stat = "identity", aes(fill = variable), width = 0.8) +
      geom_text(aes(fill = variable), size = 2, position = position_stack(vjust = 0.5)) +
      geom_text(aes(label = round(stat(y), 0), group = Sample), 
                stat = 'summary', fun = sum, vjust = -1, size = 3) + 
      facet_grid(.~Facet, scales = "free", space = "free") +
      #facet_grid(.~Type, space = "free",scales = "free", switch = "both") +
      scale_fill_manual(values = COLORS7) +
      scale_y_continuous(expand = expansion(mult = c(0,y_extension))) +
      theme_classic() +
      labs(x = "Sample", y = ylabel, fill = "Point mutation type") +
      ggtitle(title) + 
      theme(axis.text.x= element_text(angle = 45, hjust = 1), 
            strip.placement = "outside", 
            strip.background = element_rect(fill = "lightgray"), 
            plot.title = element_text(hjust = 0.5))
  }
  
}


### Read in with overview of number mutations 
# pt2322
pt2322_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/pt2322_ageline_table.txt", sep = " ", header = T)
pt2322_data$age <- rep.int(x = 9.4, times = length(rownames(pt2322_data)))
pt2322_data$nmuts_norm <- (pt2322_data$SNV_LOAD/pt2322_data$CALLABLE)*hg38_autosomal_nonN_genome_size
pt2322_data$nmuts_norm2 <- pt2322_data$SNV_LOAD*(pt2322_data$CALLABLE/hg38_autosomal_nonN_genome_size)
pt2322_data$nmuts_norm3 <- pt2322_data$SNV_LOAD*(1/(pt2322_data$CALLABLE/hg38_autosomal_nonN_genome_size))
# pt2229
pt2229_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/pt2229_ageline_table.txt", sep = " ", header = T)
pt2229_data$age <- rep.int(x = 5.5, times = length(rownames(pt2229_data)))
pt2229_data$nmuts_norm <- (pt2229_data$SNV_LOAD/pt2229_data$CALLABLE)*hg38_autosomal_nonN_genome_size
# pt2283D
pt2283D_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/pt2283D_ageline_table.txt", sep = " ", header = T)
pt2283D_data$age <- rep.int(x = 10, times = length(rownames(pt2283D_data)))
pt2283D_data$nmuts_norm <- (pt2283D_data$SNV_LOAD/pt2283D_data$CALLABLE)*hg38_autosomal_nonN_genome_size


### Read in vcf files 
# Files 
pt2322_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2322/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt2229_vcf_files <- list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/PTATO/POTATO_Dev/snvs/pt2229/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
pt2283D_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/snvs/pt2283D/",
                               pattern = "*filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
# Sample names 
pt2322_sample_names <- c("pt2322_ALLBULK", "pt2322_TALL1", "pt2322_TALL10", "pt2322_TALL11", "pt2322_TALL12", "pt2322_TALL2", "pt2322_TALL3", "pt2322_TALL4",
                         "pt2322_TALL5", "pt2322_TALL6", "pt2322_TALL7", "pt2322_TALL8", "pt2322_TALL9") 
pt2229_sample_names <- c("pt2229_ALLBULK", "pt2229_TALL1", "pt2229_TALL10", "pt2229_TALL2", "pt2229_TALL3", "pt2229_TALL4", "pt2229_TALL5", "pt2229_TALL6", 
                         "pt2229_TALL7", "pt2229_TALL8", "pt2229_TALL9")
pt2283D_sample_names <- c("pt2283D_ALLBULK", "pt2283D_TALL1", "pt2283D_TALL10", "pt2283D_TALL11", "pt2283D_TALL12", "pt2283D_TALL2", "pt2283D_TALL3", "pt2283D_TALL4", "pt2283D_TALL5", "pt2283D_TALL6", 
                          "pt2283D_TALL7", "pt2283D_TALL8", "pt2283D_TALL9")
# Create matrices 
pt2322_grl <- read_vcfs_as_granges(pt2322_vcf_files, pt2322_sample_names, ref_genome)
pt2322_mut_mat <- mut_matrix(vcf_list = pt2322_grl, ref_genome = ref_genome)
pt2229_grl <- read_vcfs_as_granges(pt2229_vcf_files, pt2229_sample_names, ref_genome)
pt2229_mut_mat <- mut_matrix(vcf_list = pt2229_grl, ref_genome = ref_genome)
pt2283D_grl <- read_vcfs_as_granges(pt2283D_vcf_files, pt2283D_sample_names, ref_genome)
pt2283D_mut_mat <- mut_matrix(vcf_list = pt2283D_grl, ref_genome = ref_genome)


#grl_list <- c(pt2322_grl, pt2322_branch_grl)#, pt2229_mut_mat, pooled_pt2229_branch_mat)
pt2322_type_occurrences <- mut_type_occurrences(pt2322_grl, ref_genome)
pt2322_type_occurrences
pt2229_type_occurrences <- mut_type_occurrences(pt2229_grl, ref_genome)
pt2229_type_occurrences
pt2283D_type_occurrences <- mut_type_occurrences(pt2283D_grl, ref_genome)
pt2283D_type_occurrences
COLORS6 <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
)
COLORS7 <- c(
  "#2EBAED", "#000000", 
  "#E98C7B", "#DE1C14","#D4D2D2", "#ADCC54",
  "#F0D0CE"
)
plot_spectrum_stacked(pt2322_type_occurrences, max_value_to_show = 10)


hg38_autosomal_nonN_genome_size <- 2745186691
pt2322_data$mut_rate <- pt2322_data$CALLABLE / hg38_autosomal_nonN_genome_size
rownames(pt2322_data) <- gsub(x = pt2322_data$SAMPLE, pattern = "-DX1BM-", replacement = "_")
pt2322_data <- pt2322_data[rownames(pt2322_type_occurrences),]
pt2322_type_occurrences * pt2322_data$mut_rate
#pt2322_type_occurrences2 <- pt2322_type_occurrences * pt2322_data$mut_rate
pt2322_type_occurrences2 <- pt2322_type_occurrences * (1 / pt2322_data$mut_rate)
pt2322_type_occurrences2[1,] <- pt2322_type_occurrences[1,]

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/NormNumbers/Original_pt2322.pdf")
plot_spectrum_stacked(pt2322_type_occurrences, max_value_to_show = 10)
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/NormNumbers/Norm_pt2322.pdf")
plot_spectrum_stacked(pt2322_type_occurrences2, max_value_to_show = 10)
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/NormNumbers/Norm2_pt2322.pdf")
plot_spectrum_stacked(pt2322_type_occurrences2, max_value_to_show = 10)
dev.off()


pt2229_data$mut_rate <- pt2229_data$CALLABLE / hg38_autosomal_nonN_genome_size
rownames(pt2229_data) <- gsub(x = pt2229_data$SAMPLE, pattern = "-DX1BM-", replacement = "_")
pt2229_data <- pt2229_data[rownames(pt2229_type_occurrences),]
pt2229_type_occurrences * pt2229_data$mut_rate
pt2229_type_occurrences2 <- pt2229_type_occurrences * (1 / pt2229_data$mut_rate)
pt2229_type_occurrences2[1,] <- pt2229_type_occurrences[1,]

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/NormNumbers/Original_pt2229.pdf")
plot_spectrum_stacked(pt2229_type_occurrences, max_value_to_show = 10)
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/NormNumbers/Norm_pt2229.pdf")
plot_spectrum_stacked(pt2229_type_occurrences2, max_value_to_show = 10)
dev.off()


pt2283D_data$mut_rate <- pt2283D_data$CALLABLE / hg38_autosomal_nonN_genome_size
rownames(pt2283D_data) <- gsub(x = pt2283D_data$SAMPLE, pattern = "-DX1BM-", replacement = "_")
pt2283D_data <- pt2283D_data[rownames(pt2283D_type_occurrences),]
pt2283D_type_occurrences * pt2283D_data$mut_rate
pt2283D_type_occurrences2 <- pt2283D_type_occurrences * (1 / pt2283D_data$mut_rate)
pt2283D_type_occurrences2[1,] <- pt2283D_type_occurrences[1,]

pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/NormNumbers/Original_pt2283D.pdf")
plot_spectrum_stacked(pt2283D_type_occurrences, max_value_to_show = 10)
dev.off()
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/NormNumbers/Norm_pt2283D.pdf")
plot_spectrum_stacked(pt2283D_type_occurrences2, max_value_to_show = 10)
dev.off()


pt2322_data[is.na(pt2322_data)] <- 1
pt2322_data[1,] <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
rownames(pt2322_data[1,]) <- c("pt2322_ALLBULK")


input_df$INDEL_LOAD_NORM <- input_df$INDEL_LOAD/input_df$CALLABLE*hg38_autosomal_nonN_genome_size






# pt2322
pt2322_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/Bulk_pt2322_ageline_table.txt", sep = " ", header = T)
pt2322_data$age <- rep.int(x = 9.4, times = length(rownames(pt2322_data)))
pt2322_data$mut_rate <- pt2322_data$CALLABLE / hg38_autosomal_nonN_genome_size
rownames(pt2322_data) <- gsub(x = pt2322_data$SAMPLE, pattern = "-DX1BM-", replacement = "-")
rownames(pt2322_data) <- gsub(x = rownames(pt2322_data), pattern = "-ALLBULK-", replacement = "-")
# pt2229
pt2229_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/Bulk_pt2229_ageline_table.txt", sep = " ", header = T)
pt2229_data$age <- rep.int(x = 5.5, times = length(rownames(pt2229_data)))
pt2229_data$mut_rate <- pt2229_data$CALLABLE / hg38_autosomal_nonN_genome_size
rownames(pt2229_data) <- gsub(x = pt2229_data$SAMPLE, pattern = "-DX1BM-", replacement = "-")
rownames(pt2229_data) <- gsub(x = rownames(pt2229_data), pattern = "-ALLBULK-", replacement = "-")
pt2229_data
# pt2283D
pt2283D_data <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/2_Code/AgeLine/Bulk_pt2283D_ageline_table.txt", sep = " ", header = T)
pt2283D_data$age <- rep.int(x = 10, times = length(rownames(pt2283D_data)))
pt2283D_data$mut_rate <- pt2283D_data$CALLABLE / hg38_autosomal_nonN_genome_size
rownames(pt2283D_data) <- gsub(x = pt2283D_data$SAMPLE, pattern = "-DX1BM-", replacement = "-")
rownames(pt2283D_data) <- gsub(x = rownames(pt2283D_data), pattern = "-ALLBULK-", replacement = "-")






### check Bulk sample 
# pt2322
pt2322_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2322/",
                               pattern = "*SMuRF.filtered.vcf$", full.names = TRUE, recursive = FALSE )
pt2322_sample_names <- c("pt2322_DN", "pt2322_DP", "pt2322_iSPCD4", "pt2322_SPCD4", "pt2322_ALLBulk", "pt2322_MONO", "pt2322_Combined")
pt2322_grl <- read_vcfs_as_granges(pt2322_vcf_files, pt2322_sample_names, ref_genome)
pt2322_type_occurrences <- mut_type_occurrences(pt2322_grl, ref_genome)
pt2322_type_occurrences <-pt2322_type_occurrences[rownames(pt2322_type_occurrences) %in% c("pt2322_DN", "pt2322_DP", "pt2322_iSPCD4", "pt2322_SPCD4", "pt2322_ALLBulk"),]
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Bulk_StackedBarCharts/WithoutCallable/pt2322.pdf")
plot_spectrum_stacked(pt2322_type_occurrences
                      , max_value_to_show = 10)
dev.off()
pt2322_type_occurrences2 <- pt2322_type_occurrences * (1 / pt2322_data$mut_rate)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Bulk_StackedBarCharts/WithCallable/pt2322.pdf")
plot_spectrum_stacked(pt2322_type_occurrences2
                      , max_value_to_show = 10)
dev.off()
# pt2229
pt2229_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2229inclDP/",
                               pattern = "*SMuRF.filtered.vcf$", full.names = TRUE, recursive = FALSE )
pt2229_sample_names <- c("pt2229-DNCD1aNeg", "pt2229-DNCD1aPos", "pt2229-DP", "pt2229-iSPCD4", "pt2229-ALLBULK", "pt2229-MSCBULK", "pt2229_Combined")
pt2229_grl <- read_vcfs_as_granges(pt2229_vcf_files, pt2229_sample_names, ref_genome)
pt2229_type_occurrences <- mut_type_occurrences(pt2229_grl, ref_genome)
pt2229_type_occurrences <- pt2229_type_occurrences[rownames(pt2229_type_occurrences) %in% c("pt2229-DNCD1aNeg", "pt2229-DNCD1aPos", "pt2229-DP", "pt2229-iSPCD4", "pt2229-ALLBULK"),]
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Bulk_StackedBarCharts/WithoutCallable/pt2229.pdf")
plot_spectrum_stacked(pt2229_type_occurrences
                      , max_value_to_show = 10)
dev.off()
pt2229_type_occurrences2 <- pt2229_type_occurrences * (1 / pt2229_data$mut_rate)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Bulk_StackedBarCharts/WithCallable/pt2229.pdf")
plot_spectrum_stacked(pt2229_type_occurrences2
                      , max_value_to_show = 10)
dev.off()
# pt2283D
pt2283D_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2283D/",
                               pattern = "*SMuRF.filtered.vcf$", full.names = TRUE, recursive = FALSE )
pt2283D_sample_names <- c("pt2283D-DN", "pt2283D-DP", "pt2283D-iSPCD4", "pt2283D-SPCD4", "pt2283D-ALLBULK", "pt2283D-MSCBULK", "pt2283D-Combined")
pt2283D_grl <- read_vcfs_as_granges(pt2283D_vcf_files, pt2283D_sample_names, ref_genome)
pt2283D_type_occurrences <- mut_type_occurrences(pt2283D_grl, ref_genome)
pt2283D_type_occurrences <- pt2283D_type_occurrences[rownames(pt2283D_type_occurrences) %in% c("pt2283D-DN", "pt2283D-DP", "pt2283D-iSPCD4", "pt2283D-SPCD4", "pt2283D-ALLBULK"),]
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Bulk_StackedBarCharts/WithoutCallable/pt2283D.pdf")
plot_spectrum_stacked(pt2283D_type_occurrences
                      , max_value_to_show = 10)
dev.off()
pt2283D_type_occurrences2 <- pt2283D_type_occurrences * (1 / pt2283D_data$mut_rate)
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Bulk_StackedBarCharts/WithCallable/pt2283D.pdf")
plot_spectrum_stacked(pt2283D_type_occurrences2
                      , max_value_to_show = 10)
dev.off()
# pt344
pt344_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt344/",
                                pattern = "*SMuRF.filtered.vcf$", full.names = TRUE, recursive = FALSE )
pt344_sample_names <- c("pt344-ALLBULK", "pt344-DNCD1aPos", "pt344-DP", "pt344-SPCD8", "pt344-MONOBULK", "pt344-Combined")
pt344_grl <- read_vcfs_as_granges(pt344_vcf_files, pt344_sample_names, ref_genome)
pt344_type_occurrences <- mut_type_occurrences(pt344_grl, ref_genome)
pt344_type_occurrences <- pt344_type_occurrences[rownames(pt344_type_occurrences) %in% c("pt344-ALLBULK", "pt344-DNCD1aPos", "pt344-DP", "pt344-SPCD8"),]
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Mutpatterns/Bulk_StackedBarCharts/WithoutCallable/pt344.pdf")
plot_spectrum_stacked(pt344_type_occurrences,
                      max_value_to_show = 10)
dev.off()

# Does this one have way too many muts? 
# pt2283R
#pt2283R_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/BulkSeq_SMuRF/pt2283R/SingleVCFs/",
#                                pattern = "*SMuRF.filtered.vcf$", full.names = TRUE, recursive = FALSE )
#pt2283R_sample_names <- c("pt2283R-DN", "pt2283R-DP", "pt2283R-iSPCD4", "pt2283R-SPCD4", "pt2283R-ALLBULK")
#pt2283R_grl <- read_vcfs_as_granges(pt2283R_vcf_files, pt2283R_sample_names, ref_genome)
#pt2283R_type_occurrences <- mut_type_occurrences(pt2283R_grl, ref_genome)
#plot_spectrum_stacked(pt2283R_type_occurrences, 
#                      max_value_to_show = 10)
