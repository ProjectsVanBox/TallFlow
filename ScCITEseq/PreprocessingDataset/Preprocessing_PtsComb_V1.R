# Load packages
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(DBI)
library(AnnotationDbi)
library(SingleCellExperiment)
library(tidyverse)
library(scales)
library(cowplot)
library(RCurl)
library(org.Hs.eg.db)
library(SeuratDisk)
library(matrixStats)
library(ggsci)
library(SingleR)
library(batchelor)
library(BiocParallel)
library(BiocNeighbors)
library(celldex)


### Change directory 
setwd("/hpc/pmc_vanboxtel/projects/TallClonal/2_Code/PreprocessingDataset/")
#setwd("/Users/ricohagelaar/Documents/Thymus_Project/SingleCellSequencing/TallClonal/2_Code/PreprocessingDataset/")


### Prepare data structure and read in info files 
inputdir = "../../1_Input/10x_Pts/"
outputdir = "../../3_Output/10x_PtsComb/Preprocessing/"
plotdir = paste0(outputdir, "Plots/")
project = "Pts_Combined"
# Metadata
meta_data = as.data.frame(readxl::read_excel(paste0(inputdir,"CohortInfoScStudy.xlsx")))


### Colors and names for coming analysis 
colors = list(SampleID = c("pt1179" = "#339999", "pt2337" = "#124D4D",
                           "pt344" = "#CC3366", "pt3291" = "#801438",
                           "pt4068" = "#99CC66", "pt5242" = "#44631F",
                           "pt2283" = "#FF9933", "pt11801" = "#8E4503",
                           "pt9160" = "#FFCC00", "pt315" = "#B29330",
                           "pt335" = "#003366", "pt9175" = "#7CE3D8",
                           "pt3045" = "#b0c4de", "pt5438" = "#6495ed",
                           "pt5676" = "#875692", "pt10138" = "#006666",
                           "Thy01" = "#ffa07a", "Thy02" = "#C2B280", "Thy03" = "#C71585"),
              PatientID = c("TALL1" = "#339999", "TALL2" = "#CC3366", 
                            "TALL3" = "#99CC66","TALL4" = "#FF9933", 
                            "TALL5" = "#FFCC00", "TALL6" = "#003366", 
                            "TALL7" = "#7CE3D8", "TALL8" = "#b0c4de",
                            "TALL9" = "#6495ed", "TALL10" = "#875692", 
                            "TALL11" = "#006666", "Thy01" = "ffa07a", 
                            "Thy02" = "#C2B280", "Thy03" = "C71585"),
              Site = c("BM" = "43C5F9", "PB" = "#FF0000"),
              Gender = c("M" = "#003366", "F" = "#CC3366"),
              Type = c("TALL" = "#052955", "THYMUS" = "#C59434"),
              Stage = c("D" = "#C59434", "R" = "#A3B7F9"),
              Clusters = c("#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5", 
                           "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", 
                           
                           "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", "#99600F", "#B3823E", "#CCAA7A", "#E6D2B8", "#54990F", "#78B33E", 
                           "#A3CC7A", "#CFE6B8", "#0F8299", "#3E9FB3", "#7ABECC", "#B8DEE6", "#3D0F99", "#653EB3", "#967ACC", "#C7B8E6", 
                           "#333333", "#666666", "#999999", "#CCCCCC"))


### Create and store count matrix 
# In inputdir store cell ranger output filtered folder with count matrix files 
# Rename folder to pt"1234"  this name will be use in the analysis 
mtx_dirs = list.dirs(path= inputdir, recursive = F, full.names = TRUE)
count_files <- list()
srat_objects <- list()
for( i in 1:length(mtx_dirs)){
  filename <-gsub(paste0(inputdir, "/"), "", mtx_dirs[i])
  count_files[[filename]] <- append(count_files, Read10X(mtx_dirs[i]))
  name <- gsub(paste0(inputdir, "/"), "", mtx_dirs[i])
  srat_objects[[name]] <- CreateSeuratObject(counts = count_files[[i]]$`Gene Expression`, 
                                             assay = "RNA", strip.suffix = T)
  srat_objects[[name]][["Protein"]] <- CreateAssayObject(counts = count_files[[i]]$`Antibody Capture`)
}


### Merge data sets to srat 
cell_id <- as.list(names(srat_objects))
srat_all <- merge(x = srat_objects[[1]], y = srat_objects[-1], 
                  add.cell.ids = cell_id, project = project)
srat_all$orig.ident <- sapply(X = strsplit(colnames(srat_all), split = "_"), FUN = "[", 1)
srat_all$SampleID <- sapply(X = strsplit(colnames(srat_all), split = "_"), FUN = "[", 1)
srat = srat_all

# Clean up global env
rm(count_files)
rm(srat_objects)
rm(srat_all)
gc()


# Add meta data
m <- match(srat$SampleID, meta_data$Sample_ID)
srat$PatientID = meta_data$Patient_ID[m]
srat$Stage = meta_data$Stage[m]
srat$Site = meta_data$Site[m]
srat$Gender = meta_data$Gender[m]
srat$Type = meta_data$Type[m]


### Mark stress or contaminating genes 
# Mitochondrial stress
srat$percent.mt <- PercentageFeatureSet(object = srat, pattern = "^MT-")
# UPR stress
GO_genes <- read.csv(file = paste0(inputdir, "GO_UPR.csv"), sep = ";")
GO_gene_list <- unique(GO_genes$Symbol)
GO_gene_list <- GO_gene_list[GO_gene_list %in% row.names(srat)]   
srat$percent.UPR <- PercentageFeatureSet(object = srat, features = GO_gene_list, assay = "RNA")
# Likely some RBCs in the data set with low n genes but high n transcripts
hb_genes <- rownames(srat[["RNA"]]@data)[grep(pattern = "^HB", rownames(srat[["RNA"]]@data))]
hb_counts <- Matrix::colSums(srat@assays$RNA@counts[hb_genes,])
srat <- AddMetaData(srat, col.name = "log2hb_genes", metadata = log2(1 +hb_counts))
srat <- AddMetaData(srat, col.name = "pct_hemo", metadata = 100 * hb_counts/srat@meta.data$nCount_RNA)


# Give an identifier
Idents(srat) <- "SampleID"


### Plot quality features 
p1 <- plot_grid(VlnPlot(srat, features=c('nCount_RNA'), pt.size = 0, y.max = 30000, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"), 
                VlnPlot(srat, features=c('nFeature_RNA'), pt.size = 0, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),
                VlnPlot(srat, features=c('percent.mt'), pt.size = 0, y.max = 40, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),
                VlnPlot(srat, features=c('percent.UPR'), pt.size = 0, y.max = 10, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),  
                VlnPlot(srat, features=c('nCount_Protein'), pt.size = 0, y.max = 25000,cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"),
                VlnPlot(srat, features=c('nFeature_Protein'), pt.size = 0, y.max = 220, cols = colors$SampleID)+
                  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                        axis.title.x = element_blank()) + guides(fill="none") + 
                  theme(axis.text.x = element_text(size = 8)) + xlab("SampleID"))
png(paste0(plotdir,"Features.png"))
p1 
dev.off()


### Plot quality 
data = srat@meta.data
# Percent mitochrondrial
p2 = ggplot(data, aes(x=nCount_RNA, y= percent.mt))+
  geom_point(shape=16, size=1)+
  scale_fill_igv()+
  coord_trans(x = "log10")+
  scale_y_continuous(breaks= scales::pretty_breaks(n=10))+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_light()
# nFeature_RNA
p3 = ggplot(data, aes(x=nCount_RNA, y= nFeature_RNA))+
  scale_fill_igv()+
  geom_point(shape=16, size=1)+
  coord_trans(x = "log10", y = "log10")+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_light()
# Plot
png(paste0(plotdir,"Quality.png"))
plot_grid(p2, p3) 
dev.off()
# Percentage hemo
v <- VlnPlot(srat, features = "pct_hemo", cols = colors$SampleID) & NoLegend()


### Filter low quality cells
srat$drop = !c(srat$nCount_RNA <= 500 | srat$nFeature_RNA <= 250 | srat$percent.mt >= 30
               | srat$percent.UPR >= 4 | srat$nCount_RNA >= 40000 | srat$nFeature_RNA >= 7000
               | srat$nCount_Protein >= 20000 | srat$pct_hemo >= 2 | srat$nFeature_Protein <= 125 )
#| srat$nFeature_Protein >= 125
data = srat@meta.data
# Percentage mitochorndrial 
p4 = ggplot(data, aes(x=nCount_RNA, y= percent.mt, color = drop))+
  geom_point(shape=16, size=1)+
  scale_fill_igv()+
  coord_trans(x = "log10")+
  scale_y_continuous(breaks= scales::pretty_breaks(n=10))+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_color_manual(values=c("red", "black")) +
  theme_light() & NoLegend()
# nFeature_RNA
p5 = ggplot(data, aes(x=nCount_RNA, y= nFeature_RNA, color = drop))+
  scale_fill_igv()+
  geom_point(shape=16, size=1)+
  coord_trans(x = "log10", y = "log10")+
  scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_color_manual(values=c("red", "black")) +
  theme_light() & NoLegend()
# Plot 
png(paste0(plotdir,"Filtered.png"))
plot_grid(p4, p5) 
dev.off()


### Filter out poor quality 
srat = subset(srat, subset = drop)


### SCT normalisation 
srat <- SCTransform(srat, vars.to.regres = NULL, verbose = FALSE)


### Filter cell cycle-,  mitochrondrial-, ribosomal-/heatshock-, gender-/sex-genes 
# Cell cycle genes 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, min.cells = 1)
prop.table(x = table(srat@meta.data$Phase))
cell_cycle_genes <- na.omit(unique(c(s.genes, g2m.genes)))
cell_cycle_genes

# Find genes which correlate strongly with cell cycle score, but do not filter yet (check if you miss biological important genes)
#cor_mat <- cor(cbind("S.score" = srat@meta.data$S.Score, "G2M.score" = srat@meta.data$G2M.Score), t(as.matrix(srat[["RNA"]]@data)))
#cor_mat <- cor_mat[ , !(is.na(colSums(cor_mat)))]
#S_S = cor_mat[ "S.score", colnames(cor_mat) %in% s.genes]
#G2M_G2M = cor_mat[ "G2M.score", colnames(cor_mat) %in% g2m.genes]
#S.cutoff <- quantile(S_S, 0.25)
#G2M.cutoff <- quantile(G2M_G2M, 0.25)
#S.genes.derived <- colnames(cor_mat)[cor_mat["S.score",] >= S.cutoff & cor_mat["S.score",] - cor_mat["G2M.score",] >= 0.05]
#G2M.genes.derived <- colnames(cor_mat)[cor_mat["G2M.score",] >= G2M.cutoff & cor_mat["G2M.score",] - cor_mat["S.score",] >= 0.05] 
# Save the derived/correlated genes
#write.table(S.genes.derived, file = paste0(outputdir, "/CC_CorGenes/S_genes_derived.txt"), sep = "\t", col.names = NA)
#write.table(G2M.genes.derived, file = paste0(outputdir, "/CC_CorGenes/G2M_genes_derived.txt"), sep = "\t", col.names = NA)

# Mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = srat[["RNA"]]@data), value = TRUE)

# Unfolder protein (UPR)
GO_UPR <- read.csv(file = paste0(inputdir,"GO_UPR.csv"), sep = ";")
GO_UPR_list <- unique(GO_UPR$Symbol)
GO_UPR_list <- GO_UPR_list[GO_UPR_list %in% row.names(srat)]

# Sex genes 
gencode_anno <- readRDS("../GENCODE26_anno_PMC.RDS")
genes_X <- c("ENSG00000229807--XIST", "ENSG00000270641--TSIX")
ENS_Y <- gencode_anno$gene_id[gencode_anno$seqnames == "chrY"]
genes_Y <- as.vector(apply(gencode_anno[gencode_anno$gene_id %in% ENS_Y, c("gene_id", "gene_name")], 1, paste, collapse = "--"))
genes_Y <- sub(".*--", "", genes_Y)
genes_X <- sub(".*--", "", genes_X)
sex_genes <- unique(c(genes_X, genes_Y))

# Function to get genes related to GO terms
reToTags = function (retok, ...) {
  require(GO.db)
  allt = dbGetQuery(GO_dbconn(), "select go_id, term from go_term")
  inds = grep(retok, as.character(allt[,"term"]), value=FALSE, ...)
  # will error if value set at call
  if (length(inds)>0) return(as.character(allt[inds,"go_id"]))
  stop("retok did not grep to any term")
}
# Get genes related to unfolded protein and the ribosome
GO_df_all <- NULL
GOs <- c("^response to unfolded protein$", "cytosolic ribosome")
for (i in 1:length(GOs)) {
  GOTags <- reToTags(GOs[i]); cat(GOTags)
  GO_df <- AnnotationDbi::select(org.Hs.eg.db, keys = GOTags,
                                 columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), keytype = "GOALL")
  GO_df_all <- rbind(GO_df_all, GO_df)
}
GO_df_all$gene_id_name <- paste(GO_df_all$ENSEMBL, GO_df_all$SYMBOL, sep = "--")
#filter_GO = GO_df_all$gene_id_name[GO_df_all$gene_id_name %in% rownames(srat)]
filter_GO = GO_df_all$SYMBOL[GO_df_all$SYMBOL %in% rownames(srat)]


### Filtering part 
all_exclud_genes <- c(cell_cycle_genes, mito.genes, GO_UPR_list, sex_genes, filter_GO)
filter_g <- VariableFeatures(srat) %in% unique(all_exclud_genes)
summary(filter_g)
VariableFeatures(srat) <- VariableFeatures(srat)[!filter_g]


### Save output for future steps 
saveRDS(srat, paste0(outputdir,"srat_10xPtsComb_Preprocessing.rds"))



