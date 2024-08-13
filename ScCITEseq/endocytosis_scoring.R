# module scoring Endocytosis

## library
library(Seurat)
library(ggplot2)
library(dplyr)

## load data 
srat_pts <- readRDS('srat_pts.rds')
# Endocytosis gene set 
endo_info <- read.delim('../../1_Input/Endocytosis_GEO/KEGG_ENDOCYTOSIS.v2023.2.Hs.tsv', sep ='\t', header = T)
endo_genes <- unlist(strsplit(endo_info$KEGG_ENDOCYTOSIS[17], ","))

## colors 
cell_types_col <- c( "#124D4D", "#222222", "#F3C300", "#875692", "#F38400", "#A1CAF1",
                     "#BE0032", "#C2B280","#F99379", "#008856", "#E68FAC", "#0067A5",
                     "#999999", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#43C5F9", "cornflowerblue",  "#003366")

# compare genelist endocytosis with genes in seurat object
endo_genes_filt <- intersect(endo_genes, rownames(srat_pts$SCT))

## add endocytosis gene score
srat_pts <- AddModuleScore(srat_pts, features = list(endocyt = endo_genes_filt), name = "Endo")

colnames(srat_pts@meta.data)

## plot scores 
plot <- VlnPlot(srat_pts, features = "Endo1", group.by = "gating")

pdf("endocytosisscore.pdf")
print(plot)
dev.off()

## per patient analysis 
df_endo <- srat_pts@meta.data[[c("orig.ident","gating", "Endo1")]]
df_endo$gating <- factor(df_endo$gating, levels= c("DNearly", "DN3", "gdTcell", "iSPCD4", "DP", "SPCD4", "SPCD8", "unknown", "Dendritic", "NK", "Bcell", "Monocyte"))

heatmap_endo <- ggplot(df_endo, aes(gating, orig.ident, fill = Endo1)) + geom_tile()

pdf('heatmap_endo.pdf')
heatmap_endo
dev.off()

saveRDS(heatmap_endo, "heatmap_endo.rds")
saveRDS(plot, "violin_endo.rds")


