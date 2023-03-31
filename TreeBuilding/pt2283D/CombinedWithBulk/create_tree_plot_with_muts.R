# Tree - script 
# Axel Rosendahl Huber & Mark van Roosmalen
# Script to generate mutations 

# ChemoEffects lineage analysis 
library("VariantAnnotation")
library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)
library(dplyr)
library(ape)
library(UpSetR)
library(ggtree)
library(dendextend)
library(phangorn)
library(ggplot2)
library(microplot)
library(ggpubr)
library(tidyr)
library(gridExtra)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/")

vcf = "./branches_vcfs/branches_merged.vcf"
#~/surfdrive/Shared/Projects/2ndCancer_Axel_Eline_Jurrian/Mutational_profile_analysis//Data/Tree_vcfs/11266.filtered_variants_snpEff_snpSift_Cosmicv89_SMuRF_tree_all.vcf"
sample_name  = gsub("\\..*$", "",basename(vcf_file))

COLORS6 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE")

# TODO make function for all samples (including sub-functions)

inset_horizontal_bar = function (tree_view, insets, width, height, hjust = 0, vjust = 0, x = "node", reverse_x = FALSE, reverse_y = FALSE) {
  #  if (width < 0 || width > 1) 
  #    stop("width should be in range of (0,1)")
  #  if (height < 0 || height > 1) 
  #    stop("height should be in range of (0,1)")
  df <- tree_view$data[as.numeric(names(insets)), ]
  x <- match.arg(x, c("node", "branch", "edge"))
  if (x == "node") {
    xx <- df$x
  }
  else {
    xx <- df$branch
  }
  yy <- df$y
  xx <- xx - hjust
  yy <- yy - vjust
  if (reverse_x) 
    xx <- -xx
  if (reverse_y) 
    yy <- -yy
  #width <- width * diff(range(tree_view$data$x, na.rm = TRUE))
  height <- height * diff(range(tree_view$data$y, na.rm = TRUE))
  subview_barplot =   tree_view + ggimage::geom_subview(subview = insets, width = width, 
                                                        height = height, x = xx, y = yy) + coord_flip()
  return(subview_barplot)
}

create_tree_table <- function( tree_vcf_name ) {
  
  sample_name  = gsub("\\..*$", "",basename(tree_vcf_name))
  vcf.raw <- readVcf(tree_vcf_name, "hg38")
  vcf <- vcf.raw
  mult_alts = elementNROWS(rowRanges(vcf)$ALT) > 1
  nr_mult_alts = sum(mult_alts)
  
  # Remove varia with multiple alternative alleles
  if (nr_mult_alts > 0){
    rowRanges(vcf) = rowRanges(vcf)[!mult_alts]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  
  
  gr = granges(vcf)

  
  tree_table <- matrix(nrow=length(rowRanges(vcf)),ncol=length(samples(header(vcf))))
  mut = paste0(gr$REF,".", as.character(unlist(gr$ALT)))
  
  mut <- recode(mut, A.C = "T.G", A.G = "T.C", A.T = "T.A", G.A = "C.T", G.T = "C.A", G.C = "C.G")
  index_dels = which(nchar(gsub("\\..*$", "", mut)) > nchar(gsub("^.*\\.", "", mut)))
  index_insertions = which((nchar(gsub("\\..*$", "", mut)) <= nchar(gsub("^.*\\.", "", mut))) & nchar(gsub("^.*\\.", "", mut)) > 1)
  mut[index_dels] = "Del"
  mut[index_insertions] = "Ins"
  # mut = mut[-c(index_dels, index_insertions)]
  
  rownames(tree_table) = paste0(seqnames(gr), "_", start(gr), "_", mut)
  colnames(tree_table) <- samples(header(vcf))
  
  for ( sample in samples(header(vcf))) {
    tmp_table <- unlist(lapply(info(vcf)$CLONAL_SAMPLE_NAMES, function(x) match(sample, x)))
    tmp_table[is.na(tmp_table)] <- 0
    tmp_table[tmp_table > 0] <- 1
    tree_table[,sample] <- tmp_table
  }
  
  tree_table <- as.data.frame(tree_table)
  
  # Deselect certain samples (shared with AML or T_ALL samples)
  # tree_table = tree_table[colnames(tree_table) != "PMC10925DX1HSCN13"]
  
  # Rename colnames
  # colnames(tree_table)[colnames(tree_table) == "PMC10925DX1HSCN13"] = "PMC10925DX1HSCN13"
  
  #Exclude MSC
  tree_table <- tree_table[,-grep("BULK",colnames(tree_table))]
  
  # SHARED MUTATIONS
  tree_table_out = list()
  tree_table_out$shared <- tree_table[rowSums(tree_table) > 1 & rowSums(tree_table) < ncol(tree_table)-1,]
  
  #UNIQUE MUTATIONS
  tree_table_out$unique <- tree_table[rowSums(tree_table) == 1,]
  
  return( tree_table_out )
}

create_shared_tree <- function( co_occur_m ) {
  clone_names = colnames(co_occur_m)
  co_occur_m = cbind(co_occur_m, "root" = c(rep(0, nrow(co_occur_m))))
  
  co_occur_m_adjusted =   mutate_all(co_occur_m, as.character)
  for (i in 1:nrow(co_occur_m)) {
    row = co_occur_m[i,]
    mut = gsub(".*_(.{3})\\.*\\d*", "\\1", rownames(row))
    row[row == "1"] = mut
    co_occur_m_adjusted[i,] = row
  }
  
  
  co_occur_m_adjusted$index_co_occurrence = "0"
  for (i in 1:nrow(co_occur_m_adjusted)) {
    mutation_row = co_occur_m_adjusted[i,]
    co_occur_m_adjusted$index_co_occurrence[i] =  paste(colnames(mutation_row)[mutation_row != "0"],collapse = ",")
  }
  
  table(co_occur_m_adjusted$index_co_occurrence)
  
  ocurrence_matrices = split(co_occur_m_adjusted, co_occur_m_adjusted$index_co_occurrence)
  
  if (nrow(co_occur_m_adjusted) == 1) {
    co_occur_m_adjusted <- rbind(co_occur_m_adjusted,0)
  }
  tree = co_occur_m_adjusted[, names(co_occur_m_adjusted) != "index_co_occurrence"] %>% t() %>% dist.gene(method="pairwise") %>% nj() 
  rooted_tree = root(tree, outgroup = "root", resolve.root = T) %>% drop.tip("root", trim.internal = T)
  
  plot = ggtree(rooted_tree)
  
  # Generate sample_type_occurrences---------------------------.
  type_ocurrences = data.frame(C.A = rep(0, length(ocurrence_matrices)), C.G = 0,  C.T = 0,  Del = 0,  Ins = 0,   T.A = 0,  T.C = 0,  T.G= 0)
  rownames(type_ocurrences) = names(ocurrence_matrices)
  for (i in names(ocurrence_matrices)) {
    if (i == "") {next}
    mat = ocurrence_matrices[[i]]
    mat = mat[,names(mat) != "index_co_ocurrence"]
    
    mut_ocurrences = type_ocurrences[i,]
    select_columns = colnames(mat)[apply(mat,2,  unique) != "0"]
    select_column = select_columns[1]
    table_mut_ocurrence = table(mat[,select_column])
    mut_ocurrences[,names(table_mut_ocurrence)] =  table_mut_ocurrence
    type_ocurrences[i,] = mut_ocurrences
  }
  
  internals = plot$data[plot$data$isTip == F,]
  descendants = phangorn::Descendants(rooted_tree)
  type_ocurrences$node = match(row.names(type_ocurrences), plot$data$label[plot$data$branch.length > 0], nomatch = NA)
  
  
  children_nodes = descendants[internals$node[internals$branch.length > 0]]
  children_names = lapply(children_nodes, function(x) plot$data$label[x])
  children_names = lapply(children_names, function(x) paste(x,collapse = ","))
  
  type_ocurrences$node[match(children_names, row.names(type_ocurrences))] = internals$node[internals$branch.length > 0]
  bars <- nodebar(type_ocurrences[,c(5,4,8:6,3:1,9)], cols=8:1, color = rev(c(COLORS6, "#fd8002","#35a02e")), position = "stack")
  bars = lapply(bars, function(x) x + scale_y_continuous(expand = c(0, 0)) + labs(x=NULL, y=NULL) + theme(plot.margin=unit(c(0,0,0,0),"mm"), axis.text = element_blank(), axis.ticks.length = unit(0, "mm")))
  
  p_plus = plot
  bar_height <- 0.36/(ncol(co_occur_m)-1)
  
  # --- Add barplots to dendrogram
  for (i in names(bars)) {
    j = match(i, type_ocurrences$node)
    i = p_plus = inset_horizontal_bar(tree_view = p_plus, bars[i], width=sum(type_ocurrences[j,1:8]), height=bar_height, x = "edge", hjust=0, reverse_x=TRUE)
  }
  
  # Futher stilistically append plot 
  plot = p_plus + theme(axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), 
                        axis.ticks.y = element_line(colour = 'black', size=0.5, linetype='solid'), 
                        axis.text.y = element_text() ) +
        theme(plot.margin = unit(c(1,0,0,0), "cm")) +
        xlab("Number of shared mutations")

  return(plot)
}

create_unique_tree <- function( co_occur_m, sample_order ) {
  clone_names = colnames(co_occur_m)
  co_occur_m = cbind(co_occur_m, "root" = c(rep(0, nrow(co_occur_m))))
  
  co_occur_m_adjusted =   mutate_all(co_occur_m, as.character)
  for (i in 1:nrow(co_occur_m)) {
    row = co_occur_m[i,]
    mut = gsub(".*_(.{3})\\.*\\d*", "\\1", rownames(row))
    row[row == "1"] = mut
    co_occur_m_adjusted[i,] = row
  }
  
  
  co_occur_m_adjusted$index_co_occurrence = "0"
  for (i in 1:nrow(co_occur_m_adjusted)) {
    mutation_row = co_occur_m_adjusted[i,]
    co_occur_m_adjusted$index_co_occurrence[i] =  paste(colnames(mutation_row)[mutation_row != "0"],collapse = ",")
  }
  bars <- co_occur_m_adjusted %>% as_tibble() %>% 
    pivot_longer(cols = colnames(co_occur_m_adjusted)[1:(ncol(co_occur_m_adjusted)-2)]) %>%
    group_by(name, value) %>% 
    count %>% 
    filter( value != '0')
  
  bars$label <- bars$name
  bars$name <- factor(bars$name, levels=sample_order)
  bars <- as.data.frame(bars)
  bars$value <- factor(bars$value, levels= c("C.A","C.G","C.T","T.A","T.C","T.G","Ins","Del"))
  bars$n <- as.numeric(bars$n)
  totals <- as.vector(by(bars$n, bars$label, sum))
  labels <- data.frame(totals=totals, label=unique(bars$label), name=factor(unique(bars$label), levels=sample_order))

  plot <- ggplot( bars, aes(x=name, y=n)) +
    geom_bar(aes(fill=value), stat="identity", width=0.3) +
    theme_classic() +
    scale_fill_manual(values=c(COLORS6,"#fd8002","#35a02e")) +
    geom_text(data=labels,aes(y=totals,label=label), fill=NA, vjust=0, hjust="right", angle=90) +
    ylab("Number of unique mutations")
  
  return(plot)
}

resize_heights <- function(g, heights = rep(1, length(idpanels))){
  idpanels <- unique(g$layout[grepl("panel",g$layout$name), "t"])
  g$heights <- grid:::unit(g$heights)
  g$heights[idpanels] <- unit(do.call(unit, list(heights, 'null')))
  g
}

align_plots1 <- function (p1, p2, plot_heights) {
  pl <- list(p1, p2)
  stopifnot(do.call(all, lapply(pl, inherits, "gg")))
  gl <- lapply(pl, ggplotGrob)
  bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
  combined <- Reduce(bind2, gl[-1], gl[[1]])
  wl <- lapply(gl, "[[", "widths")
  hl <- lapply(gl, "[[", "heights")
  combined$widths <- do.call(grid::unit.pmax, wl)
  grid::grid.newpage()
  grid::grid.draw(resize_heights(combined, plot_heights))
}

tree_table <- create_tree_table(vcf)
shared_max_y <- max(colSums(tree_table$shared))

shared_plot <- create_shared_tree( tree_table$shared ) +
scale_x_reverse(expand = c(0,0), limits=c(shared_max_y,-0.2), breaks=seq(shared_max_y,0,by=-1))

sample_order <- shared_plot$data[order(shared_plot$data[shared_plot$data$isTip,]$y),]$label
segments <- data.frame(sample=sample_order)
rownames(segments) <- segments$sample
segments$xend <- shared_max_y
segments$y <- c(1:length(segments$sample))
segments$yend <- segments$y
segments <- merge(segments,t(t(colSums(tree_table$shared))),by="row.names")
colnames(segments)[ncol(segments)] <- "x"

shared_plot <- shared_plot +
  geom_segment(data=segments, aes(x = x, y = y, xend = xend, yend = yend), linetype="dashed")
  
unique_max_y <- max(colSums(tree_table$unique))

unique_plot <- create_unique_tree( tree_table$unique, sample_order )
unique_plot <- unique_plot + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.line.x = element_blank()) +
  scale_y_reverse(expand = c(0,0), limits=c(unique_max_y*1.5,0)) + theme(plot.margin = unit(c(0,0,0,0), "cm") )

height_ratio <- shared_max_y*10/unique_max_y
if (height_ratio < 0.3 ) { height_ratio = 0.3 }
plot_heights <- c(height_ratio, 1)

pdf(paste(sample_name,".pdf",sep=""))
print( align_plots1(shared_plot, unique_plot, plot_heights) )
dev.off()

