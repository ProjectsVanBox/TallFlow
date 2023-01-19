library(VariantAnnotation)
library(dplyr)
library(ape)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2229_new/withoutBulk/")

vcf_files <- list.files("branches_vcfs/",".vcf",full.names = T)

create_tree_table <- function( vcf_fname ) {
  myvcf <- readVcf(vcf_fname, "hg38")
  mysamples <- samples(header(myvcf))
  if ( grepl("root", vcf_fname) ) {
    tree_table <- matrix(1, nrow=nrow(myvcf), ncol=length(mysamples))
    colnames(tree_table) <- mysamples
  } else {
    tree_table <- matrix(0, nrow=nrow(myvcf), ncol=length(mysamples))
    colnames(tree_table) <- mysamples
    for (mys in strsplit(basename(vcf_fname),"_")[[1]]) {
      tree_table[,grepl(paste(mys,"$",sep=""), colnames(tree_table))] <- 1
    }
  }
  return( tree_table)
}

mytree_table_l <- lapply(vcf_files, create_tree_table)
mytree_table <- do.call(rbind, mytree_table_l)

co_occur_m = cbind(mytree_table, "root" = c(rep(0, nrow(mytree_table))))
#co_occur_m[rowSums(mytree_table) == ncol(mytree_table),"root"] <- 1
co_occur_m <- co_occur_m[,-which(c("pt2229-DX1BM-ALLBULK","pt2229-DX1BM-MSCBULK") %in% colnames(co_occur_m))]

tree = co_occur_m %>% t() %>% dist.gene() %>% nj() #neighbour joining tree construction
rooted_tree = root(tree, outgroup = "root", resolve.root = T) %>% drop.tip("root", trim.internal = T)

rooted_tree$root.edge <- length(which(rowSums(mytree_table) == ncol(mytree_table)))
tip.color = tip.color <- rep("black", length(rooted_tree$tip.label))

tip.color[grep(rooted_tree$tip.label, pattern = "TALL1$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL2$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL3$")] <- "#385725"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL4$")] <- "#385725"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL5$")] <- "#385725"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL6$")] <- "#E7298A"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL7$")] <- "#E7298A"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL8$")] <- "#B17BA6"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL9$")] <- "#B17BA6"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL10$")] <- "#B17BA6"

pdf("pt2229_tree.pdf")
plot(rooted_tree, use.edge.length = T, label.offset = 1, edge.color = "black", edge.width = 1, cex=1, show.tip.label = TRUE, 
     tip.color = tip.color, root.edge = T)

edgelabels(rooted_tree$edge.length, bg="white", col="black", font=2)
dev.off()
