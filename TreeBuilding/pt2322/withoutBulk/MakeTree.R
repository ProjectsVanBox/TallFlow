library(VariantAnnotation)
library(dplyr)
library(ape)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/withoutBulk/")

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
co_occur_m <- co_occur_m[,-which(c("pt2322-DX1BM-ALLBULK","pt2322-DX1BM-MONOBULK") %in% colnames(co_occur_m))]

tree = co_occur_m %>% t() %>% dist.gene() %>% nj() #neighbour joining tree construction
rooted_tree = root(tree, outgroup = "root", resolve.root = T) %>% drop.tip("root", trim.internal = T)

rooted_tree$root.edge <- length(which(rowSums(mytree_table) == ncol(mytree_table)))
tip.color = tip.color <- rep("black", length(rooted_tree$tip.label))

pal = c("#FFA500","#33A02C","#33A02C","#E7298A", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#1965B0")
names(pal) <- c("ALLBULK", "DN", "DNCD1aNeg","iSPCD4", "DP", "MSCBULK", "MONOBulk","DNCD1aPos", "SPCD4")

tip.color[grep(rooted_tree$tip.label, pattern = "TALL1$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL2$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL3$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL4$")] <- "#E7298A"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL5$")] <- "#E7298A"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL6$")] <- "#E7298A"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL7$")] <- "#B17BA6"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL8$")] <- "#B17BA6"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL9$")] <- "#B17BA6"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL10$")] <- "#1965B0"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL11$")] <- "#1965B0"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL12$")] <- "#1965B0"

pdf("pt2322_tree.pdf")
plot(rooted_tree, use.edge.length = T, label.offset = 1, edge.color = "black", edge.width = 1, cex=1, show.tip.label = TRUE, 
     tip.color = tip.color, root.edge = T)

edgelabels(rooted_tree$edge.length, bg="white", col="black", font=2)
dev.off()
