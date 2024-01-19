library(VariantAnnotation)
library(dplyr)
library(ape)

setwd("~/hpc/projects/TallFlow/3_Output/TreeBuilding/pt1909/")

vcf_files <- list.files("branches_vcfs",".vcf",full.names = T)
#indel_vcf_files <- list.files("~/hpc/projects/TallFlow/3_Output/TreeBuilding_indel/pt2229/branches_vcfs/","branch.vcf",full.names = T)

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
#head(mytree_table_l)
#indel_mytree_table_l <- lapply(indel_vcf_files, create_tree_table)
#indel_mytree_table <- do.call(rbind, indel_mytree_table_l)

#mytree_table <- rbind(mytree_table, indel_mytree_table)

co_occur_m = cbind(mytree_table, "root" = c(rep(0, nrow(mytree_table))))
#co_occur_m[rowSums(mytree_table) == ncol(mytree_table),"root"] <- 1
co_occur_m <- co_occur_m[,!colnames(co_occor_m) %in% c("pt1909-DX1BM-ALLBULK","pt1180-DX1BM-MSCBULK")]

tree = co_occur_m %>% t() %>% dist.gene() %>% nj() #neighbour joining tree construction
rooted_tree = root(tree, outgroup = "root", resolve.root = T) %>% drop.tip("root", trim.internal = T)

rooted_tree$root.edge <- length(which(rowSums(mytree_table) == ncol(mytree_table)))
tip.color = tip.color <- rep("black", length(rooted_tree$tip.label))



# colors for plotting
#pal = c("#FFA500","#33A02C","#33A02C", "#1965B0", "#B17BA6", "#A4A4A4",  "#A4A4A4", "#385725", "#E7298A", "#E28028", "#963C25", "grey")
#names(pal) <- c("ALLBULK", "DN-like", "DNearly","SPCD4", "DP", "MSCBULK", "MONOBulk","DN3", "iSPCD4", "SPCD8", "gdTcell", "unknown")
#pal

tip.color[grep(rooted_tree$tip.label, pattern = "TALL1$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL2$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL3$")] <- "#33A02C"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL4$")] <- "#385725"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL5$")] <- "#385725"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL6$")] <- "#385725"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL7$")] <- "#1965B0"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL8$")] <- "#1965B0"
tip.color[grep(rooted_tree$tip.label, pattern = "TALL9$")] <- "#1965B0"


pdf("/Users/verapoort/Desktop/pt1909_tree.pdf")
plot(rooted_tree, use.edge.length = T, label.offset = 1, edge.color = "black", edge.width = 1, cex=1, show.tip.label = TRUE, 
     tip.color = tip.color, root.edge = T)

edgelabels(rooted_tree$edge.length, bg="white", col="black", font=2)
dev.off()
