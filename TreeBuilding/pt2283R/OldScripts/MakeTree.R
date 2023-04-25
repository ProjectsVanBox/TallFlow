library(VariantAnnotation)
library(dplyr)
library(ape)

library(pheatmap)

### Read PTATO

source("~/hpc/pmc_vanboxtel/personal/smiddelkamp/tools/MutationalPatterns_Addon.R")
source("~/hpc/pmc_vanboxtel/personal/smiddelkamp/tools/PTATOplus.R")
#PTATO_SNV_dir <- "~/hpc/pmc_vanboxtel/projects/Fanconi/3_Output/PMCID561AAD/PTATO_old/snvs/PMCID561AAD/"
#PTATO_SNV_dir <- "~/hpc/pmc_vanboxtel/projects/Fanconi/3_Output/PMCFANC06/PTATO/snvs/PMCFANC06/"
PTATO_SNV_dir <- "~/hpc/pmc_vanboxtel/projects/MutDiff/3_Output/PMC21586/PTATO/snvs/PMC21586/"

Output_dir <- "~/hpc/pmc_vanboxtel/projects/PTA_manuscript/3_Output/PMC21586/"
Individual <- tail(unlist(strsplit(PTATO_SNV_dir, split = "/")),1 )

Inputs <- data.frame(VCF_files = list.files(path = PTATO_SNV_dir, pattern = ".ptato.CALLABLE.vcf.gz$", full.names = T, recursive = TRUE),
                     VCF_files_filtered = list.files(path = PTATO_SNV_dir, pattern = ".ptato.filtered.vcf.gz$", full.names = T, recursive = T)
)
Inputs$Sample <- gsub(pattern = ".ptato.CALLABLE.vcf.gz", replacement = "", x = gsub(".*_","",Inputs$VCF_files))
Inputs$ID <-  gsub(pattern = "PMC21586-BM", replacement = "", x = Inputs$Sample)
Inputs$ID <-  gsub(pattern = "PMC21586", replacement = "", x = Inputs$ID)

Sample_order <- c("HSC-PTAH2C4","MPP-PTAH2A2","CMP-PTAH2A3", "GMP-PTAH2A1", "PROMO-PTAD2D3",
                  "PROMY-PTAD2D4","MONO-PTAD2D5" ,  "NEUT-PTAD2D6", "AML-PTAD2B1","AML-PTAD2D1", "AMLBULK")
Inputs$ID <- factor(Inputs$ID, levels = Sample_order)
Inputs <- Inputs[order(Inputs$ID ),]

# Read the sample VCFs and collect the PTAprob values + VAFs (from the unfiltered VCF) and the PTAprobCutoff (from the header of thefiltered VCF)
grl_raw <- list()
#VCF_list < list()
for(Sample in Inputs$Sample){
  print(Sample)
  vcf_file_sample <- Inputs[Inputs$Sample == Sample, "VCF_files"]
  filteredvcf_file_sample <- Inputs[Inputs$Sample == Sample, "VCF_files_filtered"]
  grl_raw[[Sample]] <- readPTATOvcf(vcf_unfiltered = vcf_file_sample, vcf_filtered = filteredvcf_file_sample, VAF_threshold = 0.15)
  #grl_raw[[Sample]]$Individual <- individual
  grl_raw[[Sample]]$ID <- Inputs[Inputs$Sample == Sample, "ID"]
}

grl <- lapply(grl_raw, function(x) x[which(x$FILTER =="PASS"),])
names(grl) = Inputs$ID
grl = GRangesList(grl)
seqlevels(grl, pruning.mode = "fine") = paste0("chr", 1:22)
#saveRDS(grl, file.path(Output_dir, "grl.rds"))
snv_grl <-  get_mut_type(grl, type = "snv")

### Tree
# Read the multi-sample VCF for the tree (doesnt contain the PTAprobs yet)
vcf_fname <- "~/hpc/pmc_vanboxtel/projects/MutDiff/3_Output/PMC21586/Tree/PMC21586_merged.vcf.gz"
vcf <- readVcf(vcf_fname, "hg38")
vcf = vcf[seqnames(vcf) %in% c(1:22),]
vcf_snv <- vcf[isSNV(vcf)]

samples <-  c("PMC21586AMLBULK" ,"PMC21586-BMCMP-PTAH2A3","PMC21586-BMMONO-PTAD2D5","PMC21586-BMPROMO-PTAD2D3","PMC21586-BMAML-PTAD2B1","PMC21586-BMGMP-PTAH2A1","PMC21586-BMMPP-PTAH2A2","PMC21586-BMPROMY-PTAD2D4","PMC21586-BMAML-PTAD2D1","PMC21586-BMHSC-PTAH2C4","PMC21586-BMNEUT-PTAD2D6")

tree_table <- matrix(nrow=length(rowRanges(vcf_snv)),ncol=length(samples))
colnames(tree_table) <- samples

for ( sample in samples) {
  tmp_table <- unlist(lapply(info(vcf_snv)$CLONAL_SAMPLE_NAMES, function(x) match(sample, x)))
  tmp_table2 <- unlist(lapply(info(vcf_snv)$FAIL_QC_SAMPLE_NAMES, function(x) match(sample, x)))
  
  tmp_table[is.na(tmp_table)] <- 0
  tmp_table[tmp_table > 0] <- 1
  # Sample that FAIL_QC get a 1 in tmp_table2. So everything that's is a 1 (and not na) should be removed
  tmp_table[!is.na(tmp_table2)] <- NA
  tree_table[,sample] <- tmp_table
}

tree_table <- as.data.frame(tree_table)
rownames(tree_table) <- rownames(vcf_snv)


# Select only the variants taht are PASS in at least one sample (eg remove the false pos that are unique for a sample)
PassVariantsPTA <- c()

for(sample in names(snv_grl)){
  print(sample)
  PassVariantsPTA <- c(PassVariantsPTA, names(snv_grl[[sample]]))
}

PassVariantsPTA <- PassVariantsPTA[!duplicated(PassVariantsPTA)]





#tree_table = tree_table[rownames(tree_table) %in% PassVariants,]


# Make tree using gt instead of clonal variants.
get_gt = function(vcf){
  gt = geno(vcf)$GT
  #gt[gt == "1/1" | gt == "0/1"] = 1 #Makes sure homozygous mutations on the x chromosome are placed in the same combi as muts on the autosomes.
  gt[gt == "0/0" | gt == "./."] = 0
  gt[gt != 0] = 1
  
  #Transform to integer. Also works for single row gt.
  col_names = colnames(gt)
  row_names = rownames(gt)
  gt = purrr::map(seq_len(nrow(gt)), function(i) as.integer(gt[i,])) %>% 
    do.call(rbind, .)
  colnames(gt) = col_names
  rownames(gt) = row_names
  return(gt)
}
gt = get_gt(vcf_snv)
#gt = gt[,colnames(gt) %in% samples]

pheatmap(gt, show_rownames = F)

gt <- gt[,-which(colnames(gt) %in% c("PMC21586AMLBULK", "PMC21586MSCBULK", "2:PMC21586MSCBULK", "2:PMC21586AMLBULK"))]
gt <- gt[rowSums(gt) > 0,]

gt_PTA <- gt[,which(colnames(gt) %in% samples)]
gt_PTA[rowSums(gt_PTA) == 8,]
gt_clone <- gt[,-which(colnames(gt) %in% samples)]
gt_clone <- gt_clone[rowSums(gt_clone) > 0,]
PassVariantsClones <- rownames(gt_clone)
PassVariants <- c(PassVariantsPTA, PassVariantsClones)
PassVariants <- PassVariants[!duplicated(PassVariants)]
gt = gt[rownames(gt) %in% PassVariants,]

# shared_vcf = vcf_snv[rownames(gt) %in% PassVariants & rowSums(gt) >= 2,]
# writeVcf(shared_vcf, file.path(Output_dir, "PMC21586_shared_FM.vcf"))
library(stringr)
#gt = gt[rownames(gt) %in% PassVariants,]
co_occur_m = cbind(gt, "root" = c(rep(0, nrow(gt))))
colnames(co_occur_m) = co_occur_m %>% 
  colnames() %>% 
  str_remove("PMC21586(-BM)*")

tree = co_occur_m %>% t() %>% dist.gene() %>% nj() #neighbour joining tree construction
rooted_tree = root(tree, outgroup = "root", resolve.root = T) %>% drop.tip("root", trim.internal = T)
tree_edge_cols = rep("black", nrow(rooted_tree$edge))

tip.color <- rep("wheat4", length(rooted_tree$tip.label))

tip.color[grep(rooted_tree$tip.label, pattern = "PTA")] <- "steelblue3"

#saveRDS(rooted_tree, file.path(Output_dir,"rooted_tree_FM.rds"))
pdf("~/Downloads/Test_tree.pdf", width = 4, height = 3, pointsize = 6)
par(oma = rep(0, 4))
my_tree_plot <- plot(rooted_tree, use.edge.length = T, label.offset = 1, edge.color = "black", edge.width = 1, cex=1, show.tip.label = TRUE, 
                     tip.color = tip.color)
axis(1,at =seq(from =0, to = 1000, by = 50),cex.axis = 1, labels = FALSE)
axis(1,at =seq(from =0, to = 1000, by = 100),cex.axis = 1)
dev.off()
plotBreakLongEdges(tree)
axis(1,at =seq(from =0, to = 1000, by = 50),cex.axis = 1, labels = FALSE)
axis(1,at =seq(from =0, to = 1000, by = 100),cex.axis = 1)
edgelabels()
dev.off()

# Check
check <- gt[,c("PMC21586MPP6", "PMC21586-BMMPP-PTAH2A2", "PMC21586MPP4")]
check[rowSums(check) > 1,]


pdf(paste(Output_dir, "PMC21586_Heatmap_sharedSNV_FM.pdf", sep = ""), width = 8, height=10)
pheatmap(gt[rowSums(gt) > 1,], show_rownames = F, na_col = "grey")
dev.off()