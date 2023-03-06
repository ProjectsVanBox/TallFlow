library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/IndelSignatures/pt2322")

vcf_files <- list.files("../../TreeBuilding_indel/pt2322/branches_vcfs/","branch.vcf",full.names = T)
sample_names <- gsub(".+/(.+)_branch.vcf","\\1",vcf_files)

filter_gt_bulk <- function( myvcf_file ) {
  myvcf <- readVcf(myvcf_file)
  myrows <- as.numeric(which(geno(myvcf)$GT[,1] != "0/0"))
  return( myrows )
}

myrows <- lapply(vcf_files, filter_gt_bulk)
names(myrows) <- sample_names


indel_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "indel")

indel_grl <- get_indel_context(indel_grl, ref_genome)
for (n in names(indel_grl)) {
  indel_grl[[n]] <- indel_grl[[n]][myrows[[n]],]
}

indel_counts <- count_indel_contexts(indel_grl)

pooled_indel_counts <- pool_mut_mat(indel_counts, grouping = c(rep("rest",8),"root"))

pdf("pt2322_signatures_overview.pdf")

plot_indel_contexts(pooled_indel_counts, condensed = TRUE)

indel_signatures = get_known_signatures(muttype="indel")
fit_res <- fit_to_signatures(pooled_indel_counts, indel_signatures)

plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)

strict_refit <- fit_to_signatures_strict(pooled_indel_counts, indel_signatures, max_delta = 0.07)
fit_res_strict <- strict_refit$fit_res
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)
dev.off()