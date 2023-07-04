library(Biostrings)
library(VariantAnnotation)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
flanking_bases <- 50

my_chr <- c(1:22,'X','Y')
my_chr <- gsub(pattern="^", replacement='chr', my_chr)

my_chr_size <- list()
for (i in my_chr){
  my_chr_size[[i]] <- length(BSgenome.Hsapiens.UCSC.hg38[[i]])
}

generateRandomChrom <- function( sv_size ) {
  if (is.na(sv_size)) {
    chroms <- my_chr
  } else {
    chroms <- names(my_chr_size[my_chr_size > sv_size])
  }
  chr.sizes <- my_chr_size[names(my_chr_size) %in% chroms]
  random_chr <- sample(x=chroms,size=1,prob=chr.sizes,replace=T)
  return( random_chr )
}

generateRandomPos <- function(chr,sv_size){
  chr_size <- as.numeric(my_chr_size[[chr]])-flanking_bases
  if (!is.na(sv_size)) {
    chr_size <- chr_size-sv_size
  }
  random_pos <- sample(chr_size,1)
  return(random_pos)
}

get_context2 <- function(chr, pos, strand, context) {
  if (strand == "+" & context == 1) {
    start = pos - flanking_bases
    end = pos
  } else if ( strand == "+" & context == 2 ) {
    start = pos
    end = pos + flanking_bases
    strand = "-"
  } else if ( strand == "-" & context == 1 ) {
    start = pos
    end = pos + flanking_bases
    strand = "+"
  } else if (strand == "-" & context == 2 ) {
    start = pos - flanking_bases
    end = pos
  }
  context = Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), names = chr, start = start , end = end, strand = strand, as.character = T )
  return(context)
}

get_context <- function(chr, pos ) {
  start = pos - flanking_bases
  end = pos + flanking_bases
  context = Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), names = chr, start = start , end = end, as.character = T )
  return(context)
}

search_motif <- function( dnaseq, motif ) {
  matches <- matchPattern(DNAString(motif), DNAString(dnaseq), max.mismatch = 1)
  motif_found <- FALSE
  if ( length(matches) > 0 ) {
    motif_found <- TRUE
  }
  return( motif_found )
}

search_motif2 <- function( dnaseq, motif ) {
  matches <- matchPattern(DNAString(motif), DNAString(dnaseq), max.mismatch = 1, fixed = F)
  motif_found <- FALSE
  if ( length(matches) > 0 ) {
    motif_found <- TRUE
  }
  return( motif_found )
}

vcf_fname <- "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones/pt2229-DX1BM-TALL1.integrated.svs.filtered_geneAnnotation_intersect_uniq.vcf"
vcf <- readVcf(vcf_fname)

my_df <- data.frame(
          chrom1=paste("chr",seqnames(vcf),sep=""),
          pos1=start(vcf),
          strand1=NA,
          chrom2 = unlist(gsub(".*(\\]|\\[)(.+):.+","chr\\2",alt(vcf))),
          pos2 = as.integer(unlist(gsub(".*(\\]|\\[).+:([0-9]+).+","\\2",alt(vcf)))),
          strand2=NA
          )

my_df[which(unlist(grepl("^\\w+\\[",alt(vcf)))),]$strand1 <- "+"
my_df[which(unlist(grepl("^\\w+\\[",alt(vcf)))),]$strand2 <- "+"

my_df[which(unlist(grepl("^\\w+\\]",alt(vcf)))),]$strand1 <- "+"
my_df[which(unlist(grepl("^\\w+\\]",alt(vcf)))),]$strand2 <- "-"

my_df[which(unlist(grepl("^\\]",alt(vcf)))),]$strand1 <- "-"
my_df[which(unlist(grepl("^\\]",alt(vcf)))),]$strand2 <- "-"

my_df[which(unlist(grepl("^\\[",alt(vcf)))),]$strand1 <- "-"
my_df[which(unlist(grepl("^\\[",alt(vcf)))),]$strand2 <- "+"

my_df$size <- NA
my_df[my_df$chrom1 == my_df$chrom2,]$size <- my_df[my_df$chrom1 == my_df$chrom2,]$pos2-my_df[my_df$chrom1 == my_df$chrom2,]$pos1
my_df$context1 <- mapply(get_context, my_df$chrom1, my_df$pos1)
my_df$context2 <- mapply(get_context, my_df$chrom2, my_df$pos2)
my_df <- my_df[my_df$pos1 < my_df$pos2,]
 
my_df_random <- data.frame(size = my_df$size)
my_df_random$chrom1 <- unlist(lapply(my_df_random$size,generateRandomChrom))
my_df_random$pos1 <- unlist(mapply(generateRandomPos, my_df_random$chrom1, my_df_random$size))
my_df_random$strand1 <- my_df$strand1
my_df_random$strand2 <- my_df$strand2
my_df_random$pos2 <- my_df_random$pos1+my_df_random$size
my_df_random$chrom2 <- my_df_random$chrom1

my_df_random[is.na(my_df_random$size),]$chrom2 <- unlist(lapply(my_df_random[is.na(my_df_random$size),]$size,generateRandomChrom))
my_df_random[is.na(my_df_random$size),]$pos2 <- unlist(mapply(generateRandomPos, my_df_random[is.na(my_df_random$size),]$chrom2, my_df_random[is.na(my_df_random$size),]$size))

my_df_random$context1 <- mapply(get_context, my_df_random$chrom1, my_df_random$pos1)
my_df_random$context2 <- mapply(get_context, my_df_random$chrom2, my_df_random$pos2)

heptamer_motif <- "CACAGTG"
nonamer_motif <- "ACAAAAACC"
heptamer_motif2 <- "CACNNNN"
nonamer_motif2 <- "NNNNAANN"

rss12_motif <- paste(c("CACNNNN",rep("N",12),"NNNNAANN"),collapse = "")
rss23_motif <- paste(c("CACNNNN",rep("N",23),"NNNNAANN"),collapse = "")


my_df$heptamer1 <- unlist(mapply(search_motif, my_df$context1, heptamer_motif2))
my_df$nonamer1 <- unlist(mapply(search_motif, my_df$context1, nonamer_motif))
my_df$heptamer2 <- unlist(mapply(search_motif, my_df$context2, heptamer_motif))
my_df$nonamer2 <- unlist(mapply(search_motif, my_df$context2, nonamer_motif))

my_df$rss12_1 <- unlist(mapply(search_motif2, my_df$context1, rss12_motif))
my_df$rss12_2 <- unlist(mapply(search_motif2, my_df$context2, rss12_motif))
my_df$rss23_1 <- unlist(mapply(search_motif2, my_df$context1, rss23_motif))
my_df$rss23_2 <- unlist(mapply(search_motif2, my_df$context2, rss23_motif))

my_df_random$heptamer1 <- unlist(mapply(search_motif, my_df_random$context1, heptamer_motif))
my_df_random$nonamer1 <- unlist(mapply(search_motif, my_df_random$context1, nonamer_motif))
my_df_random$heptamer2 <- unlist(mapply(search_motif, my_df_random$context2, heptamer_motif))
my_df_random$nonamer2 <- unlist(mapply(search_motif, my_df_random$context2, nonamer_motif))

my_df_random$rss12_1 <- unlist(mapply(search_motif2, my_df_random$context1, rss12_motif))
my_df_random$rss12_2 <- unlist(mapply(search_motif2, my_df_random$context2, rss12_motif))
my_df_random$rss23_1 <- unlist(mapply(search_motif2, my_df_random$context1, rss23_motif))
my_df_random$rss23_2 <- unlist(mapply(search_motif2, my_df_random$context2, rss23_motif))

for (dnaseq in my_df_random$context1) {
  print(dnaseq)
  matches <- matchPattern(DNAString(motif), DNAString(dnaseq), max.mismatch = 1)
  motif_found <- FALSE
  if ( length(matches) > 0 ) {
    motif_found <- TRUE
  }
}


#my_df[(my_df$heptamer1 == T & my_df$nonamer1 == T) | (my_df$heptamer2 == T & my_df$nonamer2 == T),]
#my_df_random[(my_df_random$heptamer1 == T & my_df_random$nonamer1 == T) | (my_df_random$heptamer2 == T & my_df_random$nonamer2 == T),]
  
