library(seqinr)
library(VariantAnnotation)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
# my_chr <- c(1:22,'X','Y')
# my_chr <- gsub(pattern="^", replacement='chr', my_chr)
# 
# my_chr_size <- list()
# for (i in my_chr){
#   my_chr_size[[i]] <- length(BSgenome.Hsapiens.UCSC.hg38[[i]])
# }

# vcf_fname <- "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones/pt2229-DX1BM-TALL1.integrated.svs.filtered_geneAnnotation_intersect_uniq.vcf"
vcf_files1 <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/pt2322/","uniq.vcf",full.names = T, recursive = F)
vcf_files2 <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/pt2322/",".vcf.gz",full.names = T, recursive = F)
vcf_files3 <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/pt2229//","uniq.vcf",full.names = T, recursive = F)
vcf_files4 <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/pt2229//",".vcf.gz",full.names = T, recursive = F)
vcf_files5 <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/pt2283D//","uniq.vcf",full.names = T, recursive = F)
vcf_files6 <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/pt2283D/",".vcf.gz",full.names = T, recursive = F)

vcf_files <- c(vcf_files1, vcf_files2, vcf_files3, vcf_files4, vcf_files5, vcf_files6)
flanking_bases <- 50

get_context <- function(chr, pos ) {
  start = pos - flanking_bases
  end = pos + flanking_bases
  if ( chr %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) ) {
    context = Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), names = chr, start = start , end = end, as.character = T )
  } else {
    context = rep("N",100)
  }
  return(context)
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

create_random_df <- function(my_df) {
  my_df_random <- data.frame(size = my_df$size)
  my_df_random$chrom1 <- unlist(lapply(my_df_random$size,generateRandomChrom))
  my_df_random$pos1 <- unlist(mapply(generateRandomPos, my_df_random$chrom1, my_df_random$size))
  my_df_random$pos2 <- my_df_random$pos1+my_df_random$size
  my_df_random$chrom2 <- my_df_random$chrom1
  
  my_df_random[is.na(my_df_random$size),]$chrom2 <- unlist(lapply(my_df_random[is.na(my_df_random$size),]$size,generateRandomChrom))
  my_df_random[is.na(my_df_random$size),]$pos2 <- unlist(mapply(generateRandomPos, my_df_random[is.na(my_df_random$size),]$chrom2, my_df_random[is.na(my_df_random$size),]$size))
  
  my_df_random$context1 <- mapply(get_context, my_df_random$chrom1, my_df_random$pos1)
  my_df_random$context2 <- mapply(get_context, my_df_random$chrom2, my_df_random$pos2)
  
  return( my_df_random )
  
}

for ( vcf_fname in vcf_files ) {
  print( vcf_fname )
  vcf <- readVcf(vcf_fname)
  vcf <- vcf[which(unlist(grepl(":",alt(vcf)))),]

  my_df <- data.frame(
    chrom1=paste("chr",seqnames(vcf),sep=""),
    pos1=start(vcf),
    strand1=NA,
    chrom2 = unlist(gsub(".*(\\]|\\[)(.+):.+","chr\\2",alt(vcf))),
    pos2 = as.integer(unlist(gsub(".*(\\]|\\[).+:([0-9]+).+","\\2",alt(vcf)))),
    strand2=NA
  )
  my_df$size <- NA
  my_df[my_df$chrom1 == my_df$chrom2,]$size <- my_df[my_df$chrom1 == my_df$chrom2,]$pos2-my_df[my_df$chrom1 == my_df$chrom2,]$pos1
  my_df <- my_df[my_df$pos1 < my_df$pos2,]

  my_df$context1 <- mapply(get_context, my_df$chrom1, my_df$pos1)
  my_df$context2 <- mapply(get_context, my_df$chrom2, my_df$pos2)

  sequences = as.list(c(my_df$context1,my_df$context2))
  names(sequences) <- c(paste(paste(my_df$chrom1, my_df$pos1,sep=":"),paste(my_df$chrom2, my_df$pos2,sep=":"),"1",sep="_"),paste(paste(my_df$chrom1, my_df$pos1,sep=":"),paste(my_df$chrom2, my_df$pos2,sep=":"),"2",sep="_"))
  # outfile <- paste(dirname(vcf_fname),"/",basename(dirname(vcf_fname)),".fasta",sep="")
  outfile <- gsub(".vcf(.gz)*$",".fasta",vcf_fname)
  write.fasta(sequences, names(sequences),outfile, open = "w", nbchar = 60, as.string = FALSE)

  
  # 
  # for ( x in c(1:100) ) {
  #   my_df_random <- create_random_df( my_df )
  #   sequences = as.list(c(my_df_random$context1,my_df_random$context2))
  #   names(sequences) <- c(paste(paste(my_df_random$chrom1, my_df_random$pos1,sep=":"),paste(my_df_random$chrom2, my_df_random$pos2,sep=":"),"1",sep="_"),paste(paste(my_df_random$chrom1, my_df_random$pos1,sep=":"),paste(my_df_random$chrom2, my_df_random$pos2,sep=":"),"2",sep="_"))
  #   outfile <- paste(dirname(vcf_fname),"/",basename(dirname(vcf_fname)),"_random",x,".fasta",sep="")
  #   
  #   write.fasta(sequences, names(sequences),outfile, open = "w", nbchar = 60, as.string = FALSE)
  #   
  # }
}

