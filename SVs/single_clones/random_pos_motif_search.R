my_chr <- c(1:22,'X','Y')
my_chr <- gsub(pattern="^", replacement='chr', my_chr)

my_chr_size <- list()
for (i in my_chr){
  my_chr_size[[i]] <- length(BSgenome.Hsapiens.UCSC.hg38[[i]])
}

generateRandomPos <- function(n,chr,chr.sizes,width=1){
  random_chr <- sample(x=chr,size=n,prob=chr.sizes,replace=T)
  random_pos <- sapply(random_chr,function(chrTmp){sample(chr.sizes[chr==chrTmp],1)}) 
  res <- GRanges(random_chr,IRanges(random_pos,random_pos+width))
  return(res)
}

my_gr <- generateRandomPos(6000, my_chr, unlist(my_chr_size))

get_context <- function(gr) {
  if ( length(gr) == 0 ) {
    return( gr )
  }
  chr = seqnames(gr)
  start1 = as.integer( start(gr) - flanking_bases )
  end1 = start(gr)+flanking_bases

  context1 = Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), names = chr, start =start1 , end = end1 , as.character = T )

  names(context1) <- paste(chr,start1,end1,sep="_")

  gr$context1 = context1

  return(gr)
}

my_gr_context <- get_context(my_gr)

search_motif <- function( dnaseq, motif ) {
  matches <- matchPattern(as.character(motif[1]), dnaseq, max.mismatch = 1)
  match_coords <- paste(start(matches)-50,end(matches)-50,sep=":")
}

heptamer_motif <- DNAString("CACAGTG")
nonamer_motif <- DNAString("ACAAAAACC")
heptamer_motif_rc <- reverseComplement(heptamer_motif)
nonamer_motif_rc <- reverseComplement(nonamer_motif)


context1_heptamer_motif <- mapply(search_motif, gr_context$context1, as.character(heptamer_motif))
context1_heptamer_rc_motif <- mapply(search_motif, gr_context$context1, as.character(heptamer_motif_rc))

context1_nonamer_motif <- mapply(search_motif, gr_context$context1, as.character(nonamer_motif))
context1_nonmaer_rc_motif <- mapply(search_motif, gr_context$context1, as.character(nonamer_motif_rc))

my_gr_context$context1_heptamer_motif <- context1_heptamer_motif
my_gr_context$context1_heptamer_rc_motif <- context1_heptamer_rc_motif
my_gr_context$context1_nonamer_motif <- context1_nonamer_motif
my_gr_context$context1_nonmaer_rc_motif <- context1_nonmaer_rc_motif

write.table(as.data.frame(my_gr_context),"~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones/random_pos_RSSmotif.txt",quote=T,sep="\t")

