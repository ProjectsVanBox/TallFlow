library(Biostrings)
library(VariantAnnotation)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
flanking_bases <- 50

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


search_motif <- function( dnaseq, motif ) {
  matches <- matchPattern(as.character(motif[1]), dnaseq, max.mismatch = 1)
  match_coords <- paste(start(matches)-50,end(matches)-50,sep=":")
}


vcf_fnames <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones/","*RSSmotif.vcf$",full.names = T)
for (vcf_fname in vcf_fnames) {
  vcf <- readVcf(vcf_fname)
  hm1 <- which(unlist(lapply(info(vcf)$HM1,length)) > 0)
  hmrc1 <- which(unlist(lapply(info(vcf)$HMRC1,length)) > 0)
  nm1 <- which(unlist(lapply(info(vcf)$NM1,length)) > 0)
  nmrc1 <- which(unlist(lapply(info(vcf)$NMRC1,length)) > 0)
  
  my_random_matrix <- matrix(ncol=4)[-1,]
  for ( n in c(1:100)) {
    my_random_gr <- generateRandomPos(length(vcf), my_chr, unlist(my_chr_size))
    my_random_gr_context <- get_context(my_random_gr)
    
    heptamer_motif <- DNAString("CACAGTG")
    nonamer_motif <- DNAString("ACAAAAACC")
    heptamer_motif_rc <- reverseComplement(heptamer_motif)
    nonamer_motif_rc <- reverseComplement(nonamer_motif)
    
    
    context1_heptamer_motif <- mapply(search_motif, my_random_gr_context$context1, as.character(heptamer_motif))
    context1_heptamer_rc_motif <- mapply(search_motif, my_random_gr_context$context1, as.character(heptamer_motif_rc))
    
    context1_nonamer_motif <- mapply(search_motif, my_random_gr_context$context1, as.character(nonamer_motif))
    context1_nonmaer_rc_motif <- mapply(search_motif, my_random_gr_context$context1, as.character(nonamer_motif_rc))
    
    hm1_random <- which(unlist(lapply(context1_heptamer_motif,length))>0)
    hmrc1_random <- which(unlist(lapply(context1_heptamer_rc_motif,length))>0)
    nm1_random <- which(unlist(lapply(context1_nonamer_motif,length))>0)
    nmrc1_random <- which(unlist(lapply(context1_nonmaer_rc_motif,length))>0)
    
    my_random_matrix <- rbind(my_random_matrix, c(length(hm1_random),length(hmrc1_random), length(nm1_random), length(nmrc1_random)))
  }
  
  my_random_df <- as.data.frame(my_random_matrix)
  plot(density(my_random_matrix[,1]))
  abline(v=length(hm1))
  
}

