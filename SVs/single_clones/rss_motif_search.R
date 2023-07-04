library(Biostrings)
library(VariantAnnotation)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
flanking_bases <- 50

vcf_fnames <- list.files("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/SVs/single_clones/","*uniq.vcf$",full.names = T)
for (vcf_fname in vcf_fnames) {
  vcf <- readVcf(vcf_fname)
  seqlevels(vcf) <- paste0("chr", seqlevels(vcf))
  gr = granges(vcf)
  
  get_context <- function(gr) {
    if ( length(gr) == 0 ) {
      return( gr )
    }
    chr = seqnames(gr)
    start1 = as.integer( start(gr) - flanking_bases )
    end1 = start(gr)+flanking_bases
    chr2 = unlist(gsub(".*(\\]|\\[)(.+):.+","chr\\2",gr$ALT))
    pos2 = as.integer(unlist(gsub(".*(\\]|\\[).+:([0-9]+).+","\\2",gr$ALT)))
    start2 = pos2 - flanking_bases
    end2 = pos2 + flanking_bases
    
    context1 = Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), names = chr, start =start1 , end = end1 , as.character = T )
    context2 = Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), names = chr2, start =start2 , end = end2 ,as.character = T )
    
    names(context1) <- paste(chr,start1,end1,sep="_")
    names(context2) <- paste(chr2,start2,end2,sep="_")
    
    gr$context1 = context1
    gr$context2 = context2
    
    return(gr)
  }
  
  gr_context <- get_context(gr)
  
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
  
  
  context2_heptamer_motif <- mapply(search_motif, gr_context$context2, as.character(heptamer_motif))
  context2_heptamer_rc_motif <- mapply(search_motif, gr_context$context2, as.character(heptamer_motif_rc))
  
  context2_nonamer_motif <- mapply(search_motif, gr_context$context2, as.character(nonamer_motif))
  context2_nonmaer_rc_motif <- mapply(search_motif, gr_context$context2, as.character(nonamer_motif_rc))
  
  info(header(vcf))["HM1",] = list(".","String","Heptamer motif positions (start-end) around POS")
  info(header(vcf))["HMRC1",] = list(".","String","Heptamer motif reverse complement positions (start-end) around POS")
  info(header(vcf))["NM1",] = list(".","String","Nonamer motif positions (start-end) around POS")
  info(header(vcf))["NMRC1",] = list(".","String","Nonamer motif reverse complement positions (start-end) around POS")
  
  info(header(vcf))["HM2",] = list(".","String","Heptamer motif positions (start-end) around ALT")
  info(header(vcf))["HMRC2",] = list(".","String","Heptamer motif reverse complement positions (start-end) around ALT")
  info(header(vcf))["NM2",] = list(".","String","Nonamer motif positions (start-end) around ALT")
  info(header(vcf))["NMRC2",] = list(".","String","Nonamer motif reverse complement positions (start-end) around ALT")
  
  info(vcf)$HM1 <- CharacterList(context1_heptamer_motif)
  info(vcf)$HMRC1 <- CharacterList(context1_heptamer_rc_motif)
  info(vcf)$NM1 <- CharacterList(context1_nonamer_motif)
  info(vcf)$NMRC1 <- CharacterList(context1_nonmaer_rc_motif)
  
  info(vcf)$HM2 <- CharacterList(context2_heptamer_motif)
  info(vcf)$HMRC2 <- CharacterList(context2_heptamer_rc_motif)
  info(vcf)$NM2 <- CharacterList(context2_nonamer_motif)
  info(vcf)$NMRC2 <- CharacterList(context2_nonmaer_rc_motif)
  
  outvcf <- file(gsub(".vcf$","_RSSmotif.vcf",vcf_fname), open="a")
  writeVcf(vcf, outvcf)
}
  
