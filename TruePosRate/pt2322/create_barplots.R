
library(VariantAnnotation)
library(ggplot2)
library(reshape2)

root_vcf_fname <- "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TreeBuilding/pt2322/withoutBulk/branches_vcfs/root_branch.vcf"

# sampleName <- "TALL1"

create_my_root_df <- function(sampleName ) {
  snv_ptato_vcf_fname <- paste("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/snvs/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-",sampleName,"/220919_HMFreg1764_pt2322_pt2322-DX1BM-",sampleName,".ptato.filtered.vcf.gz",sep="")
  indel_ptato_vcf_fname <- paste("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2322/indels/pt2322/220919_HMFreg1764_pt2322_pt2322-DX1BM-",sampleName,"/220919_HMFreg1764_pt2322_pt2322-DX1BM-",sampleName,".indels.ptato.filtered.vcf.gz",sep="")
  
  callable_fname <- paste("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TruePosRate/pt2322/",sampleName,"_root_Callable.txt",sep="")
  
  root_vcf <- readVcf(root_vcf_fname)
  root_vcf_df <- as.data.frame(granges(root_vcf))
  root_vcf_df$SMuRF <- "NA"
  clonal_rows <- which( unlist(lapply(grep(paste(sampleName,"$",sep=""),info(root_vcf)$CLONAL_SAMPLE_NAMES),length)) == 1)
  subclonal_rows <- which( unlist(lapply(grep(paste(sampleName,"$",sep=""),info(root_vcf)$SUBCLONAL_SAMPLE_NAMES),length)) == 1)
  absent_rows <- which( unlist(lapply(grep(paste(sampleName,"$",sep=""),info(root_vcf)$ABSENT_SAMPLE_NAMES),length)) == 1)
  root_vcf_df$SMuRF[clonal_rows] <- "CLONAL"
  root_vcf_df$SMuRF[subclonal_rows] <- "SUBCLONAL"
  root_vcf_df$SMuRF[absent_rows] <- "ABSENT"
  
  root_vcf_df$QC <- "NA"
  passQC_rows <- which( unlist(lapply(grep(paste(sampleName,"$",sep=""),info(root_vcf)$PASS_QC_SAMPLE_NAMES),length)) == 1)
  failQC_rows <- which( unlist(lapply(grep(paste(sampleName,"$",sep=""),info(root_vcf)$FAIL_QC_SAMPLE_NAMES),length)) == 1)
  root_vcf_df$QC[passQC_rows] <- "PASS"
  root_vcf_df$QC[failQC_rows] <- "FAIL"
  
  my_root_df <- root_vcf_df[,c("seqnames","start","SMuRF","QC")]
  
  snv_ptato_vcf <- readVcf(snv_ptato_vcf_fname)
  snv_ptato_df <- as.data.frame(granges(snv_ptato_vcf))
  snv_ptato_df$PTATO <- "True"
  
  indel_ptato_vcf <- readVcf(indel_ptato_vcf_fname)
  indel_ptato_df <- as.data.frame(granges(indel_ptato_vcf))
  indel_ptato_df$PTATO <- "True"
  
  my_ptato_df <- rbind(snv_ptato_df,indel_ptato_df)[,c("seqnames","start","PTATO")]

  my_root_df <- merge(my_root_df, my_ptato_df, by=c("seqnames","start"), all.x=T)
  my_root_df$PTATO[is.na(my_root_df$PTATO)] <- "False"
  
  callable_df <- read.table(callable_fname)
  colnames(callable_df) <- c("seqnames","start","chr","bin_start","bin_end","CALLABLE")
  callable_df <- callable_df[,c("seqnames","start","CALLABLE")]
  
  my_root_df <- merge(my_root_df, callable_df, by=c("seqnames","start"), all.x=T)
  my_root_df$sample <- sampleName

  return(my_root_df)
}

my_sample_list <- list("TALL1","TALL2","TALL3","TALL4","TALL5","TALL6","TALL7","TALL8","TALL9","TALL10","TALL11","TALL12")

my_root_df_list <- lapply(my_sample_list, create_my_root_df)

my_df <- do.call("rbind", my_root_df_list)

create_plot <- function( variable ) {

  my_variable_df <- as.data.frame(table(my_df$sample,my_df[[variable]]))

  p <- ggplot(data=my_variable_df, aes(x=Var1, y=Freq, fill=Var2)) +
    geom_bar(stat="identity",)
  
  print(p)
}

my_variable_list <- list("CALLABLE","SMuRF","QC","PTATO")

lapply(my_variable_list, create_plot)

my_df2 <- my_df[my_df$SMuRF == "CLONAL" & my_df$QC == "PASS",]
my_variable_df2 <- as.data.frame(table(my_df2$sample,my_df2[["PTATO"]]))

p2 <- ggplot(data=my_variable_df2, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity")

print(p2)

p3 <- ggplot(data=my_variable_df2, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity",position = "fill")

print(p3)



