
library(VariantAnnotation)
library(ggplot2)

myBulkVcf_fname <- "~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2229//intermediate/short_variants/somatic_vcfs/pt2229/220919_HMFreg1764_pt2229_pt2229-DX1BM-ALLBULK.SMuRF.filtered.vcf.gz"

myBulkVcf <- readVcf(myBulkVcf_fname)

absentTable <- as.data.frame(table(info(myBulkVcf)$ABSENT_SAMPLES))
subclonalTable <- as.data.frame(table(info(myBulkVcf)$SUBCLONAL_SAMPLES))
clonalTable <- as.data.frame(table(info(myBulkVcf)$CLONAL_SAMPLES))

ggplot() +
  geom_point(data=clonalTable, aes(x=as.numeric(Var1),y=Freq),col='black') +
  geom_point(data=subclonalTable, aes(x=as.numeric(Var1),y=Freq),col='blue') +
  geom_point(data=absentTable, aes(x=as.numeric(Var1),y=Freq),col='red')

which(info(myBulkVcf)$CLONAL_SAMPLES + info(myBulkVcf)$SUBCLONAL_SAMPLES <= 2)
which(info(myBulkVcf)$CLONAL_SAMPLES == 1)

plot_df <- data.frame( N=c(1:13), Var=0)
for ( x in c(1:13)) {
  plot_df[x,]$Var <- sum(clonalTable[1:x,]$Freq)
}
plot_df$Rate <- round((1-plot_df$Var/nrow(myBulkVcf))*100,digits=2)
plot_df$Var2 <- nrow(myBulkVcf)-plot_df$Var

pdf("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/TruePosRate/pt2229/freqplot.pdf")
ggplot(plot_df,aes(x=N,y=Var2)) +
  geom_bar(stat="identity") +
  geom_text(aes(x=N,y=Var2,label=Rate),position=position_dodge(width=0.9), vjust=-0.25)
dev.off()