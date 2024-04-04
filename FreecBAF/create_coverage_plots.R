library(ggplot2)

setwd("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/FreecBAF/coverage/")

msc_cov <- read.table("pt2229-DX1BM-MSCBULK_dedup.cov", fill = T, sep="\t")
msc_cov$TCR <- substr(msc_cov$V4,1,3)
msc_cov <- msc_cov[order(msc_cov$V1, msc_cov$V2),]
msc_cov <- msc_cov[msc_cov$TCR != "",]
msc_cov$NORM <- msc_cov$V5/mean(msc_cov$V5)

p1 <- ggplot() +
  geom_point(data=msc_cov,aes(x=V4,y=NORM)) +
  facet_wrap(TCR ~ ., scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0,2)

for ( myfile in list.files(path=".",pattern = ".cov")) {
  my_cov <- read.table(myfile, fill = T, sep="\t")
  my_cov$TCR <- substr(my_cov$V4,1,3)
  my_cov <- my_cov[order(my_cov$V1, my_cov$V2),]
  my_cov <- my_cov[my_cov$TCR != "",]
  my_cov$NORM <- my_cov$V5/mean(my_cov$V5)
  
  p2 <- p1 + geom_point(data=my_cov,aes(x=V4,y=NORM),col="red")
  
  pdf(gsub(".cov",".pdf",myfile))
  print(p2)
  dev.off()
}