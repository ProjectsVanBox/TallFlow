
library(ggplot2)
library(ggpubr)

mydf <- read.table("~/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/RSSmotifSearch/rss_count.txt")

colnames(mydf) <- c("file","rss1","rss2","rss_both","svs")

mydf$pt <- gsub("(pt.+?)-.+","\\1",mydf$file)
mydf$clone <- gsub("pt.+-.+-(T*ALL.+?)\\..+","\\1",mydf$file)
mydf$label <- "bulk"
mydf[grepl("TALL",mydf$clone),]$label <- "single_clone"
mydf$perc <- mydf$rss_both/mydf$svs


ggplot(data=mydf,aes(x=label,y=perc)) +
  geom_boxplot() +
  facet_grid(.~pt)


p <- ggboxplot(mydf, x = "label", y = "perc",
               color = "label", palette = "jco",
               add = "jitter") +
  facet_grid(.~pt)
#  Add p-value
p + stat_compare_means()
