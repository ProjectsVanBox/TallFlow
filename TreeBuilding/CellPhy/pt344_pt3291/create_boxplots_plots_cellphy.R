### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(VariantAnnotation)

setwd("/Users/v.m.poort/hpc/projects/Tsites/3_Output/TreeBuilding/PB09204/CellPhyWrapper_v2/branches_vcf/")

# colors for plotting
### blue = BM #1A3E77
## orange = LN # ED7402
## third color = purple: #D484B6
bulk_names = c("BM-LYMBULK", "LN-LYMBULK")
colors = c("#1A3E77", "#ED7402")
names(colors) = c("BM-LYMBULK", "LN-LYMBULK")
chrom_del = "^1:|^7:|X|Y" #large deletions or LOH in these chromosomes excluded for VAF calculations 
end_branch = c("D_branch", "C_branch", "B_branch", "A_branch", "H_branch", "J_branch", "U_branch", 
               "T_branch", "W_branch", "R_branch", "M_branch", "O_branch")
rootbranch = "X_branch"



combined_df <- data.frame()
for (vcf_fname in list.files(path='.',pattern='_branch.vcf',full.names = T)) {
  branch_name <- gsub(".vcf","", gsub("./", "", vcf_fname))
  print(vcf_fname)
  myvcf <- readVcf(vcf_fname)
  if (nrow(myvcf) < 5 ) {
    next
  }
  mydf <- data.frame()
  for (sname in samples(header(myvcf))) {
    ad <- sapply(geno(myvcf)$AD[,sname],"[[",2)
    dp <- geno(myvcf)$DP[,sname]
    vaf <- ad/dp
    sname <- gsub(".+DX1TLBL","\\1",sname)
    mydf <- rbind( mydf, cbind( sname, as.numeric(vaf), names(vaf)) )
  }
  print(nrow(mydf))
  colnames(mydf) <- c("bulk","vaf", "variant")
  mydf$vaf <- as.numeric(mydf$vaf)
  mydf <- mydf %>% filter(bulk %in% bulk_names)
  mydf$branch <- branch_name
  
  combined_df <- rbind(combined_df, mydf)
}

# only for the root -> per site -> check mean VAF. -> should be 0.5 -> calculate tumor percentage. 
# correct all other VAFs for this tumor percentage. 
# filter for deleted chromosomes 
combined_df2 <- combined_df[!grepl(chrom_del,combined_df$variant),]
#nrow(combined_df2)
#nrow(combined_df)

combined_df <- combined_df2

# replace all 0 with NA so that the median is not too low. 
is.na(combined_df$vaf) <- combined_df$vaf==0


median_vaf <- combined_df %>% filter(branch == rootbranch ) %>% 
  group_by(bulk) %>% summarise(median = median(vaf, na.rm = TRUE))

# correction factor? 
median_vaf$tumor_perc <- median_vaf$median/0.005
median_vaf$corf <- 100/median_vaf$tumor_perc

Bm_corf <- as.numeric(median_vaf$corf[1])
Ln_corf <- as.numeric(median_vaf$corf[2])


# correct VAFs
combined_df1 <- combined_df %>% mutate(vaf_cor = if_else(bulk == "BM-LYMBULK", vaf*Bm_corf, vaf*Ln_corf))
unique(combined_df1$branch)
# set maximus VAF value to 1
combined_df1$vaf_cor2 <- pmin(1, combined_df1$vaf_cor)

# check after VAF correction which values are > 1 -> likely only 1 allele detected so set to original vaf to corrected/2.
combined_df1$vaf_cor3 <- combined_df1$vaf_cor
combined_df1 <- na.omit(combined_df1)
#combined_df1[is.na(combined_df1)] <- 0 # use if the 0's need to be plotted
combined_df1$vaf_cor3 <- replace(combined_df1$vaf_cor3, combined_df1$vaf_cor3>=1, combined_df1$vaf_cor3/2)
head(combined_df1$vaf_cor3)

# get rid of the end branches. 
combined_df2 <- combined_df1 %>% filter(!branch %in% end_branch)
unique(combined_df2$branch)

# plot
p <- ggplot( combined_df2, aes(x=bulk, y=vaf_cor, fill=bulk)) +
  geom_violin() +
  geom_boxplot(width=0.1)+
  facet_wrap(~ branch) + 
  theme_classic() +
  scale_fill_manual(values= alpha(colors, 0.3))

p2 <- ggplot( combined_df2, aes(x=bulk, y=vaf_cor2, fill=bulk)) +
  geom_violin() +
  geom_boxplot(width=0.1)+
  facet_wrap(~ branch) + 
  theme_classic() +
  scale_fill_manual(values= alpha(colors, 0.3))

p3 <- ggplot( combined_df2, aes(x=bulk, y=vaf_cor3, fill=bulk)) +
  geom_violin() +
  geom_boxplot(width=0.1)+
  facet_wrap(~ branch) + 
  theme_classic() +
  scale_fill_manual(values= alpha(colors, 0.3))

pdf("bulkcontributions_cor.pdf")
print(p)
print(p2)
print(p3)
dev.off()


# filter branches with < 50 mutations
combined_df3 <- combined_df1 %>% filter(!branch %in% c(end_branch, "E_branch", "G_branch"))

unique(combined_df3$branch)
combined_df3$realname <- plyr::revalue(as.character(combined_df3$branch), 
              c("F_branch" = 'Branch B: 56',
                "P_branch" = 'Branch A: 176',
                "X_branch" = 'Root: 1151',
                "N_branch" = 'Branch C: 59'))

branch_order = c("Root: 1151", "Branch A: 176", "Branch B: 56", "Branch C: 59")

p4 <- ggplot( combined_df3, aes(x=bulk, y=vaf_cor3, fill=bulk)) +
  geom_violin() +
  geom_boxplot(width=0.1)+
  facet_wrap(~factor(realname, levels = branch_order), ncol = 2) + 
  theme_classic() +
  scale_fill_manual(values= alpha(colors, 0.3)) +
  ylab("Corrected VAF")

pdf("bulkcontributions_cor_forpaper.pdf")
print(p4)
dev.off()

