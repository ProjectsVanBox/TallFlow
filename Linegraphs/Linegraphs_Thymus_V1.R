#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
#BiocManager::install("")

### Load data  
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Umaps/Thymus/")
all.df <- readRDS("./CombinedThymi.rds")



median.data <- all.df %>% group_by(Celltype, SampleID) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), median)
median.data <- median.data[! median.data$Celltype %in% c("B-cell", "Monocyte"),]
plot.data <- as.data.frame(sapply(median.data[3:19], function(x) ( (x-min(x))/(max(x)-min(x))) ))
plot.data$Celltype <- median.data$Celltype
plot.data$SampleID <- median.data$SampleID

CD4.data <- plot.data[,c("CD4", "SampleID", "Celltype")]
CD4.data$Celltype <- as.factor(CD4.data$Celltype)
CD4.data_long <- melt(CD4.data, id="Celltype") 


ggplot(CD4.data_long, aes(x=Celltype, y=value)) + 
  geom_line() +
  geom_point()



#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


df2 <- data_summary(ToothGrowth, varname="len", 
                    groupnames=c("supp", "dose"))
# Convert dose to a factor variable
df2$dose=as.factor(df2$dose)

ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(0.05))


df3 <- data_summary(CD4.data, varname="CD4", 
                    groupnames=c("SampleID", "Celltype"))
# Convert dose to a factor variable
df3$Celltype=as.factor(df3$Celltype)
ggplot(df3, aes(x=Celltype, y=CD4, group=SampleID, color=SampleID)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=CD4-sd, ymax=CD4+sd), width=.2,
                position=position_dodge(0.05))



p<- ggplot(CD4.data, aes(x=Celltype, y=CD4, group=SampleID, color=SampleID)) + 
  geom_line() +
  geom_point()#+
  #geom_errorbar(aes(ymin=CD4-sd, ymax=CD4+sd), width=.2,
  #              position=position_dodge(0.05))

print(p)

#melt(setDT(wide), id.vars = c("Code","Country"), variable.name = "year")
test_long <- melt(CD4.data, id="Celltype")  # convert to long format

ggplot(data=test_long,
       aes(x=Celltype, y=value, colour=variable)) +
  geom_line()

tg <- ToothGrowth

#test <- melt(setDT(CD4.data), id.vars = c("Celltype","SampleID"), variable.name = "Marker")
ggplot(data=CD4.data, aes(x = CD4.data)) +
  geom_line()

ggplot(data=plot.data,
       aes(x=Celltype, y=CD4)) +
  geom_line()

ggplot(data=plot.data,
       aes(x=Celltype, y=CD4)) +
  geom_line(aes(group = Celltype), stat = "summary", fun = mean)

### Working one-ish
ggplot(data=plot.data,
       aes(x=Celltype, y=CD4)) +
  geom_line(aes(group = SampleID), stat = "summary", fun = median) + 
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)


ggplot(data=plot.data,
       aes(x=Celltype)) +
  geom_line(aes(y=CD4, group = SampleID), stat = "summary", fun = median) 



ggplot(plot.data, aes(Celltype, CD4)) +
  geom_line(aes(group = 1)) +
  geom_errorbar( aes(ymin = CD4-sd, ymax = CD4+sd),width = 0.2) +
  geom_point(size = 2)


ggplot(data=plot.data,
       aes(x=Celltype, y=c("CD4", "CD8a"))) +
  geom_line(aes(group = SampleID), stat = "summary", fun = median) + 
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)

ggplot(data=plot.data,
       aes(x=Celltype, y=CD4)) +
  geom_line(aes(group = Celltype), stat = "summary", fun = mean)

ggplot(data=plot.data,
       aes(x=Celltype, y=CD4, group = SampleID)) +
  geom_line() + 
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)


LowQualPlot <- plot.data[,1:17]
plot(CD4 ~ rownames(LowQualPlot), 
     data = LowQualPlot, 
     type = "b",
     xaxt = "n",
     ylim = c(min(LowQualPlot[ ,-1]), max(LowQualPlot[ ,-1])),
     lwd = 2,
     ylab = "Normalised intensity",
     xlab = NA)



long <- melt(setDT(plot.data), id.vars = c("SampleID","Celltype"), variable.name = "Marker")
ggplot(data=long,
       aes(x=Celltype, y=SampleID, colour=Marker)) +
  geom_line()
ggplot(long, aes(Celltype, SampleID, colour = Marker)) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = Celltype), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)




plot.data3 <- melt(plot.data, id="Celltype")
ggplot(data=plot.data3,
       aes(x=Celltype, y=value, colour=variable)) +
  geom_line()
ggplot(median.data, aes(x=Celltype)) + 
  geom_line(aes(y = CD4), color = "darkred") + 
  geom_line(aes(y = CD8a), color="steelblue", linetype="twodash") 


plot.data2 <- tidyr::pivot_longer(plot.data, c("SampleID"))

ggplot(plot.data2, aes(Celltype, CD4, colour = Celltype)) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = name), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)


ggplot(plot.data2, aes(Celltype, CD4, colour = Celltype)) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = name), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)


ggplot(plot.data2, aes(Celltype, CD4, colour = name)) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = name), stat = "summary", fun = mean) +
  geom_line(aes(group = name), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)
ggplot(plot.data2, aes(Celltype, CD4, colour = SampleID)) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = name), stat = "summary", fun = mean) +
  geom_line(aes(group = name), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)




ggplot(plot.data, aes(Celltype, CD4, colour = Celltype)) +
  #geom_jitter(width = 0.2) +
  geom_line(aes(group = SampleID), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)

ggplot(plot.data, aes(Celltype, CD4, colour = CD4)) +
  #geom_jitter(width = 0.2) +
  geom_line(aes(group = SampleID), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)



df <- read.table(text = "Time Replicate      Media1     Media2
1       T0        R1 0.013733333 0.03770000
2       T0        R2 0.008233333 0.03690000
3       T0        R3 0.013333333 0.03760000
4       T1        R1 0.018166667 0.03680000
5       T1        R2 0.013466667 0.03570000
6       T1        R3 0.018366667 0.03700000
7       T2        R1 0.028066667 0.04420000
8       T2        R2 0.022966667 0.04400000
9       T2        R3 0.027866667 0.04420000
10      T3        R1 0.041333333 0.05536667
11      T3        R2 0.033433333 0.05476667
12      T3        R3 0.039833333 0.05446667")
df
df <- tidyr::pivot_longer(df, c("Media1", "Media2"))
df
ggplot(df, aes(Time, value, colour = name)) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = name), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)



median.data <- all.df %>% group_by(all.df$Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), median)
colnames(median.data) <- c( "Celltype", "CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d")
median.data <- median.data[! median.data$Celltype %in% c("B-cell", "Monocyte"),]
calculation.data <- median.data %>% dplyr:: select(-starts_with("Celltype"))
plot.data <- as.data.frame(sapply(calculation.data, function(x) ( (x-min(x))/(max(x)-min(x))) ))


plot(CD4 ~ rownames(plot.data), 
     data = plot.data, 
     type = "b",
     xaxt = "n",
     ylim = c(min(plot.data[ ,-1]), max(plot.data[ ,-1])),
     lwd = 2,
     ylab = "Normalised intensity",
     xlab = NA)
axis(side = 1, at = rownames(plot.data), labels = median.data$Celltype, las = 2)
lines(rownames(plot.data), plot.data$CD8, col = "steelblue", lwd = 2, type = "b")
lines(rownames(plot.data), plot.data$CD3, col = "Pink", lwd = 2, type = "b")
lines(rownames(plot.data), plot.data$CD5, col = "Red", lwd = 2, type = "b")
lines(rownames(plot.data), plot.data$CD45, col = "Green", lwd = 2, type = "b")
lines(rownames(plot.data), plot.data$CD7, col = "Yellow", lwd = 2, type = "b")
lines(rownames(plot.data), plot.data$CD1a, col = "Orange", lwd = 2, type = "b")
lines(rownames(plot.data), plot.data$TCRa_b, col = "Purple", lwd = 2, type = "b")
lines(rownames(plot.data), plot.data$TCRg_d, col = "Gold", lwd = 2, type = "b")
legend(7.5, 1, legend=c("CD45", "CD7", "CD5", "CD1a", "CD3", "TCR_GD", "CD4", "CD8",  "TCR_AB"),
       col=c("Green", "Yellow", "Red", "Orange", "Pink", "Gold", "Black", "steelblue", "Purple"), 
       lty=1:0.5, cex=0.6, bg = NA)




















#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(plot.data, varname="CD4", 
                    groupnames=c("SampleID", "Celltype"))
# Convert dose to a factor variable
df2$Celltype=as.factor(df2$Celltype)



p<- ggplot(df2, aes(x=Celltype, y=CD4, group=SampleID, color=SampleID)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=CD4-sd, ymax=CD4+sd), width=.2,
                position=position_dodge(0.05))
print(p)



mean(plot.data[['CD4']], na.rm=TRUE)
sd(plot.data[["CD4"]], na.rm=TRUE)

df2 <- data_summary(all.df, varname="CD4", 
                    groupnames=c("SampleID", "Celltype"))
# Convert dose to a factor variable
df2$Celltype=as.factor(df2$Celltype)

p<- ggplot(df2, aes(x=Celltype, y=CD4, group=SampleID, color=SampleID)) + 
  geom_line() +
  #geom_point()+
  geom_errorbar(aes(ymin=CD4-sd, ymax=CD4+sd), width=.2,
                position=position_dodge(0.05))
print(p)

#all.df



df.summary <- all.df %>%
  group_by(SampleID, Celltype) %>%
  summarise(
    sd = sd(CD4, na.rm = TRUE),
    CD4 = median(CD4)
  )
df.summary

ggplot(df.summary, aes(Celltype, CD4)) +
  geom_line(aes(group = 1)) +
  geom_errorbar( aes(ymin = CD4-sd, ymax = CD4+sd),width = 0.2) +
  geom_point(size = 2)





df <- ToothGrowth
df$dose <- as.factor(df$dose)
head(df, 3)
df.summary <- df %>%
  group_by(dose) %>%
  summarise(
    sd = sd(len, na.rm = TRUE),
    len = mean(len)
  )
df.summary
ggplot(df.summary, aes(dose, len)) +
  geom_line(aes(group = 1)) +
  geom_errorbar( aes(ymin = len-sd, ymax = len+sd),width = 0.2) +
  geom_point(size = 2)








## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}
## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}


tg <- ToothGrowth
head(tg)
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
tgc <- summarySE(tg, measurevar="len", groupvars=c("supp","dose"))
tgc

dfw <- read.table(header=TRUE, text='
 subject pretest posttest
       1    59.4     64.5
       2    46.4     52.4
       3    46.0     49.7
       4    49.0     48.7
       5    32.5     37.4
       6    45.2     49.5
       7    60.3     59.9
       8    54.3     54.1
       9    45.4     49.6
      10    38.9     48.5
 ')

# Treat subject ID as a factor
dfw$subject <- factor(dfw$subject)
# Convert to long format
library(reshape2)
dfw_long <- melt(dfw,
                 id.vars = "subject",
                 measure.vars = c("pretest","posttest"),
                 variable.name = "condition")

dfw_long
dfwc <- summarySEwithin(dfw_long, measurevar="value", withinvars="condition",
                        idvar="subject", na.rm=FALSE, conf.interval=.95)

dfwc
# Make the graph with the 95% confidence interval
ggplot(dfwc, aes(x=condition, y=value, group=1)) +
  geom_line() +
  geom_errorbar(width=.1, aes(ymin=value-ci, ymax=value+ci)) +
  geom_point(shape=21, size=3, fill="white") +
  ylim(40,60)





med_long <- melt(median.data,
                 id.vars = "SampleID",
                 measure.vars = c("CD4","CD8a"),
                 variable.name = "condition")
med_long <- summarySEwithin(med_long, measurevar="value", withinvars="condition",
                        idvar="SampleID", na.rm=FALSE, conf.interval=.95)
ggplot(med_long, aes(x=condition, y=value, group=1)) +
  geom_line() +
  geom_errorbar(width=.1, aes(ymin=value-ci, ymax=value+ci)) +
  geom_point(shape=21, size=3, fill="white") #+
  #ylim(40,60)















######################################################################## test old fashioned summarized plots 
#median.data_old 
stage <- c("DN", "gdTcell", "ISP", "DP", "CD4", "CD8", "NK-cell", "pDC", "cDC", "Unknown")
median.data_old <- all.df %>% group_by(Celltype) %>% summarise_at(vars("CD25", "CD56", "CD4", "TCRa_b", "CD1a", "CD19", "CD14", "CD8a", "CD7", "CD123", "CD44", "CD3", "CD11c", "CD5", "CD117", "CD45", "TCRg_d"), median)
median.data_old <- median.data_old[! median.data_old$Celltype %in% c("B-cell", "Monocyte"),]
plot.data_old <- as.data.frame(sapply(median.data_old[3:18], function(x) ( (x-min(x))/(max(x)-min(x))) ))
plot.data_old$Celltype <- factor(median.data_old$Celltype, levels = stage)
#plot.data_old$SampleID <- median.data_old$SampleID
plot.data_old

#plot.data_old <- plot.data_old[match(stage, plot.data_old$Celltype), ]
#plot.data_old$Celltype <- as.factor(plot.data_old$Celltype)


colors = c("CD45"="Green", 
           "CD7"="Yellow", 
           "CD5"="Red", 
           "CD1a"="Orange", 
           "CD3"="Pink",
           "TCRg_d"="Gold",
           "CD4"="Black",
           "CD8a"="steelblue",
           "TCRa_b"="Purple")
pdf("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/Linegraphs/CombinedThymus_line_simple.pdf")
ggplot(plot.data_old, aes(x=Celltype), color = colors) + 
  geom_line(aes(y = CD4, group = 1, color = "CD4")) +
  geom_line(aes(y = CD8a, group = 1, color = "CD8a")) +
  geom_line(aes(y = CD3, group = 1, color = "CD3")) +
  geom_line(aes(y = CD5, group = 1, color = "CD5")) +
  geom_line(aes(y = CD45, group = 1, color = "CD45")) +
  geom_line(aes(y = CD7, group = 1, color = "CD7")) +
  geom_line(aes(y = CD1a, group = 1, color = "CD1a")) +
  geom_line(aes(y = TCRa_b, group = 1, color = "TCRa_b")) +
  geom_line(aes(y = TCRg_d, group = 1, color = "TCRg_d")) + 
  #geom_point() +
  labs(x = "Cell type", y = "Relative expression", color = "Legend") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()
  


colors = c("CD4" = "black",
           "CD8a" = "blue")
ggplot(plot.data_old, aes(x=Celltype), color = colors) + 
  geom_line(aes(y = CD4, group = 1, color = "CD4")) +
  geom_line(aes(y = CD8a, group = 1, color = "CD8a")) +
  scale_color_manual(values=colors)

  scale_colour_manual("", 
                      breaks = c("CD45", "CD7", "CD5", "CD1a", "CD3", "TCR_GD", "CD4", "CD8",  "TCR_AB"),
                      values=c()) #+ 
  labs(x = "Cell type", y = "Relative expression") #+
  #theme_bw()



ggplot(plot.data, aes(x=Celltype)) + 
  geom_line(aes(y = CD4, group = 1)) +
  geom_line(aes(y = CD8a, group = 1)) 


ggplot(data=plot.data,
       aes(x=Celltype, y=CD4)) +
  geom_line(aes(group = SampleID), stat = "summary", fun = median) + 
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1)
