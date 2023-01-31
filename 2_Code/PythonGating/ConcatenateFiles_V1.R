#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


### Load libraries
library(dplyr)
library(tidyr)
#library(data.table)



setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow")
### Loop over file and select columns, add sample name 
file.list <-  list.files(path = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/populationexports", pattern = "*.csv$", full.names = TRUE, recursive = TRUE)
combine.list <- list()
prev.samplename = ""
for (f in file.list){
  print(f)
  # read in file and add cell type 
  celltype <- gsub(x = gsub(x = f, pattern = ".*__", ""), pattern = ".csv", replacement = "")
  sample.name <- gsub(x = gsub(x = gsub(x = f, pattern = "/.*/", replacement = ""), pattern = "export_export_", replacement = "export_"), pattern = "_stained.*", replacement = "_stained_CD45+.csv")
  data.file <- read.table(file = f, sep = ",", header = TRUE, row.names = NULL)
  data.file$ManualGate <- rep.int(celltype, length(rownames(data.file)))
  
  # Combine the sample and write the output
  if (sample.name != prev.samplename){
    if (! prev.samplename == ""){
      combine.file <- do.call("rbind", combine.list)
      write.table(x = combine.file, paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_manGate_Thymi/", prev.samplename, ".csv"), sep = ",", col.names = NA)
      prev.samplename = sample.name
      combine.list <- list()
      combine.list[[paste0(sample.name, "_", celltype)]] <- data.file
    }else{
      prev.samplename = sample.name
      combine.list[[paste0(sample.name, "_", celltype)]] <- data.file
    }
  }else{
    combine.list[[paste0(sample.name, "_", celltype)]] <- data.file
  }
}
# Get the last one 
combine.file <- do.call("rbind", combine.list)
write.table(x = combine.file, paste0("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_manGate_Thymi/", prev.samplename, ".csv"), sep = ",", col.names = NA)

