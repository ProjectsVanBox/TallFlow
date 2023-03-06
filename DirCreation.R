#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)


setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/")
### Create main directories 
dir.create("1_Input")
dir.create("2_Code")
dir.create("3_Output")


### Create input directories
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/")
dir.create("ThresholdFiles")
dir.create("FCS_Thymus")
dir.create("FCS_Patients")


### Create code directories
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/2_Code/")
dir.create("PythonGating")
dir.create("PythonGating/Thy_scripts")
dir.create("PythonGating/Pts_scripts")
dir.create("Barcharts")
dir.create("Linegraphs")
dir.create("Umaps")
dir.create("UmapsSample")


### Create output directories
setwd("/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/")
dir.create("PythonGating")
dir.create("PythonGating/Thymus")
dir.create("PythonGating/Patients")
dir.create("Barcharts")
dir.create("Barcharts/Thymus")
dir.create("Barcharts/Patients")
dir.create("Linegraphs")
dir.create("Umaps")
dir.create("Umaps/Thymus")
dir.create("Umaps/Thymus/MarkerPlots")
dir.create("Umaps/Thymus/Consensus_plots")
dir.create("Umaps/Thymus/SingleThymusUmap")
dir.create("Umaps/Patients")
dir.create("Umaps/Patients/MarkerPlots")
dir.create("Umaps/Patients/SinglePatientUmap")
dir.create("Umaps/ThymusPatients")
dir.create("ThymusPatients/MarkerPlots")




