library(flowCore)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowAI")

### -------------------------------- ### Quality assessment ### -------------------------------- ### 
InpuDir_Thymus10X <- "/Users/verapoort/surfdrive/Shared/pmc_vanboxtel/projects/TallFlow/1_Input/Thymus_samples/20200615_VP_Thy072_COG5_R/"
setwd(InpuDir_Thymus10X)

InputDir <- InpuDir_Thymus10X

fcsfiles <- dir(InpuDir_Thymus10X, pattern="*fcs$")
for (fcsfile in fcsfiles){
  print(fcsfile)
  resQC <- flow_auto_qc(fcsfile, remove_from  = "FR_FM", output = 1) # using a character vector
  write.FCS(resQC, paste0(InputDir, "/QC_files_fcs_files"))
}
#rm(resQC)
