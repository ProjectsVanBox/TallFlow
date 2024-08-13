#!/bin/bash


### Change directory 
cd /hpc/pmc_vanboxtel/projects/TallClonal/2_Code/PreprocessingDataset/

### Load Guix
#source ~/"$GUIX_EXTRA_PROFILES"TallClonal/TallClonal/etc/profile
export PATH="$HOME/.config/guix/current/bin":$PATH
#source ~/.guix.conf 
guixr environment --pure --ad-hoc r r-ggplot2 r-seurat r-dplyr r-annotationdbi r-dbi r-singlecellexperiment r-scales r-cowplot r-rcurl r-org-hs-eg-db r-matrixstats r-ggsci r-biocparallel r-biocneighbors r-biocmanager r-celldex r-batchelor r-seuratdisk r-singler r-tidyverse python-umap-learn r-uwot r-go-db -- <<EOF
	Rscript Preprocessing_ThyComb_NoAb_V1.R
	exit
EOF
#r-umap
