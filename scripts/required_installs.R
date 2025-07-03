
# Packages from Cran
install.packages("tidyverse")
install.packages("Seurat")
install.packages("patchwork")
install.packages("hdf5r") #not really needed, only used to prepare 10Xpbmc.
install.packages("clustree")
install.packages("harmony") #formerly dev version- devtools::install_github("immunogenomics/harmony")

# packages from github
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github('satijalab/seurat-data')


# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("SingleCellExperiment")
BiocManager::install("HDF5Array")
BiocManager::install("SingleR")
BiocManager::install("celldex")  
BiocManager::install("limma")
BiocManager::install("edgeR")


#---------------------------
# Some loading of datasets into cache. Downloads data!
# Less important, only needed for cell typing section
# Run these R commands.
library(celldex)
celldex::HumanPrimaryCellAtlasData()
