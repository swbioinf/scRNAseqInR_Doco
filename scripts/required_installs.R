# For windows only:
install.packages("RTools")

# Packages from Cran
install.packages("tidyverse")
install.packages("Seurat")
install.packages("patchwork")
install.packages("hdf5r") #not really needed, only used to prepare 10Xpbmc.
install.packages("clustree")
#install.packages("harmony") #formerly dev version- devtools::install_github("immunogenomics/harmony"). Only used for harmony section, not in workshop


# packages from github
install.packages("remotes")
# If seurat-wrappers gives errors about Matrix package,
# Try to install this specific version, and retry:
#  remotes::install_version("Matrix", version = "1.6-5")
remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github('satijalab/seurat-data')
remotes::install_github('immunogenomics/presto') #optional, speeds up clustering



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
