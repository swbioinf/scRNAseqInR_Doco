
# First an up to date version of R (the most recent will do) - 4.2+


# Cran
install.packages("tidyverse")
install.packages("Seurat")
install.packages("patchwork")
install.packages("hdf5r") #not really needed, only used to prepare 10Xpbmc.
install.packages("clustree")

# github
install.packages("remotes")
install.packages("devtools")
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github('satijalab/seurat-data')

devtools::install_github("immunogenomics/harmony")
devtools::install_github("powellgenomicslab/scPred")

# And bioconductor (again, the most recent) - v 3.9 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("SingleCellExperiment")
BiocManager::install("HDF5Array")
BiocManager::install("SingleR")
BiocManager::install("celldex")   #download????
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("CHETAH") 




#---------------------------
# Some loading of datasets into cache. 
# Less important.
# Run these R commands.
library(celldex)
celldex::HumanPrimaryCellAtlasData()

# only for setup
#library(SeuratData)
#data("hcabm40k")  # COuld be replaced with a preset object





#---------------------------
# And monocle.
# We probably will *not* use monocle this time. 
# However, if its easy, it would be nice to have it there (maybe reuse the image for a future workshop)
# If its trivial, could it be installed please? 
# If it isn't trivial, then don't bother.
if (FALSE) {
  # Monocle has a depenancy on the gdal tool (via package sf), which needs to be installed before the library.
  # See:   https://cole-trapnell-lab.github.io/monocle3/docs/installation/
  
  # Install gdal
  # apt-get install libgdal-dev
  #
  # 
  # Then
  devtools::install_github('cole-trapnell-lab/monocle3')
}


