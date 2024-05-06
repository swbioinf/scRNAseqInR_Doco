
bookdown::render_book()




# Code for boxes:

# 
# #### Why do we need to do this? {- .rational}
# Low quality cells can add noise to your results leading you to the wrong biological conclusions. Using only good quality cells helps you to avoid this. Reduce noise in the data by filtering out low quality cells such as dying or stressed cells (high mitochondrial expression) and cells with few features that can reflect empty droplets.
# 
# ####  {-}
# 
# 
# #### Challenge: The meta.data slot in the Seurat object { - .challenge}
# 
# Where are QC metrics stored in Seurat?
#   
#   * The number of unique genes and total molecules are automatically calculated during `CreateSeuratObject()`
# + You can find them stored in the object meta data
# 
# 1. What do you notice has changed within the `meta.data` table now that we have calculated mitochondrial gene proportion?
#   
#   2. Imagine that this is the first of 
# several samples in our experiment. Add a `samplename` column to to the `meta.data` table.
# 
# ####  {-}
