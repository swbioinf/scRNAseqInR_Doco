library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)


##  3.1.2 Setup the Seurat Object 

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/") # This path on our VMs!
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc



##3.2 QC Filtering

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2




# Challenge: filter the cells
pbmc10k_unfiltered <- readRDS("data/10k_PBMC_v3.1ChromiumX_Intronic.rds")
VlnPlot(pbmc10k_unfiltered, features = "nCount_RNA") + scale_y_log10()
VlnPlot(pbmc10k_unfiltered, features = "nFeature_RNA") + scale_y_log10()


##Normalizing the data
pbmc <- NormalizeData(pbmc)


## Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2



## Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


## Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# Actually show the data here, 
DimPlot(pbmc, reduction = "pca")



DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

FeaturePlot(pbmc, reduction = "pca", "LYZ")


# Determine the ‘dimensionality’ of the dataset
ElbowPlot(pbmc)


## Clustering
# Skip for expeanded version after UMAP.

##Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")


## Challenge: PC genes

# To show people gene expression
FeaturePlot(pbmc, "LYZ")
# how to check if gene is in data
"LYZ" %in% rownames(pbmc)
"Lyz" %in% rownames(pbmc)



## Save
saveRDS(pbmc, file = "pbmc_tutorial_saved.rds") 



## Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)



## Clustring II

# Do clusteringat 0.1, 0.2, 0.3... 2.0
pbmc <- FindClusters(object = pbmc, reduction = "umap", resolution = seq(0.1, 2, 0.1), dims = 1:10)

names(pbmc@meta.data)
# How many clusters (and how many cells in those clusters) do we get at different resolutions?
table(pbmc$RNA_snn_res.0.1)
#> 
#>    0    1    2    3 
#> 1190  688  416  344
table(pbmc$RNA_snn_res.0.5)
#> 
#>   0   1   2   3   4   5   6   7   8 
#> 684 481 476 344 291 162 155  32  13
table(pbmc$RNA_snn_res.2)



library(clustree)
clustree(pbmc, prefix = "RNA_snn_res.") + theme(legend.key.size = unit(0.05, "cm"))


# The name of the cluster is prefixed with 'RNA_snn_res' and the number of the resolution
Idents(pbmc) <- pbmc$RNA_snn_res.0.5

DimPlot(pbmc, label = TRUE, repel = TRUE, label.box = TRUE) + NoLegend()



## Cluster Markers



# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)





# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()



##  5.2 Use markers to label or find a cluster 

genes_markers <- list(Naive_CD4_T = c("IL7R", "CCR7"))
pbmc <- AddModuleScore(object = pbmc, features = genes_markers, ctrl = 5, name = "Naive_CD4_T",
                       search = TRUE)

# label that cell type
pbmc$cell_label = NA
pbmc$cell_label[pbmc$Naive_CD4_T1 > 1] = "Naive_CD4_T"
Idents(pbmc) = pbmc$cell_label

# plot
# Using a custom colour scale 
FeaturePlot(pbmc, features = "Naive_CD4_T1", label = TRUE, repel = TRUE, ) + scale_colour_gradientn(colours = c("lightblue","beige","red"))
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the
#> existing scale.


##  5.3 Assigning cell type identity to clusters 

# Sometimes you need to change Idents, so make sure your favourite cluster is stored in its own column too!
# 5 => c5
pbmc$cluster <- factor( paste0("c", pbmc$RNA_snn_res.0.5),   levels=paste0('c', levels(pbmc$RNA_snn_res.0.5)))
Idents(pbmc) <- pbmc$cluster
levels(pbmc$cluster)
#> [1] "c0" "c1" "c2" "c3" "c4" "c5" "c6" "c7" "c8"


cluster_content <- list(
  c0 = "Naive CD4+ T",
  c1 = "CD14+ Mono",
  c2 = "Memory CD4+",  
  c3 = "B",
  c4 = "CD8+ T",
  c5 = "Mono",
  c6 = "NK",
  c7 = "DC",
  c8 = "Platelet"
)

# "c5" => "Mono" 
pbmc$celltype <- as.character(cluster_content[pbmc$cluster])

# "c5" => "Mono" 
pbmc$celltype <-factor(as.character(cluster_content[pbmc$cluster]), levels=cluster_content)

# c5 => c5: Mono
pbmc$pretty_cluster_labels <- factor (
  paste0(names(cluster_content[pbmc$cluster]), ": ", cluster_content[pbmc$cluster]) , 
  levels = paste0( names(cluster_content), ": ", cluster_content)
)

pbmc@meta.data %>% 
  group_by(RNA_snn_res.0.5, cluster, celltype, pretty_cluster_labels) %>%
  summarise(num_cells=n(), .groups="drop")  %>%
  DT::datatable()



########

## Single R 

library(SingleCellExperiment)
library(SingleR)
library(celldex)



sce <- as.SingleCellExperiment(pbmc)
ls('package:celldex')
ref.set <- celldex::HumanPrimaryCellAtlasData()



unique(ref.set$label.main)
pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.main)


lbls.keep <- table(pred.cnts$labels)>10
pbmc$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')
DimPlot(pbmc, reduction='umap', group.by='SingleR.labels')




################
# Differential Expression

kang <- readRDS("data/kang2018.rds")
head(kang@meta.data)
total_per_gene <- rowSums(GetAssayData(kang, 'RNA', layer = "counts")) #Make sure its RNA assay, layer = counts
hist(log10(total_per_gene))
kang <- kang[total_per_gene >= 50, ] 
Idents(kang) <- kang$cell

kang.celltype <- kang[, kang$cell == "CD14+ Monocytes" ]
DimPlot(kang.celltype)


## Wilcox

# Change Ident to Condition
Idents(kang.celltype) <- kang.celltype$stim

# default, wilcox test
de_result_wilcox <- FindMarkers(kang.celltype, 
                                ident.1 = 'stim',
                                ident.2 = 'ctrl',
                                logfc.threshold = 0, # Give me ALL results
                                min.pct = 0
)

# Add average expression for plotting

de_result_wilcox$AveExpr<- rowMeans(GetAssayData(kang.celltype,assay="RNA", layer = "data")[rownames(de_result_wilcox),])

p1 <- ggplot(de_result_wilcox, aes(x=AveExpr, y=avg_log2FC, col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) + 
  theme_bw() +
  ggtitle("Wilcox Test")


p2 <- ggplot(de_result_wilcox, aes(x=avg_log2FC, y=-log10(p_val), col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) + 
  theme_bw() +
  ggtitle("Wilcox Test (Volcano)")

p1 + p2




## Neg bionm



# Change Ident to Condition
Idents(kang.celltype) <- kang.celltype$stim

# default, wilcox test
de_result_negbinom <- FindMarkers(kang.celltype, 
                                  test.use="negbinom", # Choose a different test.
                                  ident.1 = 'stim',
                                  ident.2 = 'ctrl',
                                  logfc.threshold = 0, # Give me ALL results
                                  min.pct = 0
)

# Add average expression for plotting
de_result_negbinom$AveExpr<- rowMeans(GetAssayData(kang.celltype,assay="RNA", layer = "data")[rownames(de_result_negbinom),])
p1 <- ggplot(de_result_negbinom, aes(x=AveExpr, y=avg_log2FC, col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) + 
  theme_bw() +
  ggtitle("Negative Bionomial Test")


p2 <- ggplot(de_result_negbinom, aes(x=avg_log2FC, y=-log10(p_val), col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) + 
  theme_bw() +
  ggtitle("Negative Bionomial Test (Volcano)")

p1 + p2




## Pseudobulk



# Tools for bulk differential expression
library(limma)
library(edgeR)


# Change idents to ind for grouping.
kang.celltype$sample <- factor(paste(kang.celltype$stim, kang.celltype$ind, sep="_"))
Idents(kang.celltype) <- kang.celltype$sample

# THen pool together counts in those groups
# AggregateExperssion returns a list of matricies - one for each assay requested (even just requesting one)
pseudobulk_matrix_list <- AggregateExpression( kang.celltype,  assays='RNA')
pseudobulk_matrix      <- pseudobulk_matrix_list[['RNA']]
colnames(pseudobulk_matrix) <- as.character(colnames(pseudobulk_matrix)) # Changes colnames to simple text
pseudobulk_matrix[1:5,1:4]
#> 5 x 4 sparse Matrix of class "dgCMatrix"
#>          ctrl-101 ctrl-1015 ctrl-1016 ctrl-1039
#> NOC2L           2         7         .         .
#> HES4            .         3         2         1
#> ISG15          31       185       236        41
#> TNFRSF18        .         3         4         2
#> TNFRSF4         .         2         .         .

# Aggregate expression replaces _ with -! We're going to change it back (for limma.)
colnames(pseudobulk_matrix) <- gsub("-", "_",colnames(pseudobulk_matrix))



dge <- DGEList(pseudobulk_matrix)
dge <- calcNormFactors(dge)

# Remove _ and everything after it - yeilds stim group
stim <- gsub("_.*","",colnames(pseudobulk_matrix)) 

# Removing everything before the _ for the individua, then converting those numerical ind explictiy to text. Else limma will treat them as numbers!
ind  <- as.character(gsub(".*_","",colnames(pseudobulk_matrix))) 

design <- model.matrix( ~0 + stim + ind)
vm  <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(vm, design = design)

contrasts <- makeContrasts(stimstim - stimctrl, levels=coef(fit))
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)

de_result_pseudobulk <- topTable(fit, n = Inf, adjust.method = "BH")
de_result_pseudobulk <- arrange(de_result_pseudobulk , adj.P.Val)


p1 <- ggplot(de_result_pseudobulk, aes(x=AveExpr, y=logFC, col=adj.P.Val < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) + 
  theme_bw() +
  ggtitle("Pseudobulk")


p2 <- ggplot(de_result_pseudobulk, aes(x=logFC, y=-log10(P.Value), col=adj.P.Val < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) + 
  theme_bw() +
  ggtitle("Pseudobulk Test (Volcano)")

p1 + p2



VlnPlot(kang, "ISG20", group.by = 'cell')
VlnPlot(kang.celltype, "ISG20", group.by = 'stim')

FeaturePlot(kang, "ISG20", split.by = 'stim')




## Cell cycle

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Use those lists with the cell cycle scoring function in Seurat.
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(pbmc, reduction = 'umap', group.by = "Phase")


# Integration with harmony

kang <- readRDS("data/kang2018.rds")

library(harmony)

kang <- RunHarmony(kang, c("stim", "ind"), reduction.use="pca")

DimPlot(kang, reduction="pca", group.by="stim")
DimPlot(kang, reduction="harmony", group.by="stim")

kang <- RunUMAP(kang, reduction="harmony", dims=1:10, reduction.name="umap_harmony")
DimPlot(kang, reduction="umap_harmony", group.by="stim")

