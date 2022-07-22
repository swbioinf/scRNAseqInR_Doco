# (PART) Content {-} 

# Using a seurat object {#seuratobject}


Load an existing object

Plot some clusters on a umap.

Plot some gene expression on a umap

plot some gene expression by cluster.



Look at a table of cell metadata.

Look at counts data / expression data

Does gene exist?


<style type="text/css">
.watch-out {
  background-color: lightpink;
  border: 3px solid red;
  font-weight: bold;
}
</style>

Then we assign a class `watch-out` to the code chunk via the
chunk option `class.source`.


```{.r .watch-out}
mtcars[1:5, "mpg"]
#> [1] 21.0 21.0 22.8 21.4 18.7
```



## Load an existing Seurat object




```r
library(tidyverse)
library(Seurat)
```


```r
pbmc <- readRDS("data/pbmc_tutorial.rds")
```


## Plotting



```r
DimPlot(pbmc, reduction = "umap")
```

<img src="04-01-seuratobject_files/figure-html/unnamed-chunk-4-1.png" width="672" />




```r
FeaturePlot(pbmc, features =  "CD14")
```

<img src="04-01-seuratobject_files/figure-html/unnamed-chunk-5-1.png" width="672" />


```r
VlnPlot(pbmc, features = 'CD14')
```

<img src="04-01-seuratobject_files/figure-html/unnamed-chunk-6-1.png" width="672" />



## How is the data stored?


Some of the most important information for working with Seurat objects is in the metatdata.
This is cell level information.


```r
head(pbmc@meta.data)
#>                  orig.ident nCount_RNA nFeature_RNA
#> AAACATACAACCAC-1     pbmc3k       2419          779
#> AAACATTGAGCTAC-1     pbmc3k       4903         1352
#> AAACATTGATCAGC-1     pbmc3k       3147         1129
#> AAACCGTGCTTCCG-1     pbmc3k       2639          960
#> AAACCGTGTATGCG-1     pbmc3k        980          521
#> AAACGCACTGGTAC-1     pbmc3k       2163          781
#>                  percent.mt RNA_snn_res.0.5 seurat_clusters
#> AAACATACAACCAC-1  3.0177759               0               0
#> AAACATTGAGCTAC-1  3.7935958               3               3
#> AAACATTGATCAGC-1  0.8897363               2               2
#> AAACCGTGCTTCCG-1  1.7430845               5               5
#> AAACCGTGTATGCG-1  1.2244898               6               6
#> AAACGCACTGGTAC-1  1.6643551               2               2
```

That doesn't have any gene expression though, that's stored in 






<!-- Taken from MBP material,  -->

#### Discussion: The Seurat Object in R {.challenge}

Lets take a look at the seurat object we have just created in R, `pbmc`

To accomodate the complexity of data arising from a single cell RNA seq experiment, the seurat object keeps this as a container of multiple data tables that are linked.

The functions in seurat can access parts of the data object for analysis and visualisation, we will cover this later on. 

There are a couple of concepts to discuss here.
<details>
<summary>**Class**</summary>

These are essentially data containers in R as a class, and can accessed as a variable in the R environment. 

Classes are pre-defined and can contain multiple data tables and metadata. For Seurat, there are three types. 

* Seurat - the main data class, contains all the data.
* Assay - found within the Seurat object. Depending on the experiment a cell could have data on RNA, ATAC etc measured
* DimReduc - for PCA and UMAP

</details>

<details>
<summary>**Slots**</summary>

Slots are parts within a class that contain specific data. These can be lists, data tables and vectors and can be accessed with conventional R methods.
 
</details>


<details>
<summary>**Data Access**</summary>

Many of the functions in Seurat operate on the data class and slots within them seamlessly. There maybe occasion to access these separately to `hack` them, however this is an advanced analysis method. 

The ways to access the slots can be through methods for the class (functions) or with standard R accessor nomenclature.
</details>

**Examples of accessing a Seurat object**

The `assays` slot in `pbmc` can be accessed with `pbmc@assays`.

The `RNA` assay can be accessed from this with `pbmc@assays$RNA`. 

We often want to access assays, so Seurat nicely gives us a shortcut `pbmc$RNA`. You may sometimes see an alternative notation `pbmc[["RNA"]]`.

In general, slots that are always in an object are accessed with `@` and things that may be different in different data sets are accessed with `$`.

**Have a go**

Use `str` to look at the structure of the Seurat object `pbmc`.

What is in the `meta.data` slot within your Seurat object currently? What type of data is contained here?

Where is our count data within the Seurat object? 

