---
title: "metaTF: An Ensemble R package for Systematic Analysis of Transcription Factor Regulon at Single-cell level"
output: html_notebook
---
Transcription factor regulon (TFR) 

## 1. Overview
The **metaTF** is an R package used for <u>**s**</u>ingle-<u>**c**</u>ell <u>**A**</u>ctivate <u>**T**</u>ranscription <u>**F**</u>actor <u>**R**</u>egulon analysis, which combining the transcription factor regulatory information with single-cell gene expression data and evaluating the transcription factor regulon activities, identifying cell-type specific regulons and assessing the similarities between regulons and pathways. 
<u>**Conda install***</u> 

```{r}
conda create -n metaTF_env
conda activate metaTF_env
conda install -c conda-forge r-seurat -y
conda install -c conda-forge r-ppcor r-cvtools r-factoextra r-glmnet r-ksamples r-progress r-arrow r-systemfonts r-textshaping r-cairo r-ragg r-devtools -y
conda install -c bioconda bioconductor-summarizedexperiment -y
R
```
<u>**Run in R***</u> 

```{r}
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("GENIE3")
BiocManager::install("fgsea")
BiocManager::install("RcisTarget")
BiocManager::install("qvalue")
BiocManager::install("scater")
BiocManager::install("minet")
BiocManager::install("viper")
install.packages("Radviz")
install.packages("jaccard")
install.packages("scLink")
install.packages("devtools")
library(devtools)
install_github("Wanglabsmu/metatf")
```
<u>**Docker install***</u> 
<u>**A**</u>s a recommendation, we recommend using docker to install the environment in which metaTF(previous name scATFR) runs

```{r}
docker pull hobartjoe/scatfr
conda activate scATFR_env
```
## 2. Quick Start
Firstly, we load TF-target network from local data set and import expression data: 
```{r}
library(metaTF)
#load TF-target data 
data(dorothea_regulons)
regulons <- dfToList(df = dorothea_regulons$mm_regulons, tf_col = "tf",target_col = "target")
#load gene expression data and colData
expression_data <- readRDS(system.file("extdata", "mouse_HSC_formation_expression.rds", package = "metaTF"))
col_data <- readRDS(system.file("extdata", "mouse_HSC_formation_colData.rds", package = "metaTF"))
```
Then, we create the `atfr` object:
```{r}
atfr <- metaTF(exp_data = expression_data, col_data = col_data, regulons = regulons)
```
Next, we infer the Gene Regulatory Networks:
```{r}
atfr <- inferGRNs(x = atfr, method="PUIC", ncores=6)
```
To look for GRNs
```{r}
head(GRN(atfr)[,1:5])
```
After that, we filter the gene regulons using gene regulotory network:
```{r}
atfr <- filterRegulons(x=atfr, ncores=6)
plotRegulonReducedDim(object = atfr, alt_assay = "viper", dimred = "UMAP", colour_by="stage")
```
Now, we evaluate the regulon activity in each single cell:
```{r}
atfr <- regulonActivity(atfr, method="viper",ncores=6)
atfr <- regulonUMAP(x = atfr)
```
We can also get the cell-type specific regulons:
```{r}
#get one cell type specific
EC_sig <- findCellTypeSignature(atfr, pos_clust="EC",clust_name="stage")
#get all cell type specific
all_sig <- findAllCTSRegulons(atfr, clust_name="stage",ncores=6)
head(all_sig)
all_sig <- dplyr::tibble(all_sig)
```
We can compare the similarities between regulons and pathways:
```{r}
data(Reactome_pathway)
Reactome_mouse <- Reactome_pathway$mouse
regulon_pathway <- inferRegulonFun(x = atfr, pathway_list = Reactome_mouse, ncores=6)
```

### 1.2 Obtaining TFRs from DoRoThEA
Alternatively, you can also simply obtain these data from `DoRoThEA` package. The TFRs for human and mouse were deposited in latest version (1.3.2) of [dorothea](https://github.com/saezlab/dorothea) package, and can be directly extracted. Please noticed that for the standard analysis using normal tissues, `dorothea` suggest using regulons data in `dorothea_hs/dorothea_mm`, in which the tissue specific gene expression data from the GTEx portal. For abnormal tissues, the `dorothea_hs_pancancer/dorothea_mm_pancancer` is better, in which the gene expression data were  inferred from the TCGA program.

**NOTE**: only the latest version (1.3.2) include the `dorothea_hs_pancancer/dorothea_mm_pancancer` data.
```{r}
#install.packages("devtools")
#devtools::install_github("saezlab/dorothea")
library(dorothea)
mouse_regulons <- get(data("dorothea_mm", package = "dorothea"))
head(mouse_regulons)
human_abnorm_regulons <- get(data("dorothea_hs_pancancer", package = "dorothea"))
head(human_abnorm_regulons)
tfrs <- tfrFromDoRoThEA(mouse_regulons, levels=c("A","B","C"))
```
Typically, there are four columns in TFRs data:

**tf**: transcription factors in gene symbols from HGNC or MGI.

**confidence**: the confidence between TFs and their Targets. Level `A` stands for highest confidence while level `E` stands for lowest confidence.

**target**: target genes in gene symbols from HGNC or MGI.

**mor**: the mode of regulation. `1` stands for activation and `-1` stands for repression.

## Case Study

### 2. CASE1: Characterisation of cell-type specific transcription factor regulons (TFRs) during HSC formation.
In this case, we firstly performed dimension reduction and pseudotime analysis based on [Monocle3](https://cole-trapnell-lab.github.io/monocle3/) using mouse embryonic hematopoietic stem cells (HSCs) formation data ([Fan zhou, 2016](https://www.nature.com/articles/nature17997) and [Jie Zhou, 2018](https://doi.org/10.1016/j.stem.2018.11.023)). Then, we performed transcription factor regulons (TFRs) analysis such as cell-type specific TFRs identification and TFR-pathway similarity analysis using **metaTF**. Finally, we showed the visualization of the above results.

#### 2.1 Pseudotime analysis using monocle3

#### 2.1.1 Loading the data and create `cds` object
The row names of `expression_matrix` is gene symbols and the column names of `expression_matrix` is cell identifier. The `expression_colData` include three columns: sample, stage and tissue.
```{r}
library(monocle3)
#load the expression data
expression_data <- readRDS(system.file("extdata", "mouse_HSC_formation_expression.rds", package = "metaTF"))
#load the cell colData
expression_colData <- readRDS(system.file("extdata", "mouse_HSC_formation_colData.rds", package = "metaTF"))
#create cds object
cds <- new_cell_data_set(expression_data = expression_data,
                        cell_metadata = expression_colData)
```
#### 2.1.2 pre-process and dimension reduction
Then, we perform normalization and dimension reduction with object `cds`.
```{r}
#Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50, method = "PCA")
#Reduce the dimensions using UMAP
cds <- reduce_dimension(cds,reduction_method = "UMAP", preprocess_method = "PCA")
#plot_cells(cds, reduction_method = "UMAP", color_cells_by="stage")
```
#### 2.1.3 Pseudotime analysis
```{r}
# Cluster the cells
cds <- cluster_cells(cds)
# Learn a graph
cds <- learn_graph(cds,use_partition = FALSE)
# Order cells
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "stage",show_trajectory_graph = TRUE)
```
Finally, to better perform regulon analysis, we need to add pseudotime and cell clusters information into colData of cds object, this is convenient for further analysis.
```{r}
pseudotime <- pseudotime(cds)
clusters <- clusters(cds)
colData(cds)$clusters <- clusters[row.names(colData(cds))]
colData(cds)$pseudotime <- pseudotime[row.names(colData(cds))]
head(colData(cds))
```
Now, we can see pseudotime information in colData of `cds`.

#### 2.2   

```{r}
data(dorothea_regulons)
mouse_regulons <- dorothea_regulons$mm_regulons
tfrs <- tfrFromDoRoThEA(mouse_regulons, levels=c("A","B","C","D","E"))
cds <- regulonActivity(cds, gene_list = tfrs, method = "viper")
tfr_mat <- altExp(cds,"atfr_raw")
tfr_mat <- runAUCell(assay(cds),gene_list = tfrs)

```









