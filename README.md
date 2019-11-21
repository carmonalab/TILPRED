# TILPRED: Tumor-Infiltrating CD8+ Lymphocytes states Predictor

TILPRED is an R Package for the classification of tumor-infiltrating CD8+ lymphocytes (TILs) from single-cell RNA-seq data.

`TILPRED` is a logistic regression-based classifier that reads a [SingleCellExperiment](https://doi.org/doi:10.18129/B9.bioc.SingleCellExperiment) object containing (any number of) CD8 T cell profiles and assigns to each cell a probability score of belonging to any of the following reference CD8 TIL transcriptomic states: 

* _Exhausted_: a.k.a terminally exhausted cells, associated with terminal differentiation in the context of sustained antigenic stimulation. Phenotypically characterized by co-expression of inhibitory receptors PD-1(_Pdcd1_), TIM3(_Havcr2_), CD39(_Entpd1_), cytotoxic molecules (e.g. Gzmb) and lack of Tcf1 (_Tcf7_) expression
* _MemoryLike_: a.k.a progenitor exhausted cells, also associated to sustained antigenic stimulation but retain capacity to self-renew and give rise to _exhausted_ cells. Phenotypically characterized by co-expression of PD-1 (_Pdcd1_) and Tcf1 (_Tcf7_)
* _EffectorMemory_: antigen experienced T cells with effector memory features (e.g. co-expression of cytotoxicity genes such as Gzmk and Gzmb, and memory genes such as Tcf7, Lef1, Il7r and Ly6c2). These cells have low to intermediate expression of PD-1, and resemble CD8 T cells found upon acute infection (i.e. in the absence of sustained antigenic stimulation)
* _Naive_: Naive-like CD8 T cells (high expression of Tcf7, Lef1, Il7r, no expression of cytotoxicity genes or T cell activation markers such as CD44, CD69, etc.)

In addition, it predicts proliferation/cycling in each cell, independently of the CD8 TIL subtype. It was trained to classify CD8 TILs from mouse, and therefore its use with human CD8 TILs via ortholog mapping is only experimental (`TILPRED` will automatically map relevant orthologous genes between the two species).

`TILPRED` uses gene rankings only and therefore is robust to different data normalization strategies. It was tested with scRNA-seq data produced with plate-based (smart-seq2) and droplet-based (10X 5' and 3' counting) technologies. 

 
Details on the reference CD8 TIL transcriptomic states and TILPRED construction and benchmarking are available in [Carmona SJ et al BioRXiv](https://doi.org/10.1101/800847)


### Package Installation

TILPRED requires [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html), [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html) and the Bioconductor packages [AUCell](https://bioconductor.org/packages/release/bioc/html/AUCell.html) and [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)


```
install.packages(c("doParallel","doRNG"))
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("AUCell")
BiocManager::install("SingleCellExperiment")
```

To install TILPRED directly from the GitHub Repo use [devtools](https://cran.r-project.org/web/packages/devtools/index.html)

```
install.packages("devtools")
library(devtools)
install_github("carmonalab/TILPRED")
```

### Package usage

Run TILPRED on a SingleCellExperiment object containing the single-cell expression matrix of CD8 T cells (Automatic cell type detection coming soon, meanwhile make sure to remove non CD8 T cells from your matrix) 
```
sce.pred <- predictTilState(sce)
```

View output
```
table(sce.pred$predictedState)
head(sce.pred$stateProbabilityMatrix)
```

For a running example please see this [R Notebook](https://github.com/carmonalab/testTILPRED)
