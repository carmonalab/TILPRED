# `TILPRED`: Tumor-Infiltrating CD8+ Lymphocytes states Predictor

TILPRED is an R Package for the classification of tumor-infiltrating T lymphocytes (TILs) from single-cell RNA-seq data.

`TILPRED` is a logistic regression-based classifier that reads a [SingleCellExperiment](https://doi.org/doi:10.18129/B9.bioc.SingleCellExperiment) object containing CD8 T cell profiles and assigns to each cell a probability score of belonging to any of the following reference CD8 TIL transcriptomic states: 

* _Exhausted_: a.k.a terminally exhausted cells, associated with terminal differentiation in the context of sustained antigenic stimulation. Phenotypically characterized by co-expression of inhibitory receptors PD-1(_Pdcd1_), TIM3(_Havcr2_), CD39(_Entpd1_), cytotoxic molecules (e.g. _Gzmb_) and lack of Tcf1 (_Tcf7_) expression
* _MemoryLike_: a.k.a progenitor exhausted cells, also associated to sustained antigenic stimulation but retain capacity to self-renew and give rise to _exhausted_ cells. Phenotypically characterized by co-expression of PD-1 (_Pdcd1_) and Tcf1 (_Tcf7_)
* _EffectorMemory_: antigen experienced T cells with effector memory features (e.g. co-expression of cytotoxicity genes such as _Gzmk_ and _Gzmb_, and memory genes such as _Tcf7_, _Lef1_, _Il7r_ and _Ly6c2_). These cells have low to intermediate expression of PD-1, and resemble CD8 T cells found upon acute infection (i.e. in the absence of sustained antigenic stimulation)
* _Naive_: Naive-like CD8 T cells (high expression of _Tcf7_, _Lef1_, _Il7r_, no expression of cytotoxicity genes or T cell activation markers such as CD44, CD69, etc.)

In addition, it predicts proliferation/cycling in each cell, independently of the CD8 TIL subtype.  `TILPRED` uses gene rankings only and therefore is robust to different data normalization strategies. It was tested with scRNA-seq data produced with plate-based (smart-seq2) and droplet-based (10X 5' and 3' counting) technologies. 

Before computing CD8 T cell states probabilities, `TILPRED` will automatically detect non CD8 T cell types. Non CD8 T cells are classified based on curated gene signature enrichment into: _Treg_ (_Foxp3_ Regulatory T cells), _CD4T_ (non Treg CD4+ T cells), _NKT_ (NK T cells), _Tcell_unknown_ (T cells of other kinds) and _Non-Tcell_ (for cell types other than T cells, e.g. Myeloid, B cells, NKs)


Details on the reference CD8 TIL transcriptomic states and TILPRED construction and benchmarking are available in [Carmona SJ et al. 2020](https://doi.org/10.1080/2162402X.2020.1737369)

NB: Currently, TILPRED classifies CD8 TILs from mouse only. Human functionality is under development. TILPRED using parameter human=T will only classify human cell types from tumor samples into: pure T cells, contaminated T cells, NKs and other cell types (NK, macrophages, CAF, epithelial, endothelial, malignant cells, etc.)


### Package Installation

TILPRED requires [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html), [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html) and the Bioconductor packages [AUCell](https://bioconductor.org/packages/release/bioc/html/AUCell.html) and [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)


```
install.packages(c("doParallel","doRNG"))
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("AUCell","SingleCellExperiment"))
library("SingleCellExperiment")
library("AUCell")
```

To install TILPRED directly from the Git Repo use [remotes](https://cran.r-project.org/web/packages/remotes/index.html)

```
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("carmonalab/TILPRED")
library(TILPRED)
```

### Package usage

Run TILPRED on a SingleCellExperiment object containing the single-cell expression matrix of CD8 T cells 
```
data(B16CD8TIL_SCE) # example SingleCellExperiment object
sce.pred <- predictTilState(data=B16CD8TIL_SCE)
```

View output
```
table(sce.pred$predictedState)
head(sce.pred$stateProbabilityMatrix)
```

### For a running example please see this [R Notebook](https://github.com/carmonalab/testTILPRED)


### To cite TILPRED 
Santiago J. Carmona, Imran Siddiqui, Mariia Bilous, Werner Held & David Gfeller (2020) Deciphering the transcriptomic landscape of tumor-infiltrating CD8 lymphocytes in B16 melanoma tumors with single-cell RNA-Seq, OncoImmunology, 9:1, DOI: [10.1080/2162402X.2020.1737369](https://doi.org/10.1080/2162402X.2020.1737369)
