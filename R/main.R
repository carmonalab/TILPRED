#install.packages("BiocManager")
#install.packages("doParallel")
#install.packages("doRNG")
#BiocManager::install("GenomeInfoDbData",type="source")
#BiocManager::install("AUCell")
#BiocManager::install("SingleCellExperiment",type="source")
library(AUCell)
library(SingleCellExperiment)
library(doParallel)
library(doRNG)


glm_coefs <- readRDS("data/coefs.list.rds")
sigs <- readRDS("data/cellTypeSignatures_mouse.Rds")
cellTypeHumanSigs <- readRDS("data/cellTypeSignatures_human.Rds")
orthologMap <- readRDS("data/mapHsa_vs_Mmu_Orthologs.rds")



predictorSigs=unique(unlist(lapply(glm_coefs,names)))
classNames=names(glm_coefs)
