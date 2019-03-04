glm_coefs <- readRDS("data/coefs.list.rds")
sigs <- readRDS("data/CD8statesSignatures.rds")
orthologMap <- readRDS("data/mapHsa_vs_Mmu_Orthologs.rds")


predictorSigs=unique(unlist(lapply(glm_coefs,names)))
classNames=names(glm_coefs)

#my_aucs <- readRDS("data/glm.data.rds")
#my_preds <- readRDS("data/glm.data.pred.rds")
#sce <- readRDS("data/B16CD8TIL_SCE.rds")
