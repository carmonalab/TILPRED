glm_coefs <- readRDS("data/coefs.list.rds")
sigs <- readRDS("data/CD8statesSignatures.rds")

predictorSigs=unique(unlist(lapply(glm_coefs,names)))
classNames=names(glm_coefs)

#my_aucs <- readRDS("data/glm.data.rds")
#my_preds <- readRDS("data/glm.data.pred.rds")
#sce <- readRDS("data/data.sce.rds")
