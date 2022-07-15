#' Calculate CD8 TIL states enrichment scores (required to run predictTilState function)
#'
#' \code{scAUCscore} This function evaluates a logistic regression model to predict the state of individual CD8 tumor-infiltrating lymphocytes (mouse or human) based on their transcriptomes (scRNA-seq data)
#'
#' @param data single-cell expression matrix. Only gene expression ranks in each cell will be used and therefore any cell-to-cell normalization method used is not relevant (e.g. UMI counts, CPM, TPM, TMM)
#' @param nCores number of cores to use in parallel calculation of gene rankings
#' @export



scAUCscore <- function(sce, nCores=1, sigs,  aucMaxRank=1500) {

  set.seed(123)
  cells_rankings <- AUCell::AUCell_buildRankings(SingleCellExperiment::logcounts(sce),
        nCores=nCores, plotStats=F,verbose = F, splitByBlocks=TRUE)
  cells_AUC <- AUCell::AUCell_calcAUC(sigs, cells_rankings, aucMaxRank=aucMaxRank)
  aucs <- AUCell::getAUC(cells_AUC)
  return(aucs)
}



