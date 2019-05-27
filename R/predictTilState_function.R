#' Predict state of tumor-infiltrating CD8 T-cells from single-cell signature enrichment scores
#'
#' \code{predictTilState} This function evaluates a logistic regression model to predict the state of individual CD8 tumor-infiltrating lymphocytes (mouse or human) based on their transcriptomes (scRNA-seq data)
#'
#' @param data single-cell expression object of CD8 T-cells of class \emph{SingleCellExperiment} (if is.auc is FALSE) or matrix of enrichment scores (scAUCscore output, if is.auc is TRUE). For gene expression matrix, only gene expression ranks in each cell will be used and therefore any cell-to-cell normalization method used equivalent (as long as gene ranks are conserved among the top ~5"\%" genes). E.g. UMI counts, CPM, TPM, TMM are equivalent input types. Gene names correspond to mouse gene symbols (e.g. \emph{Cncb2}).
#'
#' @param is.auc logical value indicating if input matrix \emph{data} corresponds to single-cell expression matrix (FALSE, default) or matrix of enrichment scores, as calculated by scAUCscore function (TRUE)
#'
#' @param human logical value indicating if input matrix correspond to human genes (by default mouse data is expected)
#'
#' @param scoreThreshold probability threshold (0 to 1) for assigning cell states. If all state probabilities are below this threshold, 'unknown' state is assigned. Default is 0.5
#'
#' @return a two-element list containing 1) \emph{predictedState}, the predicted states
#' (naive, effector, exhausted, memoryLike, or "unknown" if no class had a score above a threshold of \emph{scoreThreshold}),
#' 2) \emph{stateProbabilityMatrix}, a matrix of number_of_cells x number_of_states (4) of probabilities of cell c belonging to class s,
#' 3) \emph{cycling}, logical vector indicating for each cell whethere there is a high cell cycle signal (independent to the cellular sub-type/state signal), and 4) \emph{scAUCscore} the matrix of AUC signature enrichment scores internally used for the prediction model (it can be used to re-run the prediction model)
#'
#'
#' @examples
##' data(B16CD8TIL_SCE)
#' x <- predictTilState(B16CD8TIL_SCE)
#' table(x$predictedState)
#' head(x$stateProbabilityMatrix)
#' @export
#'


predictTilState <- function(data,nCores=1,is.auc=F,human=F,scoreThreshold=0.5) {

  if(human) {
    sigs <- lapply(sigs,function(x) unique(orthologMap[names(orthologMap) %in% x]))
  }

  if(scoreThreshold < 0 | scoreThreshold > 1) stop("scoreThreshold is invalid (0 to 1)")

  sigGenes=unique(unlist(sigs))

  if(!is.auc){

    if(class(data)!="SingleCellExperiment") stop("input is not a SingleCellExperiment object")

    if(sum(!sigGenes %in% rownames(data)) > 0) {
      warning(paste("The following genes were not found ",paste(sigGenes[!sigGenes %in% rownames(data)],collapse=","),". Unknown prediction performance."))
      if(mean(!sigGenes %in% rownames(data)) > 0.1) stop("Too many genes not found")
    }

    aucs <- scAUCscore(data,nCores=nCores, sigs=sigs)

  } else {
    aucs <- data
  }

  odf = function(x) exp(colSums((aucs[names(x),] * x)))

  scoreM=matrix(NA,nrow=ncol(aucs),ncol=length(glm_coefs))
  for (i in seq_along(glm_coefs)){
    scoreM[,i]=odf(glm_coefs[[i]])
  }

  probM=scoreM/rowSums(scoreM)
  colnames(probM)=classNames
  rownames(probM)=colnames(aucs)
  probM.un=probM
  probM.un[probM < scoreThreshold]=NA
  cellClass=as.character(apply(probM.un,1,function(x) classNames[which.max(x)]))
  cellClass[!cellClass %in% classNames]="unknown"
  cellClass=factor(cellClass,levels=c(classNames,"unknown"))
  names(cellClass)=rownames(probM)
  return(list(predictedState=cellClass,stateProbabilityMatrix=probM,cycling=aucs["cycling",]>0.2,scAUCscore=aucs))

}



