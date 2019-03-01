#' Predict state of tumor-infiltrating CD8 T-cells from single-cell signature enrichment scores
#'
#' \code{predictTilState} This function evaluates a logistic regression model to predict the state of individual CD8 tumor-infiltrating lymphocytes (mouse or human) based on their transcriptomes (scRNA-seq data)
#'
#' @param aucs matrix of enrichment scores (scAUCscore)
#'
#' @return a two-element list containing 1) \emph{predictedState}, the predicted states
#' (naive, terminal effector, exhausted, memory-like, cycling effector or "unknown" if no class had a score above a threshold of 0.5), and
#' 2) \emph{stateProbabilityMatrix}, a matrix of number_of_cells x number_of_states (5) of probabilities of cell c belonging to class s
#'
#' @examples
#' data(B16CD8TILs_tpm)
#' x <- predictTilState(B16CD8TILs_tpm)
#' table(x$predictedState)
#' head(x$stateProbabilityMatrix)
#' @export
#'


predictTilState <- function(aucs) {

  odf = function(x) exp(colSums((aucs[names(x),] * x)))

  scoreM=matrix(NA,nrow=ncol(aucs),ncol=length(glm_coefs))
  for (i in seq_along(glm_coefs)){
    scoreM[,i]=odf(glm_coefs[[i]])
  }

  probM=scoreM/rowSums(scoreM)
  colnames(probM)=classNames
  rownames(probM)=colnames(aucs)
  probM.un=probM
  probM.un[probM<0.5]=NA
  cellClass=as.character(apply(probM.un,1,function(x) classNames[which.max(x)]))
  cellClass[cellClass=="character(0)"]="unknown"
  cellClass=factor(cellClass,levels=classNames)
  return(list(predictedState=cellClass,stateProbabilityMatrix=probM,cycling=aucs["cycling",]>0.2))

}



