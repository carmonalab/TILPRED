#' Predict state of tumor-infiltrating CD8 T-cells from single-cell signature enrichment scores
#'
#' \code{predictTilState} This function evaluates a logistic regression model to predict the state of individual CD8 tumor-infiltrating lymphocytes (mouse or human) based on their transcriptomes (scRNA-seq data)
#'
#' @param data single-cell expression object of CD8 T-cells of class \emph{SingleCellExperiment}. For gene expression matrix, only gene expression ranks in each cell will be used and therefore any cell-to-cell normalization method used equivalent (as long as gene ranks are conserved among the top ~5"\%" genes). E.g. UMI counts, CPM, TPM, TMM are equivalent input types. Gene names correspond to mouse gene symbols (e.g. \emph{Cncb2}).
#'
#' @param human logical value indicating if input matrix correspond to human genes (by default mouse data is expected)
#'
#' @param scoreThreshold probability threshold [0,1] for assigning cell states. If all state probabilities are below this threshold, 'unknown' state is assigned. Default 0.5
#'
#' @param cellCycleThreshold probability threshold [0,1] for assigning cell cyling state; Values in the range [0.1,0.2] are recommended (Default 0.2)
#'
#' @param nCores number of cores used for (AUCell) AUC scores computation (Default 1)
#'
#' @return a two-element list containing 1) \emph{predictedState}, the predicted states
#' (naive, effector memory, exhausted, memoryLike, or "unknown" if no class had a score above a threshold of \emph{scoreThreshold}),
#' 2) \emph{stateProbabilityMatrix}, a matrix of number_of_cells x number_of_states (4) of probabilities of cell c belonging to class s,
#' and 3) \emph{cycling}, logical vector indicating for each cell whethere there is a high cell cycle signal (independent to the cellular sub-type/state signal)
#'
#'
#' @examples
##' data(B16CD8TIL_SCE)
#' x <- predictTilState(B16CD8TIL_SCE)
#' table(x$predictedState)
#' head(x$stateProbabilityMatrix)
#' @export
#'


predictTilState <- function(data,nCores=1,human=F,scoreThreshold=0.5,cellCycleThreshold=0.2) {

  if(scoreThreshold < 0 | scoreThreshold > 1) stop("scoreThreshold is invalid (0 to 1)")
  if(cellCycleThreshold < 0 | cellCycleThreshold > 1) stop("cellCycleThreshold is invalid (0 to 1)")
  if(class(data)!="SingleCellExperiment") stop("input is not a SingleCellExperiment object")


  if(human) {
    #sigs <- lapply(sigs,function(x) unique(orthologMap[names(orthologMap) %in% x]))

    sigGenes=unique(unlist(cellTypeHumanSigs))

    if(sum(!sigGenes %in% rownames(data)) > 0) {
      warning(paste("The following genes were not found in the dataset provided ",paste(sigGenes[!sigGenes %in% rownames(data)],collapse=","),". Doesn't look too bad but prediction performance might be affected."))
      if(mean(!sigGenes %in% rownames(data)) > 0.1) stop("Too many genes not found")
    }

    aucs <- scAUCscore(data,nCores=nCores, sigs=cellTypeHumanSigs)
    rownames(aucs) <- paste0("AUC_",rownames(aucs))

    myCellType.AUC <- rep("nonTcell",ncol(aucs))
    myCellType.AUC[aucs["AUC_Tcell",]>0.1] <- "pureTcell"
    myCellType.AUC[aucs["AUC_Tcell",]>0.1 & (aucs["AUC_B.cell",] > 0.15 | aucs["AUC_CAF",] > 0.05 | aucs["AUC_Endo.",] > 0.1 |aucs["AUC_Macrophage",] > 0.15 | aucs["AUC_Mal",] > 0.15) ] <- "TcellDoublet"
    myCellType.AUC[aucs["AUC_NK",] > 0.2] <- "NK"

    #aucCells <- aucs[,grep("AUC_",colnames(aucs))]
    aucsMax <- (apply(aucs,2,function(x) x[which.max(x)]))
    aucsLabel <- (sub("AUC_","",apply(aucs,2,function(x) rownames(aucs)[which.max(x)])))
    myCellType.AUC[aucs["AUC_Tcell",] <= 0.1 & aucsMax > 0.1] <- aucsLabel[aucs["AUC_Tcell",] <= 0.1 & aucsMax > 0.1]
    return(list(predictedState=myCellType.AUC,stateProbabilityMatrix=NA,cycling=NA))



  } else {


  sigGenes=unique(unlist(sigs))

    if(class(data)!="SingleCellExperiment") stop("input is not a SingleCellExperiment object")

    if(sum(!sigGenes %in% rownames(data)) > 0) {
      warning(paste("The following genes were not found in the dataset provided ",paste(sigGenes[!sigGenes %in% rownames(data)],collapse=","),". Doesn't look too bad but prediction performance might be affected."))
      if(mean(!sigGenes %in% rownames(data)) > 0.1) stop("Too many genes not found")
    }

    aucs <- scAUCscore(data,nCores=nCores, sigs=sigs)
    rownames(aucs)[1:12] <- paste0("AUC_PW_DEG_",rownames(aucs)[1:12])


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
  return(list(predictedState=cellClass,stateProbabilityMatrix=probM,cycling=aucs["cycling",]>cellCycleThreshold))
  }

}



