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
#' @param filterCD8T automatic pre-filter CD8 T cells before classifing CD8 T cell states (Default TRUE). If filterCD8T is set to FALSE, all input cells are assumed to be CD8 T cells, even though they migh show no features of this cell type.
#'
#'
#' @return a two-element list containing 1) \emph{predictedState}, the predicted states for CD8 T cells
#' (naive, effector memory, exhausted, memoryLike, or "unknown" if no class had a score above a threshold of \emph{scoreThreshold}); or the predicted cell type for non CD8 T cells: Treg (Foxp3 Regulatory T cells), CD4T (non Treg CD4+ T cells), NKT (NK T cells), Tcell_unknown (T cells of other kinds) and Non-Tcell (for cell types other than T cells, e.g. Myeloid, B cells, NKs);
#' 2) \emph{stateProbabilityMatrix}, a matrix of number_of_cells x number_of_states (4) of probabilities of cell c belonging to class s, only for CD8 T cells;
#'  3) \emph{cycling}, logical vector indicating for each cell whethere there is a high cell cycle signal (independent to the cellular sub-type/state signal), and 4) \emph{cyclingScore} AUC score for the cell cycle signature
#'
#'
#' @examples
#' data(B16CD8TIL_SCE)
#' x <- predictTilState(data=B16CD8TIL_SCE)
#' table(x$predictedState)
#' head(x$stateProbabilityMatrix)
#' head(x$cyclingScore)
#' @export
#'


predictTilState <- function(data,nCores=1,human=F,scoreThreshold=0.5,cellCycleThreshold=0.2,filterCD8T=T) {

  if(scoreThreshold < 0 | scoreThreshold > 1) stop("scoreThreshold is invalid (0 to 1)")
  if(cellCycleThreshold < 0 | cellCycleThreshold > 1) stop("cellCycleThreshold is invalid (0 to 1)")
  if(class(data)!="SingleCellExperiment") stop("input is not a SingleCellExperiment object")


  if(human) {

    sigGenes=unique(unlist(cellTypeHumanSigs))

    if(sum(!sigGenes %in% rownames(data)) > 0) {
      warning(paste("The following genes were not found in the dataset provided ",paste(sigGenes[!sigGenes %in% rownames(data)],collapse=","),". Doesn't look too bad but prediction performance might be affected."))
      if(mean(!sigGenes %in% rownames(data)) > 0.1) stop("Too many genes not found")
    }

    aucs <- scAUCscore(data,nCores=nCores, sigs=cellTypeHumanSigs, aucMaxRank=1500)
    rownames(aucs) <- paste0("AUC_",rownames(aucs))

    celltype.pred <- rep("nonTcell",ncol(aucs))
    celltype.pred[aucs["AUC_Tcell",]>0.1] <- "pureTcell"
    celltype.pred[aucs["AUC_Tcell",]>0.1 & (aucs["AUC_B.cell",] > 0.15 | aucs["AUC_CAF",] > 0.05 | aucs["AUC_Endo.",] > 0.1 |aucs["AUC_Macrophage",] > 0.15 | aucs["AUC_Mal",] > 0.15) ] <- "TcellDoublet"
    celltype.pred[aucs["AUC_NK",] > 0.2] <- "NK"

    aucsMax <- (apply(aucs,2,function(x) x[which.max(x)]))
    aucsLabel <- (sub("AUC_","",apply(aucs,2,function(x) rownames(aucs)[which.max(x)])))
    celltype.pred[aucs["AUC_Tcell",] <= 0.1 & aucsMax > 0.1] <- aucsLabel[aucs["AUC_Tcell",] <= 0.1 & aucsMax > 0.1]
    return(list(predictedState=celltype.pred,stateProbabilityMatrix=NA,cycling=NA))



  } else {


    if(mean(sigs$Tcell %in% rownames(data)) < 0.7) {
      stop("T cell genes not found")
    }

    sigGenes=unique(unlist(sigs))

    if(sum(!sigGenes %in% rownames(data)) > 10) {
      warning(paste("The following genes were not found in the dataset provided ",paste(sigGenes[!sigGenes %in% rownames(data)],collapse=","),". Doesn't look too bad but prediction performance might be affected."))
      if(mean(!sigGenes %in% rownames(data)) > 0.1) stop("Too many genes not found")
    }

    aucs <- scAUCscore(data,nCores=nCores, sigs=sigs)

    celltype.pred <- rep("unknown",ncol(aucs))
    names(celltype.pred) <- colnames(aucs)
    celltype.pred[aucs["Tcell",]>0.1] <- "Tcell_unknown"
    celltype.pred[which(aucs["Tcell",]>0.1 & aucs["CD8T",]>0.1 & (aucs["CD8T",]>aucs["Cd4",])) ] <- "CD8T"
    celltype.pred[aucs["Tcell",]>0.1 & aucs["Cd4",]>0.1 & (aucs["CD8T",]<aucs["Cd4",]) ] <- "CD4T"
    celltype.pred[aucs["Tcell",]>0.1 & aucs["Treg",]>0.1 ] <- "Treg"
    celltype.pred[aucs["Tcell",]>0.1 & aucs["NK",]>0.3] <- "NKT"
    celltype.pred[apply(aucs[c("Myel","B","NK"),],2,max)>0.2] <- "Non-Tcell"


    if(filterCD8T){
      print(table(celltype.pred, exclude = NULL))
      aucsCD8T <- aucs[,celltype.pred == "CD8T"]
    } else {
      aucsCD8T <- aucs
    }

    rownames(aucsCD8T)[1:12] <- paste0("AUC_PW_DEG_",rownames(aucsCD8T)[1:12])


    odf = function(x) exp(colSums((aucsCD8T[names(x),] * x)))

    scoreM=matrix(NA,nrow=ncol(aucsCD8T),ncol=length(glm_coefs))
    for (i in seq_along(glm_coefs)){
      scoreM[,i]=odf(glm_coefs[[i]])
    }

    probM=scoreM/rowSums(scoreM)
    colnames(probM)=paste0("CD8T_",classNames)
    rownames(probM)=colnames(aucsCD8T)
    probM.un=probM
    probM.un[probM < scoreThreshold]=NA
    CD8Tstate=as.character(apply(probM.un,1,function(x) classNames[which.max(x)]))
    CD8Tstate[!CD8Tstate %in% classNames]="unknown"
    CD8Tstate=factor(CD8Tstate,levels=c(classNames,"unknown"))
    names(CD8Tstate)=rownames(probM)
    print(table(CD8Tstate))

    celltype.pred[names(CD8Tstate)] <- paste0("CD8T_",as.character(CD8Tstate))

    return(list(predictedState=celltype.pred,stateProbabilityMatrix=probM,cycling=aucs["cycling2",]>cellCycleThreshold,cyclingScore=aucs["cycling2",]))

  }

}



