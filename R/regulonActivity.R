#' @export
setGeneric("regulonActivity", function(x, ...) standardGeneric("regulonActivity"))
#' Performing regulon activity evaluation
#'
#'
#' @param x An object with gene expression data used for regulon activity evaluation,
#' which can be \code{matrix}, \code{SingleCellExperiment}, \code{Seurat} and \code{ExpressionSet}.
#' @param gene_list List. Regulons list used for activity evaluation, names of each regulons
#' should be Transcription Factors(TFs).
#' @param method Method used to calculate the active score of regulons for each sample.
#' @param normalization Logical, whether to normalize the activate score?
#' @param min_target_num Minimum number of genes in your regulons.
#' @param max_target_num Maximum number of genes in your regulons.
#' @param use_assay For object with \code{SummarizedExperiment} class, the use_assay
#'  argument can be used to select the assay containing the expression data we want as input.
#' @param use_regulon The regulon assay used for activity estimation. Default: 1, the first one.
#' @param norm_act Whether to scale the activity score? Default: FALSE.
#' @param norm_coe The coefficient used for activity score normalization. 
#' @param ncores Number of threads used for parallel calculation. Default: 1.
#' @param verbose Gives information about each calculation step. Default: TRUE.
#' @param ... Other parameters passed to \code{\link{runViper}}, \code{\link{runAUCell}} and
#' \code{\link{runGSVA}}.
#'
#' @return
#'
#'
# SingleCellExperiment
#' @rdname regulonActivity
#' @export
setMethod("regulonActivity", "SingleCellExperiment", function(x,
                                                              gene_list=NULL,
                                                              method = c("viper","aucell","ssgsea"),
                                                              normalization = TRUE,
                                                              min_target_num = 5L,
                                                              max_target_num = 500L,
                                                              use_assay = 1L,
                                                              use_regulon = 1L,
                                                              norm_act = FALSE,
                                                              norm_coe = 1,
                                                              ncores = 1,
                                                              verbose = FALSE,
                                                              ...){
  method <- match.arg(method)
  if(is.null(gene_list)){
    gene_list <- as.list(Regulon(x,use_regulon))
  }
    exp_mat <- as.matrix(assay(x,i=use_assay))
  if(method == "aucell"){
    act_mat <- runAUCell(exp_mat = exp_mat,
                      gene_list = gene_list,
                      normalization = normalization,
                      ncores = ncores,
                      verbose = verbose, ...)
  }else if(method == "viper"){
    act_mat <- runViper(exp_mat = exp_mat,
                     gene_list = gene_list,
                     ncores = ncores,
                     normalization = normalization, ...)

  }else if(method %in% c("gsva","ssgsea")){
    act_mat <- runGSVA(exp_mat = exp_mat,
                    gene_list = gene_list,
                    method = method,
                    min_target_num = min_target_num,
                    max_target_num = max_target_num,
                    ...)
  }else{
    stop("Current version of 'regulonActivity()' only support 6 methods!")
  }
  # include TF activities into SingleCellExperiment object
  if(norm_act){
    act_mat <- t(apply(act_mat,1,function(x){(x-min(x))/{max(x)-min(x)}}))*norm_coe
  }
  atfr_sce <- SingleCellExperiment::SingleCellExperiment(as.matrix(act_mat))
  SummarizedExperiment::assayNames(atfr_sce) <- "atfr"
  SummarizedExperiment::colData(atfr_sce) <- SummarizedExperiment::colData(x)
  SingleCellExperiment::altExp(x, e = method) <- atfr_sce
  #---add regulons data into object x
  x
})
# matrix input
#'
#' @rdname regulonActivity
#' @export
setMethod("regulonActivity", "matrix", function(x,
                                                gene_list,
                                                method = c("viper","ssgsea","aucell"),
                                                normalization = TRUE,
                                                min_target_num = 5,
                                                max_target_num = 500,
                                                scale_act = FALSE,
                                                scale_coe = 10,
                                                ncores = 1,
                                                verbose = FALSE,
                                                ...){
  method <- match.arg(method)
  exp_mat <- as.matrix(x)
  gene_list <- as.list(gene_list)
  if(method == "aucell"){
    act_mat <- runAUCell(exp_mat = exp_mat,
                         gene_list = gene_list,
                         normalization = normalization,
                         ncores = ncores,
                         verbose = verbose, ...)
  }else if(method == "viper"){
    act_mat <- runViper(exp_mat = exp_mat,
                        gene_list = gene_list,
                        ncores = ncores,
                        normalization = normalization, ...)

  }else if(method %in% c("gsva","ssgsea","zscore", "plage")){
    act_mat <- runGSVA(exp_mat = exp_mat,
                       gene_list = gene_list,
                       method = method,
                       min_target_num = min_target_num,
                       max_target_num = max_target_num,
                       ...)
  }else{
    stop("Current version of 'regulonActivity()' only support 6 methods!")
  }
  if(norm_act){
    act_mat <- t(apply(act_mat,1,function(x){(x-min(x))/{max(x)-min(x)}}))*norm_coe
  }
  act_mat
})

######
#Viper

#' Running Viper with Transcription Factor Regulons
#'
#' @param exp_mat A \pkg{ExpressionSet} object or a gene by sample matrix with gene expression
#' values, such as TPMs, counts or UMIs.
#' @param gene_list List, gene regulons used for activate analysis.
#' @param min_target_num Minimum number of genes in your regulons.
#' @param normalization Logical, whether to normalize the activate score?
#' @param ncores Integer, number of cores to use (only 1 in Windows-based systems)
#' @param ... Others parameters passed to \code{viper}.

#'
#' @return A regulon by sample matrix with activate score.
#' @export
#'
#' @examples
runViper <- function(exp_mat,
                     gene_list,
                     min_target_num = 5,
                     normalization = TRUE,
                     ncores = 1,
                     ... ){
  #---perpare viper regulons
  viper_regulons <- lapply(gene_list, function(x) {
    if(is.numeric(x)){
      if(is.null(names(x))) stop("The weighted regulons must have names!")
      tfmode <- stats::setNames(rep(1, length(x)), names(x))
      out_list <- list(tfmode = tfmode, likelihood = unname(x))
    }else if(is.character(x)){
      tfmode <- stats::setNames(rep(1, length(x)), x)
      out_list <- list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    }else{
      stop("The values in 'gene_list' must be gene symbols or weighted gene symbols!")
    }
    out_list
  })
  if(Sys.info()["sysname"] == "Windows") ncores = 1
  viper_res <- viper::viper(eset = exp_mat,
                            regulon = viper_regulons,
                            nes = normalization,
                            minsize = min_target_num,
                            eset.filter = TRUE,
                            cores = ncores, ...)
  viper_res
}


######
#AUCell

#' Running AUCell with Transcription Factor Regulons
#'
#' @param exp_mat An expression matrix
#' @param gene_list A list of regulons.
#' @param ncores Number of cores used.
#' @param ... Other parameters passed to \code{\link{AUCell_buildRankings}} and
#' \code{\link{AUCell_calcAUC}}
#' @param normalization
#'
#' @return A matrix with regulons by cell activity matrix
#' @export
#'
runAUCell <- function(exp_mat,
                      gene_list,
                      ncores = 1,
                      normalization = TRUE,
                      verbose = interactive(),
                      ... ){
  gene_list <- lapply(gene_list,function(x){
    if(is.numeric(x)){
      x <- names(x)
    }else if(is.character(x)){
      x <- x
    }else{
      stop("The values in 'gene_list' must be gene symbols or weighted gene symbols!")
    }
    x
  })
  cells_rankings <- AUCell::AUCell_buildRankings(exprMat = exp_mat,
                                                 nCores=ncores,
                                                 plotStats=FALSE,
                                                 verbose = verbose,
                                                 ...)
  cells_AUC <- AUCell::AUCell_calcAUC(geneSets = gene_list,
                                      rankings = cells_rankings,
                                      nCores = 1,
                                      normAUC = normalization,
                                      verbose = verbose, ...)
  assay(cells_AUC)
}

######
#GSVA ssGSEA

#' Running GSVA with Transcription Factor Regulons
#'
#' @param exp_mat An expression matrix
#' @param gene_list A list of regulons.
#' @param method Method to employ in the estimation of gene-set enrichment scores per sample.
#' @param kcdf Character string denoting the kernel to use during the non-parametric estimation
#' of the cumulative distribution function of expression levels across samples when method="gsva".
#' @param min_target_num Minimum number of genes in your regulons.
#' @param max_target_num Maximum number of genes in your regulons.
#' @param ... Other parameters passed to \code{\link{gsva()}}.
#'
#' @return A regulon by sample matrix with activate score.
#' @export
#'
#' @examples
runGSVA <- function(exp_mat,
                    gene_list,
                    method = c("gsva", "ssgsea", "zscore", "plage"),
                    kcdf = c("Poisson", "Gaussian", "none"),
                    min_target_num = 5,
                    max_target_num = 500,
                    ...){
  method <- match.arg(method)
  kcdf <- match.arg(kcdf)
  gene_list <- lapply(gene_list,function(x){
    if(is.numeric(x)){
      x <- names(x)
    }else if(is.character(x)){
      x <- x
    }else{
      stop("The values in 'gene_list' must be gene symbols or weighted gene symbols!")
    }
    x
  })
  gsva_res <- GSVA::gsva(expr = exp_mat,
                         gset.idx.list = gene_list,
                         method = method, kcdf = kcdf,
                         min.sz=min_target_num,
                         max.sz=max_target_num, ...)
  gsva_res
}


######
#ssGSEA2.0

runSsgsea2 <- function(exp_mat,gene_list
                       ){

}

#------ssgsea2.0


