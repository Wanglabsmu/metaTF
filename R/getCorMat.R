setGeneric("getCorMat", function(x, ...) standardGeneric("getCorMat"))

#' Calculating the weighted TF-Target values
#'
#' Calculating the weighted TF-Target values for \code{\link{ATFinder}} object
#' @param x An scATGR object.
#' @param use_assay The i-th assays used for calculating the co-expression network.
#' It should be an gene expression matrix with raw counts or normalized values
#' such as TPM or FPKM. Using \code{\link{assays()}} to find assays in your
#' \code{\link{scATGR}} object. Default: 1.
#' @param use_regulator Whether to use regulators for enrichment analysis?
#' @param method Methods used for weighted matrix calculation. Current version only
#' support four types of methods. The "spearman", "pearson" and "kendall" were common
#' used method in \code{\link{cor()}} function. The "genie3" for \code{\link{GENIE3}}
#' method from package \code{\link{GENIE3}}. The "pcor-spearman" for \code{\link{ppcor}}.
#' Default: "spearman".
#' @param melt_cor_mat Whether to melt cor matrix? Default:TRUE
#' @param ncores Number of cores used for weight cor matrix calculation.
#' @param verbose Whether to show messages?
#' @param ... Others parameters passed to \code{\link{GENIE3}} and \code{\link{cor}}.
#'
#' @importFrom GENIE3 GENIE3
#' @importFrom SummarizedExperiment assays
#' @examples
#'
#'
#------for objects with 'ATFinder' type.
#' @rdname getCorMat
#' @export
setMethod("getCorMat", "matrix", function(x,
                                            use_assay = 1L,
                                            method = c("spearman","genie3", "pcor-spearman", "sclink", "pearson", "kendall"),
                                            use_regulator = TRUE,
                                            melt_cor_mat = FALSE,
                                            ncores = 1,
                                            verbose = interactive(),...) {
  method <- match.arg(method)
  if(!is.integer(use_assay)){
    stop("Parameter 'use_assay' must be an integer.")
  }else if(length(assays(x))<use_assay){
    stop("Only ", length(assays(x)), " assays exists in your object, but you want to use the ",
         use_assay," assay.")
  }
  exp_mat <- assays(x)[[as.integer(use_assay)]]
  if(use_regulator){
    regulator <- TFs(x)
    regulator <- intersect(row.names(exp_mat),regulator)
    if(length(regulator)==0) stop("No regulators found in your expression matrix!")
  }else{
    regulator <- NULL
  }
  cor_mat <- calcuWeightExpMat(exp_mat,
                               use_regulator = regulator,
                               method = method,
                               ncores = ncores,
                               verbose = verbose,
                               ...)
  #------melt cor matrix
  if(melt_cor_mat) {
    cor_mat_tab <- reshape2::melt(cor_mat)
    colnames(cor_mat_tab) <- c("regulator", "target", "weight")
    cor_mat_tab <- data.table::as.data.table(cor_mat_tab)
    data.table::setkey(cor_mat_tab,regulator, target)
    cor_mat <- cor_mat_tab[,.SD[order(weight,decreasing=TRUE),],by=regulator]
  }
  x@rankAssays[['weightCorMat']] <- cor_mat
  return(x)
})

