#'  Create an object for transcription factor regulon analysis
#'
#' Create an \code{\link{SingleCellExperiment}} object for transcription factor regulon (TFR) analysis
#'
#' @param exp_data Expression matrix contains expression data for individual samples.
#' Could be non-normalized data such as counts and UMIs or normalized data such as TPM or FPKM.
#' @param col_data A \code{data.frame} contains sample annotation information, such as cell types,
#' tissues and conditions. The row names of `col_data` must be the same with the column names of
#' `exp_data`. Default: null.
#' @param row_data A \code{data.frame} contains feature annotation information. Could be gene names
#' or gene types. The row names of `row_data` must be the same with the row names of `exp_data`.
#' Default: null.
#' @param regulons A list of transcription factor regulons (TFRs).
#'
#' @return A \code{\link{SingleCellExperiment}} object.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
#' @export
#'
#' @examples
#' library(metaTF)
#' #load example data
#' expression_data <- readRDS(system.file("extdata", "mouse_HSC_formation_expression.rds", package = "metaTF"))
#' colData <- readRDS(system.file("extdata", "mouse_HSC_formation_colData.rds", package = "metaTF"))
#' #create atfr object.
#' atfr <- metaTF(exp_data = expression_data, col_data=colData)
#'
#'
metaTF <- function(exp_data, col_data=NULL, row_data=NULL, regulons=NULL){
  if(!(is.matrix(exp_data) || is(exp_data, "dgCMatrix"))){
    warnings("Your input exp_data is not matrix!")
  }
  if(is.null(col_data)){
    col_data=data.frame(sample=colnames(exp_data),
                        row.names = colnames(exp_data),
                        check.names = FALSE)
  }

  exp_data <- as.matrix(exp_data)
  if(!is.numeric(exp_data)) stop("There is non-numeric values in your input 'exp_data'!")
  if(any(is.na(exp_data))) stop("NA exists in your input 'exp_data'!")
  exp_data <- as(exp_data,"dgCMatrix")
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=exp_data),
                                                    rowData = row_data,
                                                    colData = col_data)
  metadata(sce)$metaTF <- S4Vectors::SimpleList(regulons=S4Vectors::SimpleList(),
                                                GRNs=S4Vectors::SimpleList(),
                                                atfrMetadata=S4Vectors::SimpleList())
  if(!is.null(regulons)){
    if(!is(regulons,"list")) stop("The regulons should be list!")
    Regulon(x = sce, i = "regulon") <- regulons
  }
  metadata(sce)$atfr_version <- Biobase::package.version("metaTF")
  sce
}

# import data from Seurat

#' Import single-cell data fro Seurat object
#'
#' @param object A \code{\link{Seurat}} object.
#' @param use_assay The assay used to extract the expression data. The raw counts
#' data usually stored in assay named "RNA".
#' @param add_logcounts Whether to add logcounts into object? Default: TRUE.
#' @param add_metadata Whether to add metadata into object? Default: TRUE.
#' @param add_reductions Whether to add reduced dimension into object? Such as:
#' PCA, tSNE and UMAP. Default: TRUE.
#'
#' @return An \code{\link{SingleCellExperiment}} object.
#' @importFrom Seurat Assays GetAssay Reductions
#' @importFrom SingleCellExperiment reducedDim logcounts
#' @export
#'
#' @examples
importFromSeurat <- function(object,
                             use_assay = "RNA",
                             add_logcounts = TRUE,
                             add_metadata = TRUE,
                             add_reductions = TRUE){
  if(!is(object,"Seurat")) stop("The input object is not a seurat object!")
  assay_names <- Seurat::Assays(object)
  if(!all(use_assay %in% assay_names)){
    stop("The assay: ",setdiff(use_assay, assay_names)," is not in your object!")
  }
  seu_obj <- Seurat::GetAssay(object = object, assay=use_assay)
  counts <- seu_obj@counts
  if(add_metadata){
    col_metadata <- object@meta.data
  }else{
    col_metadata <- NULL
  }
  atfr <- metaTF::metaTF(exp_data = as.matrix(counts), col_data = col_metadata)
  if(add_logcounts){
    logcounts <- seu_obj@data
    SingleCellExperiment::logcounts(atfr) <- as.matrix(logcounts)
  }
  #---add reductions
  if(add_reductions){
    reduction_names <- Seurat::Reductions(object = object)
    all_reduction <- c('pca', 'tsne', 'umap')
    if(!any(reduction_names %in% all_reduction)){
      stop("No reductions detected in your object! Please set ",
           "'add_reductions' to FALSE!")
    }
    if("pca" %in% reduction_names){
      PCA_res <- Reductions(object = object, slot = "pca")
      SingleCellExperiment::reducedDim(x = atfr, type = "PCA") <- PCA_res@cell.embeddings
      attr(SingleCellExperiment::reducedDim(atfr), "rotation") <- PCA_res@feature.loadings
    }
    if("tsne" %in% reduction_names){
      TSNE_res <- Reductions(object = object, slot = "tsne")
      SingleCellExperiment::reducedDim(x = atfr, type = "TSNE") <- TSNE_res@cell.embeddings
    }
    if("umap" %in% reduction_names){
      UMAP_res <- Reductions(object = object, slot = "umap")
      SingleCellExperiment::reducedDim(x = atfr, type = "UMAP") <- UMAP_res@cell.embeddings
    }
  }
  atfr
}








