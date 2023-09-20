#' @export
setGeneric("regulonPCA", function(x, ...) standardGeneric("regulonPCA"))

#' Perform PCA on regulon activation data
#'
#' Perform principal components analysis (PCA) on regulons by cells activation
#' matrix
#'
#' @param x A numeric \code{matrix} where rows are TF regulons and columns were
#' cells, Alternatively, a \code{\link{SingleCellExperiment}} object containing
#' such matrix in \code{altExp(x)} is also supported.
#' @param use_regulon Character vectors indicating the regulons used for PCA
#' analysis. Same as \code{subset_row} in \code{\link{calculatePCA()}}.
#' @param ncomponents Numeric scalar indicating the number of principal
#' components to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest
#' variances to use for dimensional reduction.
#' @param scale Logical scalar, should the expression values be standardized?
#' @param transposed Logical scalar, is x transposed with cells in rows?
#' @param alt_assay String or integer scalar specifying which assay in
#' \code{altExp(x)} containing the input data.
#' @param dimred String or integer scalar specifying the existing dimensional
#' reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified.
#' @param ... Other parameters passed to \code{\link{calculatePCA}} and
#' \code{\link{runPCA}}.
#'
#' @importFrom scater calculatePCA runPCA
#' @importFrom SingleCellExperiment altExp
#'
#'
#' @rdname regulonPCA
#' @export
setMethod("regulonPCA", "matrix",function(x,
                                          use_regulon = NULL,
                                          ncomponents = 50,
                                          ntop = 500,
                                          scale=FALSE,
                                          transposed = FALSE,
                                          ...){
  reduced_dims <- scater::calculatePCA(x = x,
                                       ncomponents = ncomponents,
                                       ntop = ntop,
                                       subset_row = use_regulon,
                                       scale = scale)
  reduced_dims
})

#' @rdname regulonPCA
#' @export
setMethod("regulonPCA", "SingleCellExperiment",function(x,
                                                        alt_assay = 1L,
                                                        use_regulon = NULL,
                                                        ncomponents = 50,
                                                        ntop = 500,
                                                        dimred = NULL,
                                                        n_dimred = NULL,
                                                        scale=FALSE,
                                                        ...){
  atfr <- SingleCellExperiment::altExp(x = x,e = alt_assay)
  atfr <- scater::runPCA(x = atfr,
                         ncomponents = ncomponents,
                         ntop = ntop,
                         subset_row = use_regulon,
                         scale = scale,
                         dimred = dimred,
                         n_dimred = n_dimred,
                         exprs_values = "atfr")
  SingleCellExperiment::altExp(x = x,e = alt_assay,withDimnames=TRUE) <- atfr
  x
})

#' @export
setGeneric("regulonTSNE", function(x, ...) standardGeneric("regulonTSNE"))

#' Perform t-SNE on regulon activation data
#'
#' Perform t-stochastic neighbour embedding (t-SNE) on regulons by cells
#' activation matrix
#'
#' @param x A numeric \code{matrix} where rows are TF regulons and columns were
#' cells, Alternatively, a \code{\link{SingleCellExperiment}} object containing
#' such matrix in \code{altExp(x)} is also supported.
#' @param use_regulon Character vectors indicating the regulons used for TSNE
#' analysis. Same as \code{subset_row} in \code{\link{calculateTSNE()}}.
#' @param ncomponents Numeric scalar indicating the number of principal
#' components to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest
#' variances to use for dimensional reduction.
#' @param transposed Logical scalar, is x transposed with cells in rows?
#' @param scale Logical scalar, should the expression values be standardized?
#' @param perplexity Numeric scalar defining the perplexity parameter,
#' see \code{\link{?Rtsne}} for more details.
#' @param normalize Logical scalar indicating if input values should be scaled
#' for numerical precision, see \code{\link{normalize_input}}.
#' @param theta Numeric scalar specifying the approximation accuracy of the
#' Barnes-Hut algorithm, see \code{\link{Rtsne}} for details.
#' @param external_neighbors Logical scalar indicating whether a nearest
#' neighbors search should be computed externally with \code{\link{findKNN}}.
#' @param use_fitsne Logical scalar indicating whether \code{\link{fitsne}}
#' should be used to perform t-SNE.
#' @param alt_assay String or integer scalar specifying which assay in
#' \code{altExp(x)} containing the input data.
#' @param dimred String or integer scalar specifying the existing dimensional
#' reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified.
#' @param name String specifying the name to be used to store the result in
#' the \code{\link{reducedDims()}} of the output.
#' @param ... Other parameters passed to \code{\link{calculateTSNE}} and
#' \code{\link{runTSNE}}.
#'
#' @importFrom scater calculateTSNE runTSNE
#' @importFrom SingleCellExperiment altExp
#'
#' @rdname regulonTSNE
#' @export
setMethod("regulonTSNE", "matrix",function(x,
                                           use_regulon = NULL,
                                           ncomponents = 2,
                                           ntop = 500,
                                           scale=FALSE,
                                           transposed = FALSE,
                                           perplexity = NULL,
                                           normalize = TRUE,
                                           theta = 0.5,
                                           external_neighbors = FALSE,
                                           use_fitsne = FALSE,
                                           ...){
  reduced_dims <- scater::calculateTSNE(x = x,
                                       ncomponents = ncomponents,
                                       ntop = ntop,
                                       subset_row = use_regulon,
                                       scale = scale,
                                       perplexity = perplexity,
                                       normalize = normalize,
                                       theta = theta,
                                       external_neighbors = external_neighbors,
                                       use_fitsne = use_fitsne, ...)
  reduced_dims
})

#' @rdname regulonTSNE
#' @export
setMethod("regulonTSNE", "SingleCellExperiment",function(x,
                                                         alt_assay = 1L,
                                                         use_regulon = NULL,
                                                         ncomponents = 2,
                                                         ntop = 500,
                                                         scale = FALSE,
                                                         perplexity = NULL,
                                                         normalize = TRUE,
                                                         theta = 0.5,
                                                         external_neighbors = FALSE,
                                                         use_fitsne = FALSE,
                                                         dimred = NULL,
                                                         n_dimred = NULL,
                                                         name = "TSNE", ...){
  atfr <- SingleCellExperiment::altExp(x = x,e = alt_assay)
  atfr <- scater::runTSNE(x = atfr,
                          ncomponents = ncomponents,
                          ntop = ntop,
                          subset_row = use_regulon,
                          scale = scale,
                          perplexity = perplexity,
                          normalize=normalize,
                          theta=theta,
                          external_neighbors=external_neighbors,
                          use_fitsne=use_fitsne,
                          dimred = dimred,
                          n_dimred = n_dimred,
                          exprs_values = "atfr",
                          name =name, ...)
  SingleCellExperiment::altExp(x = x,e = alt_assay,withDimnames=TRUE) <- atfr
  x
})

#' @export
setGeneric("regulonUMAP", function(x, ...) standardGeneric("regulonUMAP"))

#' Perform UMAP on regulon activation data
#'
#' Perform uniform manifold approximation and projection (UMAP) on regulons by cells
#' activation matrix
#'
#' @param x A numeric \code{matrix} where rows are TF regulons and columns were
#' cells, Alternatively, a \code{\link{SingleCellExperiment}} object containing
#' such matrix in \code{altExp(x)} is also supported.
#' @param use_regulon Character vectors indicating the regulons used for UMAP
#' analysis. Same as \code{subset_row} in \code{\link{calculateUMAP()}}.
#' @param ncomponents Numeric scalar indicating the number of principal
#' components to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest
#' variances to use for dimensional reduction.
#' @param scale Logical scalar, should the expression values be standardized?
#' @param transposed Logical scalar, is x transposed with cells in rows?
#' @param pca Integer scalar specifying how many PCs should be used as input
#' into the UMAP algorithm. By default, no PCA is performed if the input is a
#' dimensional reduction result.
#' @param n_neighbors Integer scalar, number of nearest neighbors to identify
#' when constructing the initial graph.
#' @param external_neighbors Logical scalar indicating whether a nearest
#' neighbors search should be computed externally with \code{\link{findKNN}}.
#' @param alt_assay String or integer scalar specifying which assay in
#' \code{altExp(x)} containing the input data.
#' @param dimred String or integer scalar specifying the existing dimensional
#' reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified.
#' @param name String specifying the name to be used to store the result in
#' the \code{\link{reducedDims()}} of the output.
#' @param ... Other parameters passed to \code{\link{calculateUMAP}} and
#' \code{\link{runUMAP}}.
#'
#' @importFrom scater calculateUMAP runUMAP
#' @importFrom SingleCellExperiment altExp
#'
#' @rdname regulonUMAP
#' @export
setMethod("regulonUMAP", "matrix",function(x,
                                           use_regulon = NULL,
                                           ncomponents = 2,
                                           ntop = 500,
                                           scale=FALSE,
                                           transposed = FALSE,
                                           pca = if (transposed) NULL else 50,
                                           n_neighbors=15,
                                           external_neighbors = FALSE,
                                           ...){
  reduced_dims <- scater::calculateUMAP(x = x,
                                        ncomponents = ncomponents,
                                        ntop = ntop,
                                        subset_row = use_regulon,
                                        scale = scale,
                                        pca=pca,
                                        n_neighbors=n_neighbors,
                                        external_neighbors = external_neighbors,
                                        ...)
  reduced_dims
})

#' @rdname regulonUMAP
#' @export
setMethod("regulonUMAP", "SingleCellExperiment",function(x,
                                                         alt_assay = 1L,
                                                         use_regulon = NULL,
                                                         ncomponents = 2,
                                                         ntop = 500,
                                                         scale = FALSE,
                                                         pca = if (!is.null(dimred)) NULL else 50,
                                                         n_neighbors=15,
                                                         external_neighbors = FALSE,
                                                         dimred = NULL,
                                                         n_dimred = NULL,
                                                         name = "UMAP",...){
  atfr <- SingleCellExperiment::altExp(x = x,e = alt_assay)
  atfr <- scater::runUMAP(x = atfr,
                          ncomponents = ncomponents,
                          ntop = ntop,
                          subset_row = use_regulon,
                          scale = scale,
                          external_neighbors=external_neighbors,
                          dimred = dimred,
                          n_dimred = n_dimred,
                          exprs_values = "atfr",
                          name =name, ...)
  SingleCellExperiment::altExp(x = x, e = alt_assay, withDimnames=TRUE) <- atfr
  x
})


