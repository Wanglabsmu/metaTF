% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimReduction.R
\name{regulonUMAP,matrix-method}
\alias{regulonUMAP,matrix-method}
\alias{regulonUMAP,SingleCellExperiment-method}
\title{Perform UMAP on regulon activation data}
\usage{
\S4method{regulonUMAP}{matrix}(
  x,
  use_regulon = NULL,
  ncomponents = 2,
  ntop = 500,
  scale = FALSE,
  transposed = FALSE,
  pca = if (transposed) NULL else 50,
  n_neighbors = 15,
  external_neighbors = FALSE,
  ...
)

\S4method{regulonUMAP}{SingleCellExperiment}(
  x,
  alt_assay = 1L,
  use_regulon = NULL,
  ncomponents = 2,
  ntop = 500,
  scale = FALSE,
  pca = if (!is.null(dimred)) NULL else 50,
  n_neighbors = 15,
  external_neighbors = FALSE,
  dimred = NULL,
  n_dimred = NULL,
  name = "UMAP",
  ...
)
}
\arguments{
\item{x}{A numeric \code{matrix} where rows are TF regulons and columns were
cells, Alternatively, a \code{\link{SingleCellExperiment}} object containing
such matrix in \code{altExp(x)} is also supported.}

\item{use_regulon}{Character vectors indicating the regulons used for UMAP
analysis. Same as \code{subset_row} in \code{\link{calculateUMAP()}}.}

\item{ncomponents}{Numeric scalar indicating the number of principal
components to obtain.}

\item{ntop}{Numeric scalar specifying the number of features with the highest
variances to use for dimensional reduction.}

\item{scale}{Logical scalar, should the expression values be standardized?}

\item{transposed}{Logical scalar, is x transposed with cells in rows?}

\item{pca}{Integer scalar specifying how many PCs should be used as input
into the UMAP algorithm. By default, no PCA is performed if the input is a
dimensional reduction result.}

\item{n_neighbors}{Integer scalar, number of nearest neighbors to identify
when constructing the initial graph.}

\item{external_neighbors}{Logical scalar indicating whether a nearest
neighbors search should be computed externally with \code{\link{findKNN}}.}

\item{...}{Other parameters passed to \code{\link{calculateUMAP}} and
\code{\link{runUMAP}}.}

\item{alt_assay}{String or integer scalar specifying which assay in
\code{altExp(x)} containing the input data.}

\item{dimred}{String or integer scalar specifying the existing dimensional
reduction results to use.}

\item{n_dimred}{Integer scalar or vector specifying the dimensions to use if
dimred is specified.}

\item{name}{String specifying the name to be used to store the result in
the \code{\link{reducedDims()}} of the output.}
}
\description{
Perform uniform manifold approximation and projection (UMAP) on regulons by cells
activation matrix
}
