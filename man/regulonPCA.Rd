% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimReduction.R
\name{regulonPCA,matrix-method}
\alias{regulonPCA,matrix-method}
\alias{regulonPCA,SingleCellExperiment-method}
\title{Perform PCA on regulon activation data}
\usage{
\S4method{regulonPCA}{matrix}(
  x,
  use_regulon = NULL,
  ncomponents = 50,
  ntop = 500,
  scale = FALSE,
  transposed = FALSE,
  ...
)

\S4method{regulonPCA}{SingleCellExperiment}(
  x,
  alt_assay = 1L,
  use_regulon = NULL,
  ncomponents = 50,
  ntop = 500,
  dimred = NULL,
  n_dimred = NULL,
  scale = FALSE,
  ...
)
}
\arguments{
\item{x}{A numeric \code{matrix} where rows are TF regulons and columns were
cells, Alternatively, a \code{\link{SingleCellExperiment}} object containing
such matrix in \code{altExp(x)} is also supported.}

\item{use_regulon}{Character vectors indicating the regulons used for PCA
analysis. Same as \code{subset_row} in \code{\link{calculatePCA()}}.}

\item{ncomponents}{Numeric scalar indicating the number of principal
components to obtain.}

\item{ntop}{Numeric scalar specifying the number of features with the highest
variances to use for dimensional reduction.}

\item{scale}{Logical scalar, should the expression values be standardized?}

\item{transposed}{Logical scalar, is x transposed with cells in rows?}

\item{...}{Other parameters passed to \code{\link{calculatePCA}} and
\code{\link{runPCA}}.}

\item{alt_assay}{String or integer scalar specifying which assay in
\code{altExp(x)} containing the input data.}

\item{dimred}{String or integer scalar specifying the existing dimensional
reduction results to use.}

\item{n_dimred}{Integer scalar or vector specifying the dimensions to use if
dimred is specified.}
}
\description{
Perform principal components analysis (PCA) on regulons by cells activation
matrix
}
