% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createAtfrObject.R
\name{importFromSeurat}
\alias{importFromSeurat}
\title{Import single-cell data fro Seurat object}
\usage{
importFromSeurat(
  object,
  use_assay = "RNA",
  add_logcounts = TRUE,
  add_metadata = TRUE,
  add_reductions = TRUE
)
}
\arguments{
\item{object}{A \code{\link{Seurat}} object.}

\item{use_assay}{The assay used to extract the expression data. The raw counts
data usually stored in assay named "RNA".}

\item{add_logcounts}{Whether to add logcounts into object? Default: TRUE.}

\item{add_metadata}{Whether to add metadata into object? Default: TRUE.}

\item{add_reductions}{Whether to add reduced dimension into object? Such as:
PCA, tSNE and UMAP. Default: TRUE.}
}
\value{
An \code{\link{SingleCellExperiment}} object.
}
\description{
Import single-cell data fro Seurat object
}
