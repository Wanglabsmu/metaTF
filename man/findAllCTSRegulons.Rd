% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTypeSpecific.R
\name{findAllCTSRegulons,matrix-method}
\alias{findAllCTSRegulons,matrix-method}
\alias{findAllCTSRegulons,SingleCellExperiment-method}
\title{Calculating the significant enriched regulons for each cluster}
\usage{
\S4method{findAllCTSRegulons}{matrix}(
  x,
  col_data,
  clust_name,
  fdr = 0.05,
  css = 0.3,
  scale = FALSE,
  nperm = 1000,
  ncores = 1,
  verbose = interactive()
)

\S4method{findAllCTSRegulons}{SingleCellExperiment}(
  x,
  clust_name,
  use_altExp = TRUE,
  use_assay = 1L,
  fdr = 0.05,
  css = 0.3,
  expression_per = 0.1,
  scale = FALSE,
  nperm = 1000,
  ncores = 1,
  verbose = interactive()
)
}
\arguments{
\item{x}{An object with gene expression data used for regulon activity evaluation.
Could be \code{matrix}, \code{SingleCellExperiment}.}

\item{col_data}{A \code{data.frame} with cell cluster information. Each row
represent a cell and each column represent a cell metadata. The row names of
col_data must be identical with column names of x,}

\item{clust_name}{One column name in colData(x)(for \code{SingleCellExperiment})
or col_data (for \code{matrix}) specific the clusters used for cell-type
specific score calculation.}

\item{fdr}{The FDR cutoff used for filtering significant regulons. Default:0.05.}

\item{css}{The CSS cutoff used for filtering significant regulons. Default: 0.3.}

\item{scale}{Whether to scale values for calculation. If the activate scores
have negative values, it must be scaled. Default: FALSE.}

\item{nperm}{Number of permutations used for calculating p-values against CSS.
Default: 1000.}

\item{ncores}{Number of threads used for parallel calculation. Default: 1.}

\item{verbose}{Whether to output message? Default: TRUE.}

\item{use_altExp}{Whether to use altExps in \code{SingleCellExperiment}.
Default: TRUE. The activate matrix stored in altExps(x). So when perform
cell-type specific analysis on regulon activate matrix, this parameter must be
TRUE. When performing expression level analysis, set this parameter to FALSE.}

\item{use_assay}{Setting the assays used for cell-type specific analysis. For
regulon analysis, it could be 'viper', 'aucell' or 'ssgsea'. For gene expression
level, it could be 'counts' or 'logcounts'. Default: the first assay.}

\item{expression_per}{The expression percentage cutoff for filtering
significant regulons. Default:0.1.}
}
\value{
A \code{data.frame} with differential CSS for each TFRs.
}
\description{
Calculating the significant enriched regulons for each cluster
}
