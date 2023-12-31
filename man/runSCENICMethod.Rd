% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runSCENIC.R
\name{runSCENICMethod}
\alias{runSCENICMethod}
\title{Running simplified version of SCENIC}
\usage{
runSCENICMethod(
  exprMatr,
  grn_tab,
  motif_ranking,
  motif_annotation,
  top_fraction = 0.05,
  nesThreshold = 0L,
  minSize = 5L,
  maxSize = 500L,
  verbose = interactive()
)
}
\arguments{
\item{exprMatr}{The expression matrix with normalized values where rows were
genes and columns were cells.}

\item{grn_tab}{The weighted TF regulatory matrix where rows were TF genes and 
columns were target genes.}

\item{motif_ranking}{The motif ranking object including '.feather' files 
import from \pkg{RcisTarget} database using \code{\link{importRankings()}} 
method.}

\item{motif_annotation}{The motif annotation database containing the 
annotations of the motif to transcription factors.}

\item{top_fraction}{The percentage of top ranked targets used for initial 
regulon analysis. Default: 0.05.}

\item{nesThreshold}{The normalized enrichment score for motif enrichment result. 
Default: 0, all enriched result.}

\item{minSize}{The minimum number of genes in one TF regulon.}

\item{maxSize}{The maximum number of genes in one TF regulon.}

\item{verbose}{Print message while runing?}
}
\value{
A matrix contains TF regulon activity in each single-cell.
}
\description{
This method aims to obtain the TF regulatory active matrix by running the 
RcisTarget and AUCell from SCENIC.
}
