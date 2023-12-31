% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doRadvizPlot.R
\name{RadvizFeaturePlot,data.frame-method}
\alias{RadvizFeaturePlot,data.frame-method}
\title{RadvizFeaturePlot}
\usage{
\S4method{RadvizFeaturePlot}{data.frame}(
  x,
  anchors = colnames(x),
  use_feature,
  return_rv = FALSE,
  trim_ratio = 0.005,
  norm_data = TRUE,
  optim_anchor = TRUE,
  recenter = NULL,
  text_color = "orangered4",
  text_size = 8,
  outline_circle = TRUE,
  point_size = 1,
  point_alpha = 0.8,
  point_shape = 16,
  plot_type = c("point", "hexagonal"),
  colors = scale_color_gradient(low = "grey80", high = "dodgerblue4"),
  ...
)
}
\arguments{
\item{x}{A \code{\link{data.frame}} or objects used for \code{\link{Radviz}} plot.}

\item{anchors}{Dimensions used to generate the map.}

\item{use_feature}{Features used to split the screens.}

\item{return_rv}{Only return \code{\link{getRadviz}} object, not generate figures.}

\item{optim_anchor}{Whether to optimize the anchor order?}

\item{recenter}{User defined start anchor.}

\item{text_color}{Text color, inherited from \code{\link{do.radviz}}.}

\item{text_size}{Text size.}

\item{outline_circle}{Whether to add circle in figure?}

\item{plot_type}{Figure types to plot.}

\item{...}{Other parameters passed to \code{\link{getRadviz}}.}

\item{color_by}{Groups used to label the color.}

\item{label_size}{Label size of point in \code{\link{geom_point()}}.}

\item{label_shape}{Label shape of point in \code{\link{geom_point()}}.}
}
\description{
RadvizFeaturePlot
}
