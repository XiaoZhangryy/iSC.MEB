% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization.r
\name{LowEmbedPlot}
\alias{LowEmbedPlot}
\alias{LowEmbedPlot.iSCMEBObj}
\alias{LowEmbedPlot.iSCMEBResObj}
\title{Plot low-dimensional embeddings}
\usage{
LowEmbedPlot(
  obj,
  item = c("cluster", "batch"),
  reduction = NULL,
  combine = TRUE,
  cols = NULL,
  layout.dim = c(1, 2),
  common.legend = TRUE,
  ...
)

\method{LowEmbedPlot}{iSCMEBObj}(
  obj,
  item = c("cluster", "batch"),
  reduction = NULL,
  combine = TRUE,
  cols = NULL,
  layout.dim = c(1, 2),
  common.legend = TRUE,
  ...
)

\method{LowEmbedPlot}{iSCMEBResObj}(
  obj,
  item = c("cluster", "batch"),
  reduction = NULL,
  combine = TRUE,
  cols = NULL,
  layout.dim = c(1, 2),
  common.legend = TRUE,
  ...
)
}
\arguments{
\item{obj}{A iSCMEBObj object or a iSCMEBResObj object.}

\item{item}{Which feature to be used in the plot.}

\item{reduction}{The name of embeddings to be used in the plot. If reduction is null, the last added one is used for plotting.}

\item{combine}{An indicator of whether plot all on a figure. If TRUE, all figures are plotted; otherwise, return a list with each plot as component. TRUE by default.}

\item{cols}{A vector determine the colors used in the plot.}

\item{layout.dim}{The dimension in the layout of plots when \code{combine = TRUE}.}

\item{common.legend}{An indicator of whether combine the legend of multiple plots. TRUE by default.}

\item{...}{Arguments passed to other methods.}
}
\value{
Return a ggplot2 object.

Return a ggplot2 object.
}
\description{
The function \code{LowEmbedPlot} is used to visualize the low-dimensional embeddings vs cluster or batch to demonstrate the clustering performance and batch remove performance.
}
\examples{
data(iSCMEBObj_toy)
library(Seurat)
p1 <- LowEmbedPlot(iSCMEBObj_toy, item="cluster", reduction="TSNE2")
p2 <- LowEmbedPlot(iSCMEBObj_toy, item="batch", reduction="TSNE2")
}
\seealso{
\code{\link{iSCMEBObj-class}}, \code{\link{plot_scatter}}, \code{\link{plot_RGB}}, \code{\link{fit.iscmeb}}, \code{\link{iSCMEBResObj-class}}
}
