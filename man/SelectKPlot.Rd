% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization.r
\name{SelectKPlot}
\alias{SelectKPlot}
\alias{SelectKPlot.iSCMEBObj}
\alias{SelectKPlot.iSCMEBResObj}
\title{SelectKPlot.}
\usage{
SelectKPlot(obj, criteria = c("InUse", "MBIC", "MAIC", "BIC", "AIC"))

\method{SelectKPlot}{iSCMEBObj}(obj, criteria = c("InUse", "MBIC", "MAIC", "BIC", "AIC"))

\method{SelectKPlot}{iSCMEBResObj}(obj, criteria = c("InUse", "MBIC", "MAIC", "BIC", "AIC"))
}
\arguments{
\item{obj}{An iSCMEBObj object or a iSCMEBResObj object.}

\item{criteria}{A string, specify the criteria used for selecting the number of clusters, supporting "InUse", "MBIC", "MAIC", "BIC" and "AIC" ("InUse" by default).}
}
\value{
Return a ggplot2 object.

Return a ggplot2 object.
}
\description{
The function \code{SelectKPlot} is used to demonstrate the scatter plot of the criteria used vs K for selecting the best K.
}
