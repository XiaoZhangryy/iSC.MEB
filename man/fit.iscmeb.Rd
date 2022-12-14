% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{fit.iscmeb}
\alias{fit.iscmeb}
\title{Fit a iSC.MEB model.}
\usage{
fit.iscmeb(
  VList,
  AdjList,
  K,
  beta_grid = seq(0, 5, by = 0.2),
  maxIter_ICM = 6,
  maxIter = 25,
  epsLogLik = 1e-05,
  verbose = TRUE,
  int.model = "EEE",
  init.start = 1,
  Sigma_equal = FALSE,
  Sigma_diag = TRUE,
  seed = 1,
  coreNum = 1,
  criteria = c("MBIC", "MAIC", "BIC", "AIC"),
  c_penalty = 1,
  pca.method = "PCA"
)
}
\arguments{
\item{VList}{A list of PCs of log-normalized gene expression matrix. The i-th element is a ni * npcs matrtix, where ni is the number of spots of sample i, and npcs is the number of PC. We provide this interface for those users who would like to define the PCs by their own.}

\item{AdjList}{A list of adjacency matrix with class \code{dgCMatrix}. We provide this interface for those users who would like to define the adjacency matrix by their own.}

\item{K}{An integer or integer vector, specify the candidates of number of clusters.}

\item{beta_grid}{An optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach, defualt as a sequence starts from 0, ends with 5, increase by 0.2.}

\item{maxIter_ICM}{An optional positive value, represents the maximum iterations of ICM (6 by default).}

\item{maxIter}{An optional positive value, represents the maximum iterations of EM (25 by default).}

\item{epsLogLik}{An optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.}

\item{verbose}{An optional logical value, whether output the information of the ICM-EM algorithm.}

\item{int.model}{An optional string, specify which Gaussian mixture model is used in evaluting the initial values for iSC.MEB, default as "EEE"; and see \code{Mclust} for more models' names.}

\item{init.start}{An optional number of times to calculate the initial value (1 by default). When init.start is larger than 1, initial value will be determined by log likelihood of mclust results.}

\item{Sigma_equal}{An optional logical value, specify whether Sigmaks are equal, default as \code{FALSE}.}

\item{Sigma_diag}{An optional logical value, specify whether Sigmaks are diagonal matrices, default as \code{TRUE}.}

\item{seed}{An optional integer, the random seed in fitting iSC.MEB model.}

\item{coreNum}{An optional positive integer, means the number of thread used in parallel computating (1 by default).}

\item{criteria}{A string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "MAIC", "BIC" and "AIC" ("MBIC" by default).}

\item{c_penalty}{An optional positive value, the adjusted constant used in the MBIC criteria (1 by default).}

\item{pca.method}{A string specify the dimension reduction method used to generate \code{VList}. "PCA" by default.}
}
\value{
Returns a iSCMEBResObj object which contains all model results.
}
\description{
The function \code{fit.iscmeb} is used to fit a iSC.MEB model.
}
\details{
iSCMEBResObj is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the body of our algorithm.
}
\examples{
data(iSCMEBObj_toy)
library(Seurat)
XList <- lapply(iSCMEBObj_toy@seulist, function(x) Matrix::t(x@assays$RNA@data))
VList <- runPCA(XList)
posList <- lapply(iSCMEBObj_toy@seulist, function(x) x@meta.data[,c("row", "col")])
AdjList <- CreateNeighbors(posList, platform = "Visium")
resList <- fit.iscmeb(VList, AdjList, K=7, maxIter=10)
}
\seealso{
\code{\link{runPCA}}, \code{\link{CreateNeighbors}}, \code{\link{iSCMEBResObj-class}}
}
