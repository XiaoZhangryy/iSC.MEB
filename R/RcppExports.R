# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

SepSpatialClusterCpp <- function(vList, Adjlist, yList_int, Mu_int, Sigma_int, Psi_int, beta_int, beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose, homo, Sigma_equal, Sigma_diag, Sp_embed, maxK, minK, coreNum) {
    .Call(`_iSC_MEB_SepSpatialClusterCpp`, vList, Adjlist, yList_int, Mu_int, Sigma_int, Psi_int, beta_int, beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose, homo, Sigma_equal, Sigma_diag, Sp_embed, maxK, minK, coreNum)
}

getneighborhood_fast <- function(x, radius) {
    .Call(`_iSC_MEB_getneighborhood_fast`, x, radius)
}

wpcaCpp <- function(X, nPCs, weighted = TRUE) {
    .Call(`_iSC_MEB_wpcaCpp`, X, nPCs, weighted)
}

