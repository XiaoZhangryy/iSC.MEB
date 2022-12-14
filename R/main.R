#' @importFrom mclust Mclust mclustBIC 
mycluster <- function(Z, G, int.model='EEE', verbose=FALSE){
    # require(mclust)
    mclus2 <- mclust::Mclust(Z, G=G,modelNames =int.model ,verbose=verbose)
    return(mclus2)
}

parafun_int <- function(k, Z, Sigma_equal, Sigma_diag, seed=1, init.start=5, int.model='EEE', verbose=FALSE){
    loglik0 = -1e20
    for (i in 1:init.start) {
        set.seed(seed+(i-1)*10)
        mclus0 <- mycluster(Z, G=k, int.model =int.model ,verbose=verbose)
        if (mclus0$loglik > loglik0) {
            loglik0 <- mclus0$loglik
            mclus2 <- mclus0
        }
    }
    
    Mu0k <- t(mclus2$parameters$mean)
    Sigmak <- mclus2$parameters$variance$sigma
    if(Sigma_diag){
        Sigma0k <- array(0, dim=dim(Sigmak))
        for(kk in 1:k){
            diag(Sigma0k[,,kk]) <- diag(Sigmak[,,kk])
        }
    } else {
        Sigma0k <- Sigmak
    }
    if(Sigma_equal){
        for(kk in 1:k){
            Sigma0k[,,kk] <- apply(Sigma0k, c(1,2), mean)
        }
    }
    Pi0k <- mclus2$parameters$pro
    return(list(y0k = mclus2$classification, Mu0k=Mu0k, Sigma0k=Sigma0k, Pi0k=Pi0k))
}

#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster
get_initial_value4seqK <- function(
    Kseq, hZ, Sigma_equal, Sigma_diag, seed=1, init.start=5, int.model='EEE', verbose=FALSE, coreNum){
    nK = length(Kseq)
    if (nK>1 & coreNum>1) {
        ## at most use 80% CPU
        cores <- min(c(nK, 10, parallel::detectCores()*0.8, coreNum))
        if (Sys.info()[1]=="Windows") {
            cl <- parallel::makeCluster(cores) # memory can not be shared type in Windows.
        } else {
            cl <- parallel::makeCluster(cores, type='FORK') # memory-shared type in linux or Mac.
        }
    
        message("Starting parallel computing initial values...")
        parallel::clusterExport(cl, varlist = c("Mclust", "mclustBIC"), envir=environment())
        # Run
        intList <- parallel::parLapply(
            cl, X=Kseq, parafun_int, Z=hZ, Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag, seed=seed,
            init.start=init.start, int.model=int.model, verbose=verbose)
        parallel::stopCluster(cl)
    } else {
        intList <- list(
            parafun_int(Kseq, Z=hZ, Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag, seed=seed,
            init.start=init.start, int.model=int.model, verbose=verbose))
    }

    Mu0k     = lapply(intList, function(intListK) intListK$Mu0k)
    Sigma0k  = lapply(intList, function(intListK) intListK$Sigma0k)
    y0k      = lapply(intList, function(intListK) intListK$y0k)
    Pi0k     = lapply(intList, function(intListK) intListK$Pi0k)

    return(list(y0k = y0k, Mu0k=Mu0k, Sigma0k=Sigma0k, Pi0k=Pi0k))
}

#' Each iSCMEBResObj object has a number of slots which store information.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description Each iSCMEBResObj object has a number of slots which store information. Key slots to access are listed below.
#' \itemize{
#'   \item \code{posList} - The position matrix list for a iSCMEBResObj object.
#'   \item \code{paramList} - The model parameter settings for a iSCMEBResObj object.
#'   \item \code{fitList} - The details of iSC.MEB models.
#'   \item \code{project} - Name of the project.
#'   \item \code{reduction} - The dimension reduction result of iSC.MEB model.
#'   \item \code{idents} - The clustering result of iSC.MEB model.
#' }
setClass("iSCMEBResObj", slots=list(
  posList = "ANY", 
  paramList= "list", 
  fitList = "ANY",
  project = "character",
  reduction = "list",
  idents = "ANY"
) )

#' Fit a iSC.MEB model.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{fit.iscmeb} is used to fit a iSC.MEB model.
#' @export
#' @param VList A list of PCs of log-normalized gene expression matrix. The i-th element is a ni * npcs matrtix, where ni is the number of spots of sample i, and npcs is the number of PC. We provide this interface for those users who would like to define the PCs by their own. 
#' @param AdjList A list of adjacency matrix with class \code{dgCMatrix}. We provide this interface for those users who would like to define the adjacency matrix by their own.
#' @param K An integer or integer vector, specify the candidates of number of clusters.
#' @param beta_grid An optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach, defualt as a sequence starts from 0, ends with 5, increase by 0.2. 
#' @param maxIter_ICM An optional positive value, represents the maximum iterations of ICM (6 by default).
#' @param maxIter An optional positive value, represents the maximum iterations of EM (25 by default).
#' @param epsLogLik An optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.
#' @param verbose An optional logical value, whether output the information of the ICM-EM algorithm.
#' @param init.start An optional number of times to calculate the initial value (1 by default). When init.start is larger than 1, initial value will be determined by log likelihood of mclust results. 
#' @param int.model An optional string, specify which Gaussian mixture model is used in evaluting the initial values for iSC.MEB, default as "EEE"; and see \code{Mclust} for more models' names.
#' @param Sigma_equal An optional logical value, specify whether Sigmaks are equal, default as \code{FALSE}.
#' @param Sigma_diag An optional logical value, specify whether Sigmaks are diagonal matrices, default as \code{TRUE}.
#' @param seed An optional integer, the random seed in fitting iSC.MEB model.
#' @param coreNum An optional positive integer, means the number of thread used in parallel computating (1 by default).
#' @param criteria A string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "MAIC", "BIC" and "AIC" ("MBIC" by default).
#' @param c_penalty An optional positive value, the adjusted constant used in the MBIC criteria (1 by default).
#' @param pca.method A string specify the dimension reduction method used to generate \code{VList}. "PCA" by default. 
#'
#' @return Returns a iSCMEBResObj object which contains all model results.
#' @details iSCMEBResObj is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the body of our algorithm. 
#' @seealso \code{\link{runPCA}}, \code{\link{CreateNeighbors}}, \code{\link{iSCMEBResObj-class}}
#'
#' @importFrom stats cov
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' XList <- lapply(iSCMEBObj_toy@seulist, function(x) Matrix::t(x@assays$RNA@data))
#' VList <- runPCA(XList)
#' posList <- lapply(iSCMEBObj_toy@seulist, function(x) x@meta.data[,c("row", "col")])
#' AdjList <- CreateNeighbors(posList, platform = "Visium")
#' resList <- fit.iscmeb(VList, AdjList, K=7, maxIter=10)
fit.iscmeb <- function(
    VList, AdjList, K, beta_grid=seq(0, 5, by=0.2), maxIter_ICM=6, maxIter=25, 
    epsLogLik=1e-5, verbose=TRUE, int.model="EEE", init.start=1, 
    Sigma_equal = FALSE, Sigma_diag = TRUE, seed=1, coreNum=1, 
    criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty=1, pca.method="PCA") 
{
    error_heter <- TRUE
    M <- length(VList)
    q <- ncol(VList[[1]])
    K.order = order(K, decreasing = T)
    K = sort(K, decreasing = T)
    
    Psi_int <- array(dim=c(q,q, M))
    for( j in 1: M) Psi_int[,,j] <- cov(VList[[j]]/4)
    
    message("Evaluate initial values...")
    intlist <- get_initial_value4seqK(
        K, matlist2mat(VList), Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag, seed=seed, 
        init.start=init.start, int.model=int.model, verbose=verbose, coreNum=coreNum) 
    
    Mu_int    <- intlist$Mu0k 
    Sigma_int <- intlist$Sigma0k
    y_int     <- intlist$y0k
    
    message("Fit SC-MEB2...")
    Sp_embed <- TRUE
    result <- SepSpatialClusterCpp(
        VList, AdjList, y_int, Mu_int, Sigma_int, Psi_int, 1.5, 
        beta_grid, maxIter_ICM, maxIter, epsLogLik, verbose, 
        !error_heter, Sigma_equal, Sigma_diag, Sp_embed, max(K), min(K), coreNum) 
    
    result <- result[order(K.order)]
    K      <- K[order(K.order)]

    output <- new(
        Class = "iSCMEBResObj",
        posList = NULL,
        paramList= list(),
        fitList = result,
        project = "iSC.MEB",
        reduction = list(),
        idents = NULL
    )

    n  <- sum(sapply(VList, nrow))
    nT <- length(VList)

    output@paramList$dfs = sapply(K, function(K0) degree.freedom(K0, q, nT, Sigma_equal, Sigma_diag, Sp_embed))
    output@paramList$log_likes = unlist(lapply(result, function(fit) fit$loglik))
    output@paramList$MAIC <- -2.0*output@paramList$log_likes + output@paramList$dfs*2*log(log(q+n))*c_penalty
    output@paramList$AIC  <- -2.0*output@paramList$log_likes + output@paramList$dfs*2
    output@paramList$MBIC <- -2.0*output@paramList$log_likes + output@paramList$dfs*log(n)*log(log(q+n))*c_penalty
    output@paramList$BIC  <- -2.0*output@paramList$log_likes + output@paramList$dfs*log(n)

    criteria = match.arg(criteria)
    Mycriteria <- switch(
        criteria, 
        MAIC = output@paramList$MBIC,
        AIC  = output@paramList$AIC,
        MBIC = output@paramList$MBIC, 
        BIC  = output@paramList$BIC 
    )

    output@paramList$opt  = which.min(Mycriteria)
    output@paramList$optK = K[output@paramList$opt]
    if (is.null(names(VList))) {
        output@paramList$sample_name = paste0("Sample", c(1:nT))
    } else {
        output@paramList$sample_name = names(VList)
    }
    output@paramList$K  = K
    output@paramList$n  = sapply(VList, nrow)
    output@paramList$q  = q
    output@paramList$Sigma_diag  = Sigma_diag
    output@paramList$Sigma_equal = Sigma_equal
    output@paramList$Sp_embed    = Sp_embed
    output@paramList$nT = nT
    output@paramList$modelselect = switch(
        criteria, 
        MAIC = paste0("MAIC_", c_penalty),
        AIC  = "AIC",
        MBIC = paste0("MBIC_", c_penalty), 
        BIC  = "BIC" 
    )
    output@paramList$pca.method <- pca.method
    
    output@reduction[[pca.method]] <- VList
    output@reduction$iSCMEB = output@fitList[[output@paramList$opt]]$hZ
    output@idents = lapply(output@fitList[[output@paramList$opt]]$cluster, as.vector)
    
    return(output)
}

degree.freedom <- function(K, q, nT, Sigma_equal, Sigma_diag, Sp_embed) {
    if (Sigma_diag) {
        # message("Sigma is set to diagonal matrices.\n")
        df_Sigma = q
    } else {
        # message("Sigma is set to dense matrices.\n")
        df_Sigma = q*(q+1)/2.0
    }
    if (Sigma_equal) {
        df_Sigma = df_Sigma
    } else {
        df_Sigma = df_Sigma*K
    }

    if (Sp_embed) {
        df_psi = q*(q+1)/2.0
    } else {
        df_psi = 0
    }
    df_psi = df_psi*nT

    # Mu + Sigma + Psi + beta
    dfree <- K*q + df_Sigma + df_psi + nT
    
    return(dfree)
}

#' Fit a iSC.MEB model.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{iSCMEB} is used to fit a iSC.MEB model for a iSCMEBObj object.
#' @export
#' @param iSCMEBObj A iSCMEBObj object. 
#' @param K An optional integer or integer vector, specify the candidates of number of clusters. if \code{K=NULL}, it will be set to 5~12.
#'
#' @details The model results are saved in the slot of resList.
#'
#' @return Returns a revised iSCMEBObj object.
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{runPCA}}, \code{\link{CreateNeighbors}}, \code{\link{SetModelParameters}}, \code{\link{iSCMEBResObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' iSCMEBObj_toy <- runPCA(iSCMEBObj_toy)
#' iSCMEBObj_toy <- CreateNeighbors(iSCMEBObj_toy, platform = "Visium")
#' iSCMEBObj_toy <- SetModelParameters(iSCMEBObj_toy, maxIter=10)
#' iSCMEBObj_toy <- iSCMEB(iSCMEBObj_toy, K = 7)
iSCMEB <- function(iSCMEBObj, K=NULL){
    if(!inherits(iSCMEBObj, "iSCMEBObj")) 
        stop("iSCMEB: Check the argument: iSCMEBObj! iSCMEBObj must be a iSCMEBObj object.")
    
    if(is.null(K)) K <- 5:12
    if(is.null(iSCMEBObj@seulist)) stop("The slot seulist in iSCMEBObj is NULL!")
    
    ## Get pcs
    reduction.name <- iSCMEBObj@parameterList$reduction.name
    VList <- lapply(iSCMEBObj@seulist, function(seu) seu[[reduction.name]]@cell.embeddings)
    names(VList) = names(iSCMEBObj@seulist)
  
    # get parameters
    beta_grid <- iSCMEBObj@parameterList$beta_grid
    maxIter_ICM <- iSCMEBObj@parameterList$maxIter_ICM
    maxIter <- iSCMEBObj@parameterList$maxIter
    epsLogLik <- iSCMEBObj@parameterList$epsLogLik
    verbose <- iSCMEBObj@parameterList$verbose
    int.model <- iSCMEBObj@parameterList$int.model
    init.start <- iSCMEBObj@parameterList$init.start
    Sigma_equal <- iSCMEBObj@parameterList$Sigma_equal 
    Sigma_diag <- iSCMEBObj@parameterList$Sigma_diag
    seed <- iSCMEBObj@parameterList$seed
    coreNum <- iSCMEBObj@parameterList$coreNum
    criteria <- iSCMEBObj@parameterList$criteria
    c_penalty <- iSCMEBObj@parameterList$c_penalty
    pca.method <- iSCMEBObj@parameterList$pca.method
  
    ## Centering
    iSCMEBObj@resList <- fit.iscmeb(
        VList, 
        AdjList = iSCMEBObj@AdjList, 
        K=K, 
        beta_grid= beta_grid,
        maxIter_ICM=maxIter_ICM,
        maxIter=maxIter,
        epsLogLik=epsLogLik,
        verbose=verbose,
        int.model=int.model,
        init.start=init.start,
        Sigma_equal=Sigma_equal, 
        Sigma_diag=Sigma_diag,
        seed=seed,
        coreNum=coreNum,
        criteria=criteria,
        c_penalty=c_penalty,
        pca.method=pca.method)
    if (!is.null(iSCMEBObj@posList)) iSCMEBObj@resList@posList = iSCMEBObj@posList
    
    return(iSCMEBObj)
}

