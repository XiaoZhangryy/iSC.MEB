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

#' Fit a ISCMEB model.
#'
#' @useDynLib SC.MEB2, .registration = TRUE
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
#' @param int.model An optional string, specify which Gaussian mixture model is used in evaluting the initial values for ISCMEB, default as "EEE"; and see \code{Mclust} for more models' names.
#' @param Sigma_equal An optional logical value, specify whether Sigmaks are equal, default as \code{FALSE}.
#' @param Sigma_diag An optional logical value, specify whether Sigmaks are diagonal matrices, default as \code{TRUE}.
#' @param seed An optional integer, the random seed in fitting ISCMEB model.
#' @param coreNum An optional positive integer, means the number of thread used in parallel computating (1 by default).
#' @param criteria A string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "MAIC", "BIC" and "AIC" ("MBIC" by default).
#' @param c_penalty An optional positive value, the adjusted constant used in the MBIC criteria (1 by default).
#'
#' @return Returns a SCMEB_Result_Object object which contains all model results.
#' @seealso \code{\link{runPCA}}, \code{\link{CreateNeighbors}}
#'
#' @importFrom stats cov
#'
#' @examples
#' data(ISCMEBObj_toy)
#' library(Seurat)
#' XList <- lapply(ISCMEBObj_toy@seulist, function(x) Matrix::t(x@assays$RNA@data))
#' VList <- runPCA(XList)
#' posList <- lapply(ISCMEBObj_toy@seulist, function(x) x@meta.data[,c("row", "col")])
#' AdjList <- CreateNeighbors(posList, platform = "Visium")
#' resList <- SC.MEB2(VList, AdjList, K=7, maxIter=10)
SC.MEB2 <- function(
    VList, AdjList, K, beta_grid=seq(0, 5, by=0.2), maxIter_ICM=6, maxIter=25, 
    epsLogLik=1e-5, verbose=TRUE, int.model="EEE", init.start=1, 
    Sigma_equal = FALSE, Sigma_diag = TRUE, seed=1, coreNum=1, 
    criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty=1) 
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

    output = vector("list")
    output$fit = result

    n  <- sum(sapply(VList, nrow))
    nT <- length(VList)

    output$dfs = sapply(K, function(K0) degree.freedom(K0, q, nT, Sigma_equal, Sigma_diag, Sp_embed))
    output$log_likes = unlist(lapply(result, function(fit) fit$loglik))
    output$MAIC <- -2.0*output$log_likes + output$dfs*2*log(log(q+n))*c_penalty
    output$AIC  <- -2.0*output$log_likes + output$dfs*2
    output$MBIC <- -2.0*output$log_likes + output$dfs*log(n)*log(log(q+n))*c_penalty
    output$BIC  <- -2.0*output$log_likes + output$dfs*log(n)

    criteria = match.arg(criteria)
    Mycriteria <- switch(
        criteria, 
        MAIC = output$MBIC,
        AIC  = output$AIC,
        MBIC = output$MBIC, 
        BIC  = output$BIC 
    )

    output$opt  = which.min(Mycriteria)
    output$optK = K[output$opt]
    output$optSolution = output$fit[[output$opt]]

    param <- list(
        K  = K,
        n  = sapply(VList, nrow),
        q  = q,
        Sigma_diag  = Sigma_diag,
        Sigma_equal = Sigma_equal,
        Sp_embed    = Sp_embed,
        nT = nT
    )
    param$modelselect = switch(
        criteria, 
        MAIC = paste0("MAIC_", c_penalty),
        AIC  = "AIC",
        MBIC = paste0("MBIC_", c_penalty), 
        BIC  = "BIC" 
    )
    output$param <- param

    class(output) <- "SCMEB_Result_Object"
    
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

#' Fit a ISCMEB model.
#'
#' @useDynLib SC.MEB2, .registration = TRUE
#' @export
#' @param ISCMEBObj A ISCMEBObj object. 
#' @param K An optional integer or integer vector, specify the candidates of number of clusters. if \code{K=NULL}, it will be set to 4~12.
#'
#' @details The model results are saved in the slot of resList.
#'
#' @return Returns a revised ISCMEBObj object.
#' @seealso \code{\link{ISCMEBObj-class}}, \code{\link{runPCA}}, \code{\link{CreateNeighbors}}, \code{\link{SetModelParameters}}
#'
#' @examples
#' data(ISCMEBObj_toy)
#' library(Seurat)
#' ISCMEBObj_toy <- runPCA(ISCMEBObj_toy)
#' ISCMEBObj_toy <- CreateNeighbors(ISCMEBObj_toy, platform = "Visium")
#' ISCMEBObj_toy <- SetModelParameters(ISCMEBObj_toy, maxIter=10)
#' ISCMEBObj_toy <- ISCMEB(ISCMEBObj_toy, K = 7)
ISCMEB <- function(ISCMEBObj, K=NULL){
    if(!inherits(ISCMEBObj, "ISCMEBObj")) 
        stop("ISCMEB: Check the argument: ISCMEBObj! ISCMEBObj must be a ISCMEBObj object.")
    
    if(is.null(K)) K <- 4:12
    if(is.null(ISCMEBObj@seulist)) stop("The slot seulist in ISCMEBObj is NULL!")
    
    ## Get pcs
    reduction.name <- ISCMEBObj@parameterList$reduction.name
    VList <- lapply(ISCMEBObj@seulist, function(seu) seu[[reduction.name]]@cell.embeddings)
  
    # get parameters
    beta_grid <- ISCMEBObj@parameterList$beta_grid
    maxIter_ICM <- ISCMEBObj@parameterList$maxIter_ICM
    maxIter <- ISCMEBObj@parameterList$ maxIter
    epsLogLik <- ISCMEBObj@parameterList$epsLogLik
    verbose <- ISCMEBObj@parameterList$verbose
    int.model <- ISCMEBObj@parameterList$int.model
    init.start <- ISCMEBObj@parameterList$init.start
    Sigma_equal <- ISCMEBObj@parameterList$Sigma_equal 
    Sigma_diag <- ISCMEBObj@parameterList$Sigma_diag
    seed <- ISCMEBObj@parameterList$seed
    coreNum <- ISCMEBObj@parameterList$coreNum
    criteria <- ISCMEBObj@parameterList$criteria
    c_penalty <- ISCMEBObj@parameterList$c_penalty
  
    ## Centering
    ISCMEBObj@resList <- SC.MEB2(
        VList, 
        AdjList = ISCMEBObj@AdjList, 
        K=K, 
        beta_grid= beta_grid,
        maxIter_ICM=maxIter_ICM,
        maxIter= maxIter,
        epsLogLik=epsLogLik,
        verbose=verbose,
        int.model=int.model,
        init.start=init.start,
        Sigma_equal=Sigma_equal, 
        Sigma_diag=Sigma_diag,
        seed=seed,
        coreNum=coreNum,
        criteria=criteria,
        c_penalty=c_penalty)
    
    return(ISCMEBObj)
}

