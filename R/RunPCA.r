#' @importFrom Matrix colSums
approxPCA <- function(X, q){ ## speed the computation for initial values.
    # require(irlba)
    n <- nrow(X)
    svdX <- irlba::irlba(A =X, nv = q)
    PCs <- svdX$u %*% diag(svdX$d[1:q])
    loadings <- svdX$v
    dX <- PCs %*% t(loadings) - X
    Lam_vec <- colSums(dX^2)/n
    return(list(PCs = PCs, loadings = loadings, Lam_vec = Lam_vec))
}

mat2list <- function(z_int, nvec){
    zList_int <- list()
    istart <- 1
    for(i in 1:length(nvec)){
        zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
        istart <- istart + nvec[i]
    }
    return(zList_int)
}

matlist2mat <- function (XList) {
    r_max <- length(XList)
    X0 <- XList[[1]]
    if (r_max > 1) {
        for (r in 2:r_max) {
            X0 <- rbind(X0, XList[[r]])
        }
    }
    return(X0)
}

#' Run a PCA dimensionality reduction.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @export
#' @param obj A seuList or iSCMEBObj object or matrix list. 
#' @param npcs Total Number of PCs to compute and store (15 by default). 
#' @param pca.method A string, specify which kind of PCA to be used. Supporting "APCA" (the approximate PCA), "PCA" (the classical PCA) and "WPCA" (the weighted PCA). Default method is APCA. 
#' @param reduction.name Dimensional reduction name, pca by default
#'
#' @details seuList is a \link{list} with Seurat object as component, and each Seurat object includes the raw expression count matrix, spatial coordinates and meta data for each data batch, where the spatial coordinates information must be saved in the metadata of Seurat, named "row" and "col" for each data batch. matrix list is a list of log-transformed expression matrix, element of which has same columns. 
#'
#' @seealso \code{\link{iSCMEBObj-class}}
#'
#' @importFrom irlba irlba
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' ## For convenience, we show the iSCMEBObj for perform dimension reduction. 
#' ## Users can use PCA method or WPCA.
#' iSCMEBObj_toy2 <- runPCA(iSCMEBObj_toy)
#' ## seulist <- iSCMEBObj_toy$seulist
#' ## seulist <- runPCA(seulist)
runPCA <- function(obj, npcs = 15, pca.method = c("APCA", "PCA", "WPCA"), reduction.name = "pca") UseMethod("runPCA")

#' @return Returns a revised iSCMEBObj object when input a iSCMEBObj.
#' @rdname runPCA
#' @method runPCA iSCMEBObj
#' @export
runPCA.iSCMEBObj <- function(obj, npcs = 15, pca.method = c("APCA", "PCA", "WPCA"), reduction.name = "pca") {
    pca.method = match.arg(pca.method)
    obj@seulist <- runpca(obj@seulist, npcs=npcs, pca.method=pca.method, reduction.name=reduction.name)
    obj@parameterList$npcs = npcs
    obj@parameterList$pca.method = pca.method
    obj@parameterList$reduction.name = reduction.name
    obj
}

#' @return Returns a revised list when input a seuList or matrix list.
#' @rdname runPCA
#' @method runPCA list
#' @export
runPCA.list <- function(obj, npcs = 15, pca.method = c("APCA", "PCA", "WPCA"), reduction.name = "pca") {
    pca.method = match.arg(pca.method)
    if (all(sapply(obj, function(x) inherits(x, "Seurat")))) {
        obj <- runpca(obj, npcs=npcs, pca.method=pca.method, reduction.name=reduction.name)
        return( obj )
    } else if (all(sapply(obj, function(x) inherits(x, "matrix") | inherits(x, "dgCMatrix") ))) {
        obj <- runpca_mat(obj, npcs=npcs, pca.method=pca.method, reduction.name=reduction.name)
        return( obj )
    } else {
        stop("runPCA: check the argument: obj! When apply runPCA to a list, each component of obj must be a Seurat object or matrix or sparse matrix.")
    }
}

runpca_mat <- function(matlist, npcs = 15, pca.method = c("APCA", "PCA", "WPCA"), reduction.name = "pca") {
    pca.method = match.arg(pca.method)
    Xmat <- lapply(matlist, scale, scale=F) %>% matlist2mat() %>% as.matrix()
    if (pca.method == "APCA") {
        princ <- approxPCA(Xmat, q=npcs)
    } else if (pca.method == "PCA") {
        princ <- wpcaCpp(Xmat, npcs, FALSE)
    } else if (pca.method == "WPCA") {
        princ <- wpcaCpp(Xmat, npcs, TRUE)
    } else {
        stop("CreateiSCMEBObject: undefined PCA method.")
    }

    VList <- mat2list(princ$PCs, nvec=sapply(matlist, nrow))
    VList
}

#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject
#' @importFrom dplyr %>%
runpca <- function(seulist, npcs = 15, pca.method = c("APCA", "PCA", "WPCA"), reduction.name = "pca") {
    pca.method = match.arg(pca.method)
    Xmat <- lapply(seulist, function(x) scale(Matrix::t(x@assays$RNA@data), scale=F)) %>% matlist2mat() %>% as.matrix()
    if (pca.method == "APCA") {
        princ <- approxPCA(Xmat, q=npcs)
    } else if (pca.method == "PCA") {
        princ <- wpcaCpp(Xmat, npcs, FALSE)
    } else if (pca.method == "WPCA") {
        princ <- wpcaCpp(Xmat, npcs, TRUE)
    } else {
        stop("CreateiSCMEBObject: undefined PCA method.")
    }

    VList <- mat2list(princ$PCs, nvec=sapply(seulist, ncol))
    for (i in 1:length(VList)) {
        hZ = VList[[i]]
        rownames(hZ) <- colnames(seulist[[i]])
        colnames(hZ) <- paste0(pca.method, "_", 1:ncol(hZ))
        seulist[[i]]@reductions[[reduction.name]] <- Seurat::CreateDimReducObject(embeddings = hZ)
    }
    seulist
}