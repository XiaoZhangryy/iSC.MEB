#' @importFrom Matrix colSums
getAdj_auto <- function(pos, lower.med=4, upper.med=6, radius.upper= 100){
  if (!inherits(pos, "matrix"))
    stop("method is only for  matrix object!")
  
  radius.lower <- 1
  Adj_sp <- getneighborhood_fast(pos, radius=radius.upper)
  Med <- summary(Matrix::rowSums(Adj_sp))['Median']
  if(Med < lower.med) stop("The radius.upper is too smaller that cannot find median neighbors greater than 4.")
  start.radius <- 1
  Med <- 0
  message("Find the adjacency matrix by bisection method...")
  maxIter <- 30
  k <- 1
  while(!(Med >= lower.med && Med <=upper.med)){ # ensure that each spot has about 4~6 neighborhoods in median.
    
    Adj_sp <- getneighborhood_fast(pos, radius=start.radius)
    Med <- summary(Matrix::rowSums(Adj_sp))['Median']
    if(Med < lower.med){
      radius.lower <- start.radius
      start.radius <- (radius.lower + radius.upper)/2
    }else if(Med >upper.med){
      radius.upper <- start.radius
      start.radius <- (radius.lower + radius.upper)/2
    }
    message("Current radius is ", round(start.radius, 2))
    message("Median of neighborhoods is ", Med)
    if(k > maxIter) {
      message("Reach the maximum iteration but can not find a proper radius!")
      break;
    }
    k <- k + 1
  }
  
  return(Adj_sp)
}

#' @importFrom purrr discard keep
find_neighbors <- function(pos, platform=c('ST', "Visium")) {
  # require(purrr)
  # require(S4Vectors)
  if (tolower(platform) == "visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (tolower(platform) == "st") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  } else {
    stop(".find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  colnames(pos) <- c("row", "col")
  # pos <- DataFrame(pos) # reduce the dependency on S4Vector
  pos <- as.data.frame(pos)
  spot.positions <- pos
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  
  ## Compute coordinates of each possible spot neighbor
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
  
  ## Select spots that exist at neighbor coordinates
  neighbors <- merge(as.data.frame(neighbor.positions), 
                     as.data.frame(spot.positions), 
                     by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                     suffixes=c(".primary", ".neighbor"),
                     all.x=TRUE)
  
  ## Shift to zero-indexing for C++
  #neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
  
  ## Group neighbor indices by spot 
  ## (sort first for consistency with older implementation)
  neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)
  
  ## Discard neighboring spots without spot data
  ## This can be implemented by eliminating `all.x=TRUE` above, but
  ## this makes it easier to keep empty lists for spots with no neighbors
  ## (as expected by C++ code)
  ## df_j <- map(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  df_j <- lapply(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  
  ## Log number of spots with neighbors
  n_with_neighbors <- length(keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ",
          nrow(pos), " spots.")
  
  n <- length(df_j) 
  
  D <- matrix(0,  nrow = n, ncol = n)
  for (i in 1:n) {
    if(length(df_j[[i]]) != 0)
      D[i, df_j[[i]]] <- 1
  }
  ij <- which(D != 0, arr.ind = T)
  ij
}

#' @importFrom Matrix sparseMatrix
getAdj_reg <- function(pos, platform= "Visium"){
  ij <- find_neighbors(pos, platform)
  # library(Matrix)
  n <- nrow(pos)
  Adj_sp <- Matrix::sparseMatrix(ij[,1], ij[,2], x = 1, dims=c(n, n))
  return(Adj_sp)
}

#' Calculate adjacency matrix form spatial coordinates.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @export
#' @param obj A seuList or iSCMEBObj object or posList. 
#' @param platform A string, specify the platform of the provided data, default as "Other", and only support "ST" and "Visium" platform. When platform is specified as "Other", adjacency matrix will built based on L2 distance. 
#' @param lower.med An option parameter determine the lower bound of median of number of neighborhoods. Should be provided when platform is specified as "Other". 
#' @param upper.med An option parameter determine the upper bound of median of number of neighborhoods. Should be provided when platform is specified as "Other". 
#' @param radius.upper An option parameter determine the search radius upper bound of neighborhoods. Should be provided when platform is specified as "Other". 
#'
#' @details seuList is a \link{list} with Seurat object as component, and each Seurat object includes the raw expression count matrix, spatial coordinates and meta data for each data batch, where the spatial coordinates information must be saved in the metadata of Seurat, named "row" and "col" for each data batch. posList is a list of position information, each element of which can be a data.frame or matrix. When the element of posList is data.frame, those data.frame should contain columns with names "row" and "col". When the element of posList is matrix, those matrices should be two column matrices. 
#'
#' @seealso \code{\link{iSCMEBObj-class}}
#'
#' @importFrom stats median
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' iSCMEBObj_toy2 <- CreateNeighbors(iSCMEBObj_toy, platform = "Visium")
#' ## iSCMEBObj_toy2 <- CreateNeighbors(iSCMEBObj_toy, radius.upper = 10)
CreateNeighbors <- function(obj, platform=c("Other", "ST", "Visium"), lower.med=4, upper.med=6, radius.upper= 100) UseMethod("CreateNeighbors")

#' @return Returns a revised iSCMEBObj object.
#' @rdname CreateNeighbors
#' @method CreateNeighbors iSCMEBObj
#' @export
CreateNeighbors.iSCMEBObj <- function(obj, platform=c("Other", "ST", "Visium"), lower.med=4, upper.med=6, radius.upper= 100) {
    posList <- lapply(obj@seulist, function(x) as.matrix(x@meta.data[,c("row", "col")]))
    platform = match.arg(platform)
    AdjList <- runneighbors(posList, platform=platform, lower.med=lower.med, upper.med=upper.med, radius.upper=radius.upper)
    obj@AdjList <- AdjList

    obj@parameterList$platform = platform
    if (platform == "Other") {
        obj@parameterList$lower.med = lower.med
        obj@parameterList$upper.med = upper.med
        obj@parameterList$radius.upper = radius.upper
    }
    obj
}

#' @return Returns a revised list.
#' @rdname CreateNeighbors
#' @method CreateNeighbors list
#' @export
CreateNeighbors.list <- function(obj, platform=c("Other", "ST", "Visium"), lower.med=4, upper.med=6, radius.upper= 100) {
    if (all(sapply(obj, function(x) inherits(x, "Seurat")))) {
        # Check spatial coordinates for each object.
        flag_spa <- sapply(obj, function(seu) all(c("row", "col") %in% colnames(seu@meta.data)))
        if(any(!flag_spa)) stop("CreateNeighbors: check the argument: obj! Each Seurat object in obj must include the spatial coordinates saved in the meta.data, named 'row' and 'col'!")
        posList <- lapply(obj, function(x) x@meta.data[,c("row", "col")])
    } else if (all(sapply(obj, function(x) inherits(x, "matrix") ))) {
        if (any(sapply(obj, nrow) != 2)) stop("CreateNeighbors: check the argument: obj! Each matrix object in obj must be a two columns matrix, corresponding to 'row' and 'col'!")
        posList <- lapply(obj, function(x) x@meta.data[,c("row", "col")])
    } else if (all(sapply(obj, function(x) inherits(x, "data.frame") ))) {
        flag_spa <- sapply(obj, function(seu) all(c("row", "col") %in% colnames(seu)))
        if(any(!flag_spa)) stop("CreateNeighbors: check the argument: obj! Each data.frame object in obj must include the spatial coordinates, named 'row' and 'col'!")
        posList <- lapply(obj, function(x) x[,c("row", "col")])
    } else {
        stop("CreateNeighbors: check the argument: obj! When apply CreateNeighbors to a list, each component of obj must be a Seurat object or matrix or data.frame.")
    }
    
    platform = match.arg(platform)
    AdjList <- runneighbors(posList, platform=platform, lower.med=lower.med, upper.med=upper.med, radius.upper=radius.upper)
    AdjList
}

runneighbors <- function(posList, platform=c("Other", "ST", "Visium"), lower.med=4, upper.med=6, radius.upper= 100) {
    posList <- lapply(posList, as.matrix)
    if (platform == "ST" | platform == "Visium") {
        AdjList <- lapply(posList, getAdj_reg, platform = platform)
    } else if (platform == "Other") {
        AdjList <- lapply(posList, getAdj_auto, lower.med=lower.med, upper.med=upper.med, radius.upper=radius.upper)
    } else {
        stop("CreateNeighbors: undefined platform!")
    }

    AdjList
}

#' Set model parameters for ISCMEB method
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @export
#' @param SCMEBObj An iSCMEBObj object. 
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
#' @return Returns ISCMEB object. 
#' @seealso \code{\link{iSCMEBObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' iSCMEBObj_toy <- SetModelParameters(iSCMEBObj_toy)
SetModelParameters <- function(SCMEBObj, beta_grid=seq(0, 5, by=0.2), maxIter_ICM=6, maxIter=25, 
    epsLogLik=1e-5, verbose=TRUE, int.model="EEE", init.start=1, 
    Sigma_equal = FALSE, Sigma_diag = TRUE, seed=1, coreNum=1, 
    criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty=1) {
    criteria = match.arg(criteria)
    para_settings <- list(
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
    SCMEBObj@parameterList <- c(SCMEBObj@parameterList, para_settings)
    SCMEBObj
}