#' Select best iSC.MEB model from candidated models. 
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{SelectModel} is used to select an optimal solution by given criteria or number of clusters.
#' @export
#' @param obj A iSCMEBObj object or iSCMEBResObj object. 
#' @param criteria A string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "MAIC", "BIC" and "AIC" ("MBIC" by default).
#' @param c_penalty An optional positive value, the adjusted constant used in the MBIC criteria (1 by default).
#' @param K An optional number of cluster. When K is not null, the iSC.MEB solution with K clusters will be selected, default is \code{NULL}. 
#'
#' @details iSCMEBResObj is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the body of our algorithm. 
#' 
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}, \code{\link{iSCMEBResObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' iSCMEBObj_toy <- SelectModel(iSCMEBObj_toy)
SelectModel <- function(obj, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) UseMethod("SelectModel")

#' @return Returns a revised iSCMEBObj object.
#' @rdname SelectModel
#' @method SelectModel iSCMEBObj
#' @export
SelectModel.iSCMEBObj <- function(obj, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) {
    obj@resList <- selectmodel(obj@resList, criteria=criteria, c_penalty = c_penalty, K = K)
    obj
}

#' @return Returns a revised list
#' @rdname SelectModel
#' @method SelectModel iSCMEBResObj
#' @export
SelectModel.iSCMEBResObj <- function(obj, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) {
    obj <- selectmodel(obj, criteria=criteria, c_penalty = c_penalty, K = K)
    obj
}

selectmodel <- function(resList, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) {
    if (is.null(K)) {
        criteria = match.arg(criteria)
        d <- resList@paramList$q
        n <- sum(resList@paramList$n)
        Kseq <- resList@paramList$K
        resList@paramList$MAIC <- -2.0*resList@paramList$log_likes + resList@paramList$dfs*2*log(log(d+n))*c_penalty
        resList@paramList$MBIC <- -2.0*resList@paramList$log_likes + resList@paramList$dfs*log(n)*log(log(d+n))*c_penalty
        
        Mycriteria <- switch(
            criteria, 
            MAIC = resList@paramList$MAIC,
            AIC  = resList@paramList$AIC,
            MBIC = resList@paramList$MBIC, 
            BIC  = resList@paramList$BIC 
        )
        resList@paramList$opt  = which.min(Mycriteria)
        resList@paramList$optK = Kseq[resList@paramList$opt]
        resList@reduction$iSCMEB = resList@fitList[[resList@paramList$opt]]$hZ
        resList@idents = lapply(resList@fitList[[resList@paramList$opt]]$cluster, as.vector)
        resList@paramList$modelselect = switch(
            criteria, 
            MAIC = paste0("MAIC_", c_penalty),
            AIC  = "AIC",
            MBIC = paste0("MBIC_", c_penalty), 
            BIC  = "BIC" 
        )
    } else {
        Kseq <- resList@paramList$K
        if (!(K %in% Kseq)) stop("SelectModel: check argument: K! K is not in Kseq!")
        resList@paramList$optK = K
        resList@paramList$opt  = which(Kseq == K)
        resList@reduction$iSCMEB = resList@fitList[[resList@paramList$opt]]$hZ
        resList@idents = lapply(resList@fitList[[resList@paramList$opt]]$cluster, as.vector)
        resList@paramList$modelselect = paste0("Given K=", K)
    }
    
    return(resList)
}

#' Get the identity of iSC.MEB
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{idents} is used to get the clustering result of iSC.MEB model.
#' @export
#' @param obj An iSCMEBObj object or iSCMEBResObj object. 
#'
#' @details iSCMEBResObj is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the main function of our package. 
#' 
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}
#' @details iSCMEBResObj is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the body of our algorithm. 
#'
#' @examples
#' data(iSCMEBObj_toy)
#' IdentsList <- idents(iSCMEBObj_toy)
idents <- function(obj) UseMethod("idents")

#' @return Returns a list of identity, whose i-th element is the identity of i-th tissue section.
#' @rdname idents
#' @method idents iSCMEBObj
#' @export
idents.iSCMEBObj <- function(obj) {
    yList2factor(obj@resList@idents)
}

#' @return Returns a list of identity, whose i-th element is the identity of i-th tissue section.
#' @rdname idents
#' @method idents iSCMEBResObj
#' @export
idents.iSCMEBResObj <- function(obj) {
    yList2factor(obj@idents)
}

#' Get the identity of iSC.MEB
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{embeddings} is used to get the low dimensional embeddings of iSC.MEB model.
#' @export
#' @param obj An iSCMEBObj object or iSCMEBResObj object. 
#'
#' @details iSCMEBResObj is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the main function of our package. 
#' 
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}
#' @details iSCMEBResObj is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the body of our algorithm. 
#'
#' @examples
#' data(iSCMEBObj_toy)
#' EmbeddingsList <- embeddings(iSCMEBObj_toy)
embeddings <- function(obj) UseMethod("embeddings")

#' @return Returns a list of embedding, whose i-th element is the embedding of i-th tissue section.
#' @rdname embeddings
#' @method embeddings iSCMEBObj
#' @export
embeddings.iSCMEBObj <- function(obj) {
    obj@resList@reduction$iSCMEB
}

#' @return Returns a list of embedding, whose i-th element is the embedding of i-th tissue section.
#' @rdname embeddings
#' @method embeddings iSCMEBResObj
#' @export
embeddings.iSCMEBResObj <- function(obj) {
    obj@reduction$iSCMEB
}

getFeature <- function(obj, feature.name) UseMethod("getFeature")

getFeature.iSCMEBObj <- function(obj, feature.name) {
    getFeature(obj@resList, feature.name)
}

getFeature.iSCMEBResObj <- function(obj, feature.name) {
    if (feature.name == "VList") {
        pca.method = obj@paramList$pca.method
        obj@reduction[[pca.method]]
    } else {
        opt = obj@paramList$opt
        obj@fitList[[opt]][[feature.name]]
    }
}
