#' Select best ISCMEB model from candidated models. 
#'
#' @useDynLib SC.MEB2, .registration = TRUE
#' @export
#' @param obj A ISCMEBObj object or SCMEB_Result_Object object. 
#' @param criteria A string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "MAIC", "BIC" and "AIC" ("MBIC" by default).
#' @param c_penalty An optional positive value, the adjusted constant used in the MBIC criteria (1 by default).
#' @param K An optional number of cluster. When K is not null, the ISCMEB solution with K clusters will be selected, default is \code{NULL}. 
#'
#' @details SCMEB_Result_Object is an object that contains all ISCMEB solution information. It is the output of function SC.MEB2, which is the main function of our package. 
#' 
#' @seealso \code{\link{ISCMEBObj-class}}, \code{\link{SC.MEB2}}
#'
#' @examples
#' data(ISCMEBObj_toy)
#' library(Seurat)
#' ## For convenience, we show the ISCMEBObj for perform dimension reduction. 
#' ## Users can use PCA method or WPCA.
#' ISCMEBObj_toy2 <- runPCA(ISCMEBObj_toy)
#' ## seulist <- ISCMEBObj_toy$seulist
#' ## seulist <- runPCA(seulist)
SelectModel <- function(obj, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) UseMethod("SelectModel")

#' @return Returns a revised ISCMEBObj object.
#' @rdname SelectModel
#' @method SelectModel ISCMEBObj
#' @export
SelectModel.ISCMEBObj <- function(obj, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) {
    obj@resList <- selectmodel(obj@resList, criteria=criteria, c_penalty = c_penalty, K = K)
    obj
}

#' @return Returns a revised list
#' @rdname SelectModel
#' @method SelectModel SCMEB_Result_Object
#' @export
SelectModel.SCMEB_Result_Object <- function(obj, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) {
    obj <- selectmodel(obj, criteria=criteria, c_penalty = c_penalty, K = K)
    obj
}

selectmodel <- function(resList, criteria=c("MBIC", "MAIC", "BIC", "AIC"), c_penalty = 1, K = NULL) {
    if (is.null(K)) {
        criteria = match.arg(criteria)
        d <- resList$param$q
        n <- sum(resList$param$nt)
        Kseq <- resList$param$Kseq
        resList$MAIC <- -2.0*resList$log_likes + resList$dfs*2*log(log(d+n))*c_penalty
        resList$MBIC <- -2.0*resList$log_likes + resList$dfs*log(n)*log(log(d+n))*c_penalty
        
        Mycriteria <- switch(
            criteria, 
            MAIC = resList$MAIC,
            AIC  = resList$AIC,
            MBIC = resList$MBIC, 
            BIC  = resList$BIC 
        )
        resList$opt  = which.min(Mycriteria)
        resList$optK = Kseq[resList$opt]
        resList$optSolution = resList$fit[[resList$opt]]
        resList$param$modelselect = switch(
            criteria, 
            MAIC = paste0("MAIC_", c_penalty),
            AIC  = "AIC",
            MBIC = paste0("MBIC_", c_penalty), 
            BIC  = "BIC" 
        )
    } else {
        Kseq <- resList$param$Kseq
        if (!(K %in% Kseq)) stop("K is not in Kseq!")
        resList$optK = K
        resList$opt  = which(Kseq == K)
        resList$optSolution = resList$fit[[resList$opt]]
        resList$param$modelselect = paste0("Given K=", K)
    }
    
    return(resList)
}