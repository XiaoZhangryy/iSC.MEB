#' Select best iSC.MEB model from candidated models. 
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @export
#' @param obj A iSCMEBObj object or SCMEB_Result_Object object. 
#' @param criteria A string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "MAIC", "BIC" and "AIC" ("MBIC" by default).
#' @param c_penalty An optional positive value, the adjusted constant used in the MBIC criteria (1 by default).
#' @param K An optional number of cluster. When K is not null, the iSC.MEB solution with K clusters will be selected, default is \code{NULL}. 
#'
#' @details SCMEB_Result_Object is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the main function of our package. 
#' 
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}
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

#' Get the identity of iSC.MEB
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @export
#' @param obj An iSCMEBObj object or SCMEB_Result_Object object. 
#'
#' @details SCMEB_Result_Object is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the main function of our package. 
#' 
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}
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
    obj@resList$optSolution$cluster
}

#' @return Returns a list of identity, whose i-th element is the identity of i-th tissue section.
#' @rdname idents
#' @method idents SCMEB_Result_Object
#' @export
idents.SCMEB_Result_Object <- function(obj) {
    obj$optSolution$cluster
}

#' Get the identity of iSC.MEB
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @export
#' @param obj An iSCMEBObj object or SCMEB_Result_Object object. 
#'
#' @details SCMEB_Result_Object is an object that contains all iSC.MEB solution information. It is the output of function \code{fit.iscmeb}, which is the main function of our package. 
#' 
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}
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
    obj@resList$optSolution$hZ
}

#' @return Returns a list of embedding, whose i-th element is the embedding of i-th tissue section.
#' @rdname embeddings
#' @method embeddings iSCSCMEB_Result_ObjectMEBObj
#' @export
embeddings.iSCSCMEB_Result_ObjectMEBObj <- function(obj) {
    obj$optSolution$hZ
}