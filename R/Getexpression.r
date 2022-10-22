firstup <- function(x) {
    ## First letter use upper capital
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

#' @importFrom pbapply pblapply
#' @importFrom stats coef lm
get_correct_exp <- function(XList, RfList, housekeep, q_unwanted=10){
    if(!all(sapply(XList, is.matrix))){
        XList <- lapply(XList, as.matrix)
    }
    MList <- pbapply::pblapply(XList, function(x) {
        wpcaCpp(scale(x[,housekeep], scale=F), q_unwanted, FALSE)$PCs
    })
    M0 <- matlist2mat(MList)
    rm(MList)
    Rf <- matlist2mat(RfList)
    rm(RfList)
    XList <- lapply(XList, scale, scale=FALSE)
    X0 <- matlist2mat(XList)
    rm(XList)
    lm1 <- lm(X0~ 0+cbind(Rf, M0))
    hK <- ncol(Rf)
    coefmat <- coef(lm1)[-c(1:hK),]
    rm(lm1)
    hX <- X0 - M0 %*% coefmat
    return(hX)
}

get_correct_mean_exp <- function(XList, VList, UList) {
    X0 <- lapply(XList, as.matrix) %>% lapply(scale, scale=FALSE) %>% matlist2mat()
    rm(XList)
    V0 <- matlist2mat(VList)
    rm(VList)
    U0 <- matlist2mat(UList)
    rm(UList)
    W = solve(t(V0) %*% V0) %*% t(V0) %*% X0 # r by p
    X0 - U0 %*% W
}

#' Integrate multiple SRT expression data
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{IntegrateSpaData} is used to integrate multiple SRT data based on the iSCMEBObj by iSC.MEB model fitting.
#' @export
#' @param obj A iSCMEBObj object. 
#' @param species An optional string, one of 'Human', 'Mouse' and 'Unknown', specify the species of the SRT data to help choose the housekeeping genes. 'Unknown' means only using the PRECAST results reconstruct the alligned gene expression.
#' @param custom_housekeep User-specified housekeeping genes.
#'
#' @seealso \code{\link{iSCMEBObj-class}}
#' @importFrom Seurat "Idents<-"
IntegrateSpaData <- function(obj, species=c("human", "mouse", "unknown"), custom_housekeep=NULL) {
    # suppressMessages(require(Matrix))
    # suppressMessages(require(Seurat))
    
    if(!inherits(obj, "iSCMEBObj")) 
        stop("IntegrateSpaData: Check the argument: obj! obj must be a iSCMEBObj object.")
    if(is.null(obj@seulist)) 
        stop("IntegrateSpaData: Check the argument: obj! The slot seulist in iSCMEBObj is NULL!")
    species = match.arg(species)
    
    XList <- lapply(obj@seulist, function(x) Matrix::t(x[["RNA"]]@data))
    n_r <- length(XList)
    for (r in 1:n_r) colnames(XList[[r]]) <- firstup(colnames(XList[[r]]))
    barcodes_all <- lapply(XList, row.names)
    if (any(duplicated(unlist(barcodes_all)))) {
        for (r in 1:n_r) row.names(XList[[r]]) <- paste0(row.names(XList[[r]]), r)
    }
    genelist <- colnames(XList[[1]])
    lower_species <- tolower(species) 
    houseKeep <- switch (lower_species,
        human = {
            # data(Human_HK_genes)
             intersect((genelist), Mouse_HK_genes$Gene)
            },
        mouse={
            #data(Mouse_HK_genes)
            intersect((genelist), Mouse_HK_genes$Gene)
        },
        unknown={
            character()
        }
    )
    houseKeep <- c(houseKeep, custom_housekeep)
    if(length(houseKeep) < 5){
        message("Using only iSC.MEB results to obtain the batch corrected gene expressions since species is unknown or the genelist in obj has less than 5 overlapp with the housekeeping genes of given species.")
        message("Users can specify the custom_housekeep by themselves to use the housekeeping genes based methods.")
        hX <- get_correct_mean_exp(XList, getFeature(obj, "VList"), getFeature(obj, "hV"))
    }else{
        message("Using bouth housekeeping gene and iSC.MEB results to obtain the batch corrected gene expressions.")
        hX <- get_correct_exp(XList, getFeature(obj, "Rf"), houseKeep, q_unwanted=min(10, length(houseKeep)))
    }
    meta_data <- data.frame(batch=unlist(getitem(obj, "batch")), cluster=unlist(getitem(obj, "cluster")))
    row.names(meta_data) <- row.names(hX)
    rm(XList)
    count <- sparseMatrix(i=1,j=1, x=0, dims=dim(t(hX)))
    row.names(count) <- colnames(hX)
    colnames(count) <- row.names(hX)
    seuInt <- CreateSeuratObject(counts = count, assay = 'iSC_MEB', meta.data=meta_data)
    seuInt[['iSC_MEB']]@data <- t(hX)
    
    Idents(seuInt) <- factor(meta_data$cluster)
    return(seuInt)
}

# utils::globalVariables(c("Human_HK_genes", "Mouse_HK_genes", "makeCluster", "stopCluster"))
utils::globalVariables(c("Human_HK_genes", "Mouse_HK_genes"))

#' Run DEG for a Seurat object.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{doDEG} is used to run DEG for a Seurat object
#' @export
#' @param seu A Seurat object. 
#' @param topn The number of features to be extracted when \code{features} is not provided. 10 by default.
#'
#' @seealso \code{\link{IntegrateSpaData}}
#' @importFrom Seurat FindAllMarkers
#' @importFrom dplyr group_by top_n
doDEG <- function(seu, topn=10) {
    # top n DEG
    top10 <- FindAllMarkers(seu) %>% group_by(cluster) %>% top_n(n = topn, wt = avg_log2FC)
    top10
}

#' Heatmap for spots-by-feature matrix
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{doHeatmap} is used to plot heatmap for a Seurat object with expressioin data.
#' @export
#' @param seu A Seurat object. 
#' @param features The features to be plotted.
#' @param topn The number of features to be extracted when \code{features} is not provided. 10 by default.
#' @param cell_label A string specify the name of legend. \code{"Domain"} by default. 
#' @param grp_label An indicator of whether display the group names. \code{FALSE} by default. 
#' @param pt_size The point size used in the plot. 4 by default. 
#' @param grp_color The colors to use for the group color bar. \code{NULL} by default. 
#' @param ncol.legend The number of columns of legend. 
#' @param ... Arguments passed to other methods.
#'
#' @seealso \code{\link{IntegrateSpaData}}
#' @importFrom Seurat DoHeatmap ScaleData
doHeatmap <- function(seu, features=NULL, topn=10, cell_label='Domain', grp_label = FALSE, pt_size=4, grp_color=NULL, ncol.legend=1, ...) {
    if (is.null(features)) features = doDEG(seu, topn)$gene
    seu <- ScaleData(seu)
    seu <- subset(seu, downsample = 400)

    ngrp <- nlevels(factor(seu@meta.data$cluster))
    if(is.null(grp_color)){
        grp_color <- gg_color_hue(ngrp)
    }
    
    Seurat::DoHeatmap(object = seu, features=features, group.colors = grp_color[1:ngrp], label = grp_label, ...) +
        guides(
            color = guide_legend(
                title=cell_label, 
                ncol = ncol.legend,
                override.aes = list(stroke = 1, alpha = 1, shape = 16, size = pt_size, color=grp_color[1:ngrp])
            ), alpha =  "none") + 
        theme(
            legend.text = element_text(size = 10), 
            legend.title = element_text(size = 13, face = "bold"), 
            axis.text.y = element_text(size = 5, face = "italic", family = "serif")
        )
}






