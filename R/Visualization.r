
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Scatter plot for two-dimensional embeddings
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{plot_scatter} is used to draw a scatter plot.
#' @export
#' @param position A two columns matrix or data.fram, provides the coordinates information.
#' @param labels A vector provides the label information. 
#' @param label_name A string provides the name of labels. 
#' @param axis_names Axis names. 
#' @param cols A vector determine the colors used in the plot.
#' @param no_guides A indicator of whether remove the legend. FALSE by default. 
#' @param do_density A indicator of whether plot the density. FALSE by default. 
#' @param no_axis_name A indicator of whether remove axis names. FALSE by default. 
#' @param point_size The point size of scatter plot. 2 by default. 
#' @param point_alpha The transparency of scatter plot. 1 by default. 
#' @param nrow.legend An integer specify row number in legend. 
#' @param legend.position The position of legend. "bottom" by default. 
#' @param no_axis A indicator of whether remove axis labels and names.
#' @param ... Arguments passed to other methods.
#'
#' @return Returns a scatter plot.
#'
#' @import ggplot2
#' @importFrom grDevices hcl rgb
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' position <- iSCMEBObj_toy@posList[[1]]
#' labels <- idents(iSCMEBObj_toy)[[1]]
#' p1 <- plot_scatter(position, labels, "clusters")
plot_scatter <- function (
    position, labels, label_name, axis_names=c('tSNE1', 'tSNE2'), cols = NULL, 
    no_guides = FALSE, do_density = FALSE, no_axis_name = FALSE, point_size = 2, 
    point_alpha = 1, nrow.legend = NULL, legend.position="bottom", no_axis = FALSE, ...) {
    plt_df <- as.data.frame(position)
    colnames(plt_df) = c("x", "y")
    plt_df[label_name] = labels
    
    if(is.null(cols)){
        cluster <- as.vector(plt_df[[label_name]])
        ngrp <- length(unique(cluster))
        cols    <- gg_color_hue(ngrp)
    }
    
    plt <- plt_df %>% ggplot(aes_string("x", "y", col = label_name, fill = label_name)) + 
        theme(plot.title = element_text(size=15,face = "bold",hjust=0.5),
              axis.text.x = element_text(size=10,face = "bold"),
              axis.title.x = element_text(size=20,face = "bold"),
              text = element_text(size=20),
              axis.text.y = element_text(size=10,face = "bold"),
              axis.title.y = element_text(size=20,face = "bold"),
              legend.position=legend.position,
              legend.key.size = unit(20, "pt"), 
              strip.background = element_rect(linetype = 'solid', color='gray3'),
              strip.text = element_text(size=15,face = "bold"),
              legend.title=element_text(size=20,face = "bold"),
              legend.text=element_text(size=20,face = "bold"),
              panel.background= element_rect(fill = 'white', color='gray')) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 4)), alpha = "none") + 
        geom_point( size = point_size, alpha = point_alpha ) + 
        scale_color_manual(values = cols) + 
        scale_fill_manual(values = cols)
    if (!is.null(nrow.legend))
        plt <- plt + guides(color = guide_legend(nrow = nrow.legend))
    if (no_axis) {
        plt <- plt + theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank()
            )
    } else if (!no_axis_name) {
        plt <- plt + labs(x = axis_names[1], y = axis_names[2])
    }
    if (do_density) 
        plt <- plt + geom_density_2d()
    if (no_guides) 
        plt <- plt + guides(col = 'none', fill = 'none', alpha = 'none')
    if (legend.position == "bottom") {
        plt <- plt + theme(legend.direction = "horizontal")
    } else if (legend.position == "right") {
        plt <- plt + theme(legend.direction = "vertical")
    }
    
    return(plt)
}

#' Plot spatial RGB heatmap.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{plot_RGB} is used to plot spatial RGB heatmap.
#' @export
#' @param position A two columns matrix or data.fram, provides the coordinates information.
#' @param embed_3d A three columns  matrix or data.fram, provides the rgb information.
#' @param pointsize The size of point in the scatter plot.
#' @param ... Arguments passed to other methods.
#'
#' @return Returns a RGB plot.
#'
#' @import ggplot2
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' position <- iSCMEBObj_toy@posList[[1]]
#' tsne3 <- iSCMEBObj_toy@resList@reduction$TSNE3[[1]]
#' p1 <- plot_RGB(position, tsne3)
plot_RGB <- function(position, embed_3d, pointsize=2, ...){
    # suppressMessages(require(ggplot2))
    info = as.data.frame(position)
    colnames(info) = c("sdimx","sdimy")
    
    r = (embed_3d[,1]-min(embed_3d[,1]))/(max(embed_3d[,1])-min(embed_3d[,1]))
    g = (embed_3d[,2]-min(embed_3d[,2]))/(max(embed_3d[,2])-min(embed_3d[,2]))
    b = (embed_3d[,3]-min(embed_3d[,3]))/(max(embed_3d[,3])-min(embed_3d[,3]))
    x =  info$sdimx
    y =  info$sdimy
    dat = data.frame(x,y,r,g,b)
    p1 = ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
        geom_point(size=pointsize) +
        scale_color_identity()+
        theme_void()+
        theme(plot.title = element_text(size=25,face = "bold",hjust=0.5),
              axis.text.x = element_text(size=35,face = "bold"),
              axis.title.x = element_text(size=35,face = "bold"),
              text = element_text(size=20),
              # axis.title.y = element_text(size=30,face = "bold"),
              axis.text.y = element_text(size=35,face = "bold"),
              axis.title.y = element_text(size=35,face = "bold"),
              # axis.text.y = element_blank(),
              legend.position="bottom",
              legend.key.size = unit(50, "pt"), 
              # legend.direction = "horizontal",
              #strip.background = element_blank(),
              # strip.background = element_rect(fill = "grey"),
              strip.text = element_text(size=25,face = "bold"),
              legend.title=element_text(size=35,face = "bold"),
              legend.text=element_text(size=35,face = "bold"))
    p1
}

#' Run UMAP dimensionality reduction
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{CalculateUMAP} is used to run UMAP dimensionality reduction. 
#' @export
#' @param obj A iSCMEBObj object or a iSCMEBResObj object or a list of low dimension embeddings. 
#' @param reduction The name of embeddings to be used to calculate UMAP. If reduction is null, the last added one is used for plotting.
#' @param n_comp An optional positive integer, specify the number of features to be extracted.
#' @param seed A random seed to be used.
#' @param toList An indicator of whether convert the UMAP representation to a list. 
#'
#' @return Returns a UMAP representation.
#'
#' @importFrom scater calculateUMAP
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}, \code{\link{iSCMEBResObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' iSCMEBObj_toy <- CalculateUMAP(iSCMEBObj_toy, reduction="iSCMEB")
#' ## resList <- CalculateUMAP(iSCMEBObj_toy@resList, reduction="iSCMEB")
CalculateUMAP <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE) UseMethod("CalculateUMAP")

#' @return Returns a revised iSCMEBObj object.
#' @rdname CalculateUMAP
#' @method CalculateUMAP iSCMEBObj
#' @export
CalculateUMAP.iSCMEBObj <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE){
    obj@resList <- CalculateUMAP(obj@resList, reduction=reduction, n_comp=n_comp, seed=seed, toList=toList)
    return(obj)

    # set.seed(seed)
    # embeddings = obj@resList$optSolution$hZ
    # nvec = sapply(embeddings, nrow)
    # embeddings = matlist2mat(embeddings)
    # myUMAP <- scater::calculateUMAP(t(embeddings), ncomponents=n_comp)
    # if (toList) myUMAP <- mat2list(myUMAP, nvec)
    # obj@resList[paste0("UMAP", n_comp)] <- myUMAP
    # return(obj)
}

#' @return Returns a revised iSCMEBResObj object.
#' @rdname CalculateUMAP
#' @method CalculateUMAP iSCMEBResObj
#' @export
CalculateUMAP.iSCMEBResObj <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE){
    set.seed(seed)
    if (is.null(reduction)) reduction = names(obj@reduction)[length(obj@reduction)]
    embeddings = obj@reduction[[reduction]]
    nvec = sapply(embeddings, nrow)
    embeddings = matlist2mat(embeddings)
    myUMAP <- scater::calculateUMAP(t(embeddings), ncomponents=n_comp)
    if (toList) myUMAP <- mat2list(myUMAP, nvec)
    obj@reduction[[paste0("UMAP", n_comp)]] <- myUMAP
    return(obj)
}

#' @return Returns the UMAP matrix.
#' @rdname CalculateUMAP
#' @method CalculateUMAP list
#' @export
CalculateUMAP.list <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE){
    set.seed(seed)
    nvec = sapply(obj, nrow)
    embeddings = matlist2mat(obj)
    myUMAP <- scater::calculateUMAP(t(embeddings), ncomponents=n_comp)
    if (toList) myUMAP <- mat2list(myUMAP, nvec)
    return(myUMAP)
}

#' Run t-SNE dimensionality reduction.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{CalculateTSNE} is used to run t-SNE dimensionality reduction. 
#' @export
#' @param obj A iSCMEBObj object or a iSCMEBResObj object or a list of low dimension embeddings. 
#' @param reduction The name of embeddings to be used to calculate t-SNE. If reduction is null, the last added one is used for plotting.
#' @param n_comp An optional positive integer, specify the number of features to be extracted.
#' @param seed A random seed to be used.
#' @param toList An indicator of whether convert the UMAP representation to a list. 
#'
#' @importFrom scater calculateTSNE
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{fit.iscmeb}}, \code{\link{iSCMEBResObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' iSCMEBObj_toy <- CalculateTSNE(iSCMEBObj_toy, reduction="iSCMEB")
#' ## resList <- CalculateTSNE(iSCMEBObj_toy@resList, reduction="iSCMEB")
CalculateTSNE <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE) UseMethod("CalculateTSNE")

#' @return Returns a revised iSCMEBObj object.
#' @rdname CalculateTSNE
#' @method CalculateTSNE iSCMEBObj
#' @export
CalculateTSNE.iSCMEBObj <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE){
    obj@resList <- CalculateTSNE(obj@resList, reduction=reduction, n_comp=n_comp, seed=seed, toList=toList)
    return(obj)

    # set.seed(seed)
    # embeddings = obj@resList$optSolution$hZ
    # nvec = sapply(embeddings, nrow)
    # embeddings = matlist2mat(embeddings)
    # myTSNE <- scater::calculateTSNE(t(embeddings), ncomponents=n_comp)
    # if (toList) myTSNE <- mat2list(myTSNE, nvec)
    # obj@resList[paste0("TSNE", n_comp)] <- myTSNE
    # return(obj)
}

#' @return Returns a revised iSCMEBResObj object.
#' @rdname CalculateTSNE
#' @method CalculateTSNE iSCMEBResObj
#' @export
CalculateTSNE.iSCMEBResObj <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE){
    set.seed(seed)
    if (is.null(reduction)) reduction = names(obj@reduction)[length(obj@reduction)]
    embeddings = obj@reduction[[reduction]]
    nvec = sapply(embeddings, nrow)
    embeddings = matlist2mat(embeddings)
    myTSNE <- scater::calculateTSNE(t(embeddings), ncomponents=n_comp)
    if (toList) myTSNE <- mat2list(myTSNE, nvec)
    obj@reduction[[paste0("TSNE", n_comp)]] <- myTSNE
    return(obj)
}

#' @return Returns the UMAP matrix.
#' @rdname CalculateTSNE
#' @method CalculateTSNE list
#' @export
CalculateTSNE.list <- function(obj, reduction=NULL, n_comp=3, seed=1, toList = TRUE){
    set.seed(seed)
    nvec = sapply(obj, nrow)
    embeddings = matlist2mat(obj)
    myTSNE <- scater::calculateTSNE(t(embeddings), ncomponents=n_comp)
    if (toList) myTSNE <- mat2list(myTSNE, nvec)
    return(myTSNE)
}

#' Add position information to fitted model
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{Addposition} is used to add position information to the fitted model. 
#' @export
#' @param obj A iSCMEBResObj object. 
#' @param posList A list of two columns matrices, each element contains the position information of i-th tissue section. 
#'
#' @seealso \code{\link{fit.iscmeb}}, \code{\link{iSCMEBResObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' posList <- iSCMEBObj_toy@posList
#' resList <- iSCMEBObj_toy@resList
#' resList <- Addposition(resList, posList)
Addposition <- function(obj, posList) {
    nT = obj@paramList$nT
    if (nT != length(posList)) stop("Addposition: the length of posList is not equal to the number of tissue sections!")
    if (any(obj@paramList$n != sapply(posList, nrow))) stop("Addposition: sample size not match!")
    obj@posList = posList
    obj
}


yList2factor <- function(yList0) {
    nT <- length(yList0)
    n <- sapply(yList0, length)
    sampleID <- vector("list", nT)
    yList <- vector("list", nT)
    for (r in 1:nT) {
        yList[[r]] <- yList0[[r]]
        sampleID[[r]] <- rep(r, n[r])
    }
    sampleID <- unlist(sampleID)
    yVec = factor(unlist(yList))

    IntyList <- vector("list", nT)
    istart <- 1
    for (r in 1:nT) {
        IntyList[[r]] <- yVec[istart:(istart+n[r]-1)]
        istart <- istart+n[r]
    }

    return(IntyList)
}

getitem <- function(obj, item) UseMethod("getitem")

getitem.iSCMEBObj <- function(obj, item) {
    getitem(obj@resList, item)
}

getitem.iSCMEBResObj <- function(obj, item) {
    if (item == "cluster") {
        cluster = obj@idents
    } else {
        nT = obj@paramList$nT
        cluster <- vector("list", nT)
        for(i in 1:nT) cluster[[i]] <- rep(i, obj@paramList$n[i])
    }
    yList2factor(cluster)
}

#' Plot spatial heatmap
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{SpaHeatMap} is used to visualize the physical location vs cluster or batch.  
#' @export
#' @param obj A iSCMEBObj object or a iSCMEBResObj object. 
#' @param item Which feature to be used in the plot.
#' @param plot_type A string shows the type of plot.
#' @param title_name A string shows the title of plot.
#' @param combine An indicator of whether plot all on a figure. If TRUE, all figures are plotted; otherwise, return a list with each plot as component. TRUE by default. 
#' @param layout.dim The dimension in the layout of plots when \code{combine = TRUE}.
#' @param common.legend An indicator of whether combine the legend of multiple plots. TRUE by default. 
#' 
#' @param ... Arguments passed to other methods.
#'
#' @importFrom ggpubr ggarrange
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{plot_scatter}}, \code{\link{plot_RGB}}, \code{\link{fit.iscmeb}}, \code{\link{iSCMEBResObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' p1 <- SpaHeatMap(iSCMEBObj_toy, item="cluster", plot_type="Scatter", common.legend=TRUE)
#' p2 <- SpaHeatMap(iSCMEBObj_toy, item="cluster", plot_type="RGB_TSNE", common.legend=TRUE)
SpaHeatMap <- function(
    obj, item=c("cluster", "batch"), plot_type = c("Scatter", "RGB_TSNE", "RGB_UMAP"), title_name = NULL,
    combine = TRUE, layout.dim = c(1, 2), common.legend = TRUE, ...) UseMethod("SpaHeatMap")

#' @return Return a ggplot2 object or list of ggplots objects.
#' @rdname SpaHeatMap
#' @method SpaHeatMap iSCMEBObj
#' @export
SpaHeatMap.iSCMEBObj <- function(
    obj, item=c("cluster", "batch"), plot_type = c("Scatter", "RGB_TSNE", "RGB_UMAP"), title_name = NULL,
    combine = TRUE, layout.dim = c(1, 2), common.legend = TRUE, ...) {

    SpaHeatMap(
        obj@resList, item=item, plot_type=plot_type, title_name=title_name, 
        combine=combine, layout.dim=layout.dim, common.legend=common.legend, ...)
}

#' @return Return a ggplot2 object or list of ggplots objects.
#' @rdname SpaHeatMap
#' @method SpaHeatMap iSCMEBResObj
#' @export
SpaHeatMap.iSCMEBResObj <- function(
    obj, item=c("cluster", "batch"), plot_type = c("Scatter", "RGB_TSNE", "RGB_UMAP"), title_name = NULL,
    combine = TRUE, layout.dim = c(1, 2), common.legend = TRUE, ...) {
    if (is.null(obj@posList)) stop("SpaHeatMap: check argument: obj! posList is not found. Please use Addposition function add position information first!")

    plot_type = match.arg(plot_type)
    sample_name = obj@paramList$sample_name
    nT = obj@paramList$nT
    pList <- vector("list", nT)
    item = match.arg(item)

    if (plot_type == "Scatter") {
        cluster = getitem(obj, item)
        for (i in 1:nT) {
            p1 <- plot_scatter(obj@posList[[i]], cluster[[i]], label_name=item, axis_names=c("row", "col"), ...)
            if(!is.null(title_name)){
                if (length(title_name) == nT) {
                    p1 <- p1 + ggtitle(label=paste0(sample_name[i], ": ", title_name[i]))
                } else {
                    p1 <- p1 + ggtitle(label=paste0(title_name, ": ", sample_name[i]))
                }
            }
            pList[[i]] <- p1
        }
    } else if (plot_type == "RGB_TSNE" | plot_type == "RGB_UMAP") {
        if (plot_type == "RGB_TSNE") {
            if (is.null(obj@reduction$TSNE3)) obj = CalculateTSNE(obj)
            if (is.matrix(obj@reduction$TSNE3)) {
                embed_3d <- mat2list(obj@reduction$TSNE3, obj@paramList$n)
            } else {
                embed_3d <- obj@reduction$TSNE3
            }
        } else {
            if (is.null(obj@reduction$UMAP3)) obj = CalculateUMAP(obj)
            if (is.matrix(obj@reduction$UMAP3)) {
                embed_3d <- mat2list(obj@reduction$UMAP3, obj@paramList$n)
            } else {
                embed_3d <- obj@reduction$UMAP3
            }
        }
        for (i in 1:nT) {
            p1 <- plot_RGB(obj@posList[[i]], embed_3d[[i]], ...)
            if(!is.null(title_name)){
                p1 <- p1 + ggtitle(label=paste0(title_name, sample_name[i]))
            }
            
            pList[[i]] <- p1
        }
    }

    if(combine){
        # pList <- patchwork::wrap_plots(pList, ncol=ncol)
        pList <- ggpubr::ggarrange(plotlist=pList, ncol=layout.dim[2], nrow=layout.dim[1], legend = "bottom", common.legend=common.legend)
    }
    return(pList)
}

#' @return Return a ggplot2 object or list of ggplots objects.
#' @rdname SpaHeatMap
#' @method SpaHeatMap list
#' @export
SpaHeatMap.list <- function(
    obj, item=c("cluster", "batch"), plot_type = c("Scatter", "RGB_TSNE", "RGB_UMAP"), title_name = NULL,
    combine = TRUE, layout.dim = c(1, 2), common.legend = TRUE, ...) {
    if (length(obj) != 2) stop("SpaHeatMap: check argument: obj! The list must constain tow elements. The first is a list of two column matrices to specify the position information, and the second specify the features information!")  
    if (is.list(obj[[1]])) stop("SpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")
    if (any(sapply(obj[[1]], ncol) != 2)) stop("SpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")

    plot_type = match.arg(plot_type)
    nT = length(obj[[1]])
    sample_name = names(obj[[1]])
    if (is.null(sample_name)) sample_name = paste0("Sample", c(1:nT))
    pList <- vector("list", nT)

    if (plot_type == "Scatter") {
        if (is.list(obj[[2]])) stop("SpaHeatMap: check argument: obj! For scatter plot, the second element of list must be a list of vector which specify the features information!")
        if (any(!sapply(obj[[2]], is.vector))) stop("SpaHeatMap: check argument: obj! For scatter plot, the second element of list must be a list of vector which specify the features information!")
        cluster = obj[[2]]
        label_name = names(obj)[2]
        if (is.null(label_name)) label_name = "cluster"
        for (i in 1:nT) {
            p1 <- plot_scatter(obj[[1]][[i]], cluster[[i]], label_name, ...)
            if(!is.null(title_name)){
                p1 <- p1 + ggtitle(label=paste0(title_name, sample_name[i]))
            }
            
            pList[[i]] <- p1
        }
    } else if (plot_type == "RGB_TSNE" | plot_type == "RGB_UMAP") {
        if (is.list(obj[[2]])) stop("SpaHeatMap: check argument: obj! For rgb plot, the second element of list must be a list of three column matrices which specify the features information!")
        if (any(sapply(obj[[1]], ncol) != 3)) stop("SpaHeatMap: check argument: obj! For rgb plot, the second element of list must be a list of three column matrices which specify the features information!")
        if (plot_type == "RGB_TSNE") {
            if (names(obj)[2] != "TSNE") stop("SpaHeatMap: check argument: obj! For RGB_TSNE plot, the second element of list must be named TSNE!")
            embed_3d = obj$TSNE
        } else if (plot_type == "RGB_UMAP") {
            if (names(obj)[2] != "UMAP") stop("SpaHeatMap: check argument: obj! For RGB_TSNE plot, the second element of list must be named UMAP!")
            embed_3d = obj$UMAP
        }
        for (i in 1:nT) {
            p1 <- plot_RGB(obj[[1]][[i]], embed_3d[[i]], ...)
            if(!is.null(title_name)){
                p1 <- p1 + ggtitle(label=paste0(title_name, sample_name[i]))
            }
            
            pList[[i]] <- p1
        }
    }

    if(combine){
        # pList <- patchwork::wrap_plots(pList, ncol=ncol)
        pList <- ggpubr::ggarrange(plotlist=pList, ncol=layout.dim[2], nrow=layout.dim[1], legend = "bottom", ...)
    }
    return(pList)
}




#' Plot low-dimensional embeddings
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{LowEmbedPlot} is used to visualize the low-dimensional embeddings vs cluster or batch to demonstrate the clustering performance and batch remove performance.  
#' @export
#' @param obj A iSCMEBObj object or a iSCMEBResObj object. 
#' @param item Which feature to be used in the plot.
#' @param reduction The name of embeddings to be used in the plot. If reduction is null, the last added one is used for plotting.
#' @param combine An indicator of whether plot all on a figure. If TRUE, all figures are plotted; otherwise, return a list with each plot as component. TRUE by default. 
#' @param cols A vector determine the colors used in the plot.
#' @param layout.dim The dimension in the layout of plots when \code{combine = TRUE}.
#' @param common.legend An indicator of whether combine the legend of multiple plots. TRUE by default. 
#' @param ... Arguments passed to other methods.
#'
#' @importFrom ggpubr ggarrange
#' @seealso \code{\link{iSCMEBObj-class}}, \code{\link{plot_scatter}}, \code{\link{plot_RGB}}, \code{\link{fit.iscmeb}}, \code{\link{iSCMEBResObj-class}}
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' p1 <- LowEmbedPlot(iSCMEBObj_toy, item="cluster", reduction="TSNE2")
#' p2 <- LowEmbedPlot(iSCMEBObj_toy, item="batch", reduction="TSNE2")
LowEmbedPlot <- function(obj, item=c("cluster", "batch"), reduction=NULL, combine=TRUE, cols=NULL, layout.dim=c(1,2), common.legend = TRUE, ...) UseMethod("LowEmbedPlot")

#' @return Return a ggplot2 object.
#' @rdname LowEmbedPlot
#' @method LowEmbedPlot iSCMEBObj
#' @export
LowEmbedPlot.iSCMEBObj <- function(obj, item=c("cluster", "batch"), reduction=NULL, combine=TRUE, cols=NULL, layout.dim=c(1,2), common.legend = TRUE, ...) {
    LowEmbedPlot(obj@resList, item=item, reduction=reduction, combine=combine, cols=cols, layout.dim=layout.dim, common.legend=common.legend, ...)
}

#' @return Return a ggplot2 object.
#' @rdname LowEmbedPlot
#' @method LowEmbedPlot iSCMEBResObj
#' @export
LowEmbedPlot.iSCMEBResObj <- function(obj, item=c("cluster", "batch"), reduction=NULL, combine=TRUE, cols=NULL, layout.dim=c(1,2), common.legend = TRUE, ...) {
    if (is.null(reduction)) reduction = names(obj@reduction)[length(obj@reduction)]
    item = match.arg(item)
    cluster = getitem(obj, item)
    
    if(is.null(cols)){
        sort_id <- sort(unique(as.numeric(unlist(cluster))))
        cols <- gg_color_hue(length(sort_id))[sort_id]
    }
    if(!is.vector(cols)) stop("dimPlot: check argument: cols! it must be a vector object.")
    
    embed_use <- obj@reduction[[reduction]]
    if (combine) {
        cluster = unlist(cluster)
        embed_use = matlist2mat(embed_use)[,1:2]
        p1 <- plot_scatter(embed_use, cluster, label_name=item, cols=cols, ...)
        p1 <- p1 + mytheme_graybox(base_size=16, base_family="serif", bg_fill ="white",border_color = "gray10")
    } else {
        nT = obj@paramList$nT
        pList <- vector("list", nT)
        for (i in 1:nT) {
            pList[[i]] <- plot_scatter(embed_use[[i]][,1:2], cluster[[i]], label_name=item, cols=cols, ...)
        }
        p1 <- ggpubr::ggarrange(plotlist=pList, ncol=layout.dim[2], nrow=layout.dim[1], legend = "bottom", common.legend=common.legend)
        p1 <- p1 + mytheme_graybox(base_size=16, base_family="serif", bg_fill ="white",border_color = "gray10")
    }
 
    return(p1)
}

mytheme_graybox <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                             base_rect_size = base_size/22, border_color = 'gray10', bg_fill='white')
{
    half_line <- base_size/2
    t <- theme(
        panel.background = element_rect(fill = bg_fill, colour = NA), 
        panel.border = element_rect(fill = NA, colour = border_color),
        #line = element_blank(), #rect = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(
            family = base_family, face = "plain", colour = "black", size = base_size, 
            lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = FALSE), 
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
        axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
        axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
        axis.ticks.length.y.right = NULL, legend.box = NULL,
        legend.key.size = unit(1.2, "lines"), legend.position = "right",
        legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
        strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,"pt"), 
        strip.switch.pad.wrap = unit(half_line/2, "pt"), 
        panel.ontop = FALSE,
        panel.spacing = unit(half_line/2, "pt"), plot.margin = unit(rep(0.2,4), "lines"),
        plot.title = element_text(size = rel(1.2), hjust = 0, vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
        plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
        plot.caption = element_text(size = rel(0.8), hjust = 1, vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
        plot.tag = element_text(size = rel(1.2), hjust = 0.5, vjust = 0.5), 
        plot.tag.position = "topleft", complete = TRUE)
    #ggplot2:::ggplot_global$theme_all_null %+replace% t
    t
}

#' SelectKPlot.
#' 
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description
#' The function \code{SelectKPlot} is used to demonstrate the scatter plot of the criteria used vs K for selecting the best K.  
#'
#' @param obj An iSCMEBObj object or a iSCMEBResObj object. 
#' @param criteria A string, specify the criteria used for selecting the number of clusters, supporting "InUse", "MBIC", "MAIC", "BIC" and "AIC" ("InUse" by default).
#' @export
SelectKPlot <- function(obj, criteria=c("InUse", "MBIC", "MAIC", "BIC", "AIC")) UseMethod("SelectKPlot")

#' @return Return a ggplot2 object.
#' @rdname SelectKPlot
#' @method SelectKPlot iSCMEBObj
#' @export
SelectKPlot.iSCMEBObj <- function(obj, criteria=c("InUse", "MBIC", "MAIC", "BIC", "AIC")) {
    SelectKPlot(obj@resList, criteria=criteria)
}

#' @return Return a ggplot2 object.
#' @rdname SelectKPlot
#' @method SelectKPlot iSCMEBResObj
#' @export
SelectKPlot.iSCMEBResObj <- function(obj, criteria=c("InUse", "MBIC", "MAIC", "BIC", "AIC")) {
    if (criteria == "InUse") {
        criteria = obj@paramList$modelselect
        if (grepl( "Given", criteria, fixed = TRUE)) stop("selectKPlot: the optimal solution is selected by given K!")
        if (grepl( "MBIC", criteria, fixed = TRUE)) {
            mbic = obj@paramList$MBIC
            mycriteria = "MBIC"
        } else if (grepl( "MAIC", criteria, fixed = TRUE)) {
            mbic = obj@paramList$MAIC
            mycriteria = "MAIC"
        } else if (grepl( "BIC", criteria, fixed = TRUE)) {
            mbic = obj@paramList$BIC
            mycriteria = "BIC"
        } else if (grepl( "AIC", criteria, fixed = TRUE)) {
            mbic = obj@paramList$AIC
            mycriteria = "AIC"
        } else {
            stop("selectKPlot: unknown criteria!")
        }
    } else {
        mbic <- switch(
            criteria, 
            MAIC = obj@paramList$MBIC,
            AIC  = obj@paramList$AIC,
            MBIC = obj@paramList$MBIC, 
            BIC  = obj@paramList$BIC 
        )
        mycriteria = criteria
    }
    
    Kseq = obj@paramList$K

    data = data.frame(Kseq=Kseq, mbic=mbic)   
    p1 <- ggplot(data, aes(x=Kseq, y=mbic)) + 
        geom_line(size=0.5) + 
        geom_point(size=2) + 
        theme_classic() + 
        labs(x = "K", y = mycriteria) + 
        scale_x_continuous( breaks=Kseq) + 
        theme(plot.title = element_text(size=15,face = "bold",hjust=0.5),
              axis.text.x = element_text(size=10),
              axis.title.x = element_text(size=20),
              text = element_text(size=20),
              axis.text.y = element_text(size=10),
              axis.title.y = element_text(size=20),
              strip.background = element_rect(linetype = 'solid', color='gray3'),
              strip.text = element_text(size=15,face = "bold"),
              legend.title=element_text(size=20,face = "bold"),
              legend.text=element_text(size=20,face = "bold"),
              panel.background= element_rect(fill = 'white', color='gray'))
    p1
}

