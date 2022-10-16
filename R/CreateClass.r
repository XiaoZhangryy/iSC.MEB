filter_spot <- function(seu, min_feature=0, assay=NULL){ # each spots at least include 1 non-zero features
    if(is.null(assay)) assay <- DefaultAssay(seu)
    col_name <- paste0("nFeature_",assay)
    idx <- seu@meta.data[,col_name] > min_feature
    seu[, idx]
    # subset(seu, subset = nFeature_RNA > min_feature)
}

#' @importFrom Matrix rowSums
filter_gene <- function(seu, min_spots=20, assay= NULL){
    if(is.null(assay)) assay <- DefaultAssay(seu)
    if(sum(dim(seu[[assay]]@counts))!=0){
        gene_flag <- Matrix::rowSums(seu[[assay]]@counts>0)>min_spots
        return(seu[names(gene_flag), ])
    }else if(sum(dim(seu[[assay]]@data))!=0){
        gene_flag <- Matrix::rowSums(seu[[assay]]@data>0)>min_spots
        return(seu[names(gene_flag), ])
    }else{
        stop("filter_gene: Seuat object must provide slots count or data in assay!")
    }
}

#' @importFrom stats sd
#' @importFrom pbapply pbapply
selectIntFeatures <- function(seulist, spaFeatureList, IntFeatures=2000){
  ## This function is used for selecting common informative features
  if(length(seulist) != length(spaFeatureList)) stop("The length of suelist and spaFeatureList must be equal!")
  if(length(seulist) ==1){
    if(length(spaFeatureList[[1]]) >= IntFeatures){
      genelist <- spaFeatureList[[1]][1:IntFeatures]
    }else{
      genelist <- spaFeatureList[[1]]
      warning("The IntFeature is larger than the  number of elements in FeatureList!")
    }
    return(genelist)
  } 
  if(any(sapply(spaFeatureList, length)< IntFeatures))
    stop("Feature list exists number of features less than IntFeatures!")
  geneUnion <- unique(unlist(lapply(spaFeatureList, function(x) x[1:IntFeatures])))
  ## ensure each seuobject has the genes in geneUnion
  gene_delete <- unique(unlist(lapply(seulist, function(x) setdiff(geneUnion, row.names(x)))))
  geneUnion <- setdiff(geneUnion, gene_delete)
  
  
  # Remove zero-variance genes
  genes_zeroVar <- unique(unlist(lapply(seulist, function(x) 
    geneUnion[pbapply::pbapply(x@assays$RNA@counts[geneUnion,],1, sd)==0])))
  gene_Var <- setdiff(geneUnion, genes_zeroVar)
  
  # sort by number of datasets that identified this gene as SVG.
  nsample <- length(seulist)
  numVec <- rep(0, length(gene_Var))
  rankMat <-matrix(NA,length(gene_Var), nsample)
  row.names(rankMat) <- gene_Var
  for(i in 1:length(gene_Var)){
    for(j in 1:nsample){
      if(is.element(gene_Var[i], spaFeatureList[[j]])){
        numVec[i] <- numVec[i] +1
        rank1 <- which(spaFeatureList[[j]]==gene_Var[i])
        rankMat[i, j] <- rank1
      }
    }
    
  }
  
  cutNum <- sort(numVec, decreasing = T)[min(IntFeatures, length(numVec))]
  if(max(numVec)> cutNum){
    genelist1 <- gene_Var[numVec>cutNum]
  }else{
    genelist1 <- NULL
  }
  num_rest_genes <- min(IntFeatures, length(numVec)) - length(genelist1)
  
  gene2 <- gene_Var[numVec==cutNum]
  ### select top 2000 genes that rank 
  rankMat2 <- rankMat[gene2, ]
  rowMedian <- function(xmat, na.rm=TRUE){
    apply(xmat, 1, median, na.rm=na.rm)
  }
  genes1rank <- gene2[order(rowMedian(rankMat2, na.rm=T))[1:num_rest_genes]]
  genelist <- c(genelist1, genes1rank)
  
  return(genelist)
}

#' Each iSCMEBObj object has a number of slots which store information.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @description Each iSCMEBObj object has a number of slots which store information. Key slots to access are listed below.
#' \itemize{
#'   \item \code{seuList} - A list with Seurat object as component, representing the raw expression count matrix, spatial coordinates and meta data for each data batch, where the spatial coordinates information is saved in the metadata of Seurat, named "row" and "col" for eahc data batch.
#'   \item \code{seulist} - A Seurat list after the preprocessing step in preparation for ISCMEB model.
#'   \item \code{AdjList} - The adjacency matrix list for a iSCMEBObj object.
#'   \item \code{parameterList} - The model parameter settings for a iSCMEBObj object.
#'   \item \code{resList} - The results after fitting ISCMEB models.
#'   \item \code{project} - Name of the project.
#' }
setClass("iSCMEBObj", slots=list(
  seuList = "ANY",
  seulist = "ANY",
  AdjList = "ANY", 
  parameterList= "list", 
  resList = "ANY",
  project = "character"
) )

#' Create the ISCMEB object with preprocessing step.
#'
#' @useDynLib iSC.MEB, .registration = TRUE
#' @export
#' @param seuList A list consisting of Seurat objects, where each object is a SRT data batch. The default assay of each Seurat object will be used for data preprocessing and followed model fitting. The specified format about seuList argument can be referred to the details and example. 
#' @param project An optional string, name of the project, default as "ISCMEB".
#' @param gene.number An optional integer, the number of top spatially variable genes (SVGs) or highly variable  genes (HVGs) to be chosen.
#' @param selectGenesMethod An optional integer, the method to select genes for each sample. It supports 'SPARK-X' and 'HVGs' to select genes now. Users can provide self-selected genes using customGenelist argument.
#' @param numCores_sparkx An optional integer, specify the number of CPU cores in SPARK package to use when selecting spatial genes.
#' @param customGenelist An optional string vector, the list of user specified genes to be used for PRECAST model fitting. If this argument is given, SVGs/HVGs will not be selected.
#' @param premin.spots An optional integer, the features (genes) are retained in raw data filtering step with at least premin.spots number of spots, default is 20.
#' @param premin.features An optional integer, the locations are retained in raw data filtering step with at least premin.features number of  nonzero-count features (genes), default is 20.
#' @param postmin.spots An optional integer, the features (genes) are retained in filtering step after common genes selected among all data batches with at least premin.spots number of spots, default is 15.
#' @param postmin.features An optional integer, the locations are retained in filtering step after common genes selected among all data batches  with at least premin.features number of  nonzero-count features (genes), default is 15.
#' @param rawData.preserve An optional logical value, whether preserve the raw seuList data.
#' @param verbose An indictor of whether display the message in the creating process.
#'
#' @details seuList is a \link{list} with Seurat object as component, and each Seurat object includes the raw expression count matrix, spatial coordinates and meta data for each data batch, where the spatial coordinates information must be saved in the metadata of Seurat, named "row" and "col" for each data batch.
#'
#' @return Returns ISCMEB object prepared for ISCMEB model fitting.
#' @seealso \code{\link{iSCMEBObj-class}}
#'
#' @importFrom pbapply pbapply
#' @importFrom DR.SC FindSVGs topSVGs
#' @importFrom methods new
#' @importFrom Seurat NormalizeData CreateSeuratObject FindVariableFeatures DefaultAssay
#' @import gtools
#'
#' @examples
#' data(iSCMEBObj_toy)
#' library(Seurat)
#' seuList <- iSCMEBObj_toy@seulist
#' ## Check the input of seuList for create ISCMEB object.
#' ## Check the default assay for each data batch
#' sapply(seuList, DefaultAssay)
#' ## Check the spatial coordinates in the meta data named "row" and "col".
#' colnames(seuList[[1]]@meta.data)
#' ## Then create ISCMEB object using this seuList. 
#' ## For convenience, we show the  user-specified genes' list for creating ISCMEB object. 
#' ## Users can use SVGs from SPARK-X or HVGs.
#' Genelist <- row.names(seuList[[1]])
#' iSCMEBObj_toy2 <- CreateiSCMEBObject(seuList, customGenelist=Genelist, verbose=FALSE)
CreateiSCMEBObject <- function(seuList, project = "ISCMEB", gene.number=2000, 
    selectGenesMethod='SPARK-X', numCores_sparkx=1, customGenelist=NULL, 
    premin.spots = 20, premin.features=20, postmin.spots=15, postmin.features=15,
    rawData.preserve=FALSE, verbose=TRUE ){
    if(!inherits(seuList, "list")) stop("CreateiSCMEBObject: check the argument: seuList! it must be a list.")
    
    # Check Seurat object
    flag <- sapply(seuList, function(x) !inherits(x, "Seurat"))
    if(any(flag)) stop("CreateiSCMEBObject: check the argument: seuList! Each component of seuList must be a Seurat object.")
    
    # Check spatial coordinates for each object.
    exist_spatial_coods <- function(seu){
        flag_spatial <- all(c("row", "col") %in% colnames(seu@meta.data))
        return(flag_spatial)
    }
    flag_spa <- sapply(seuList,    exist_spatial_coods)
    if(any(!flag_spa)) stop("CreateiSCMEBObject: check the argument: seuList! Each Seurat object in seuList must include the spatial coordinates saved in the meta.data, named 'row' and 'col'!")
    
    # Check cores 
    if(numCores_sparkx<0) 
        stop("CreateiSCMEBObject: check the argument: numCores_sparkx! It must be a positive integer.")
 
    # Check customGenelist
    if(!is.null(customGenelist) && (!is.character(customGenelist))) 
        stop("CreateiSCMEBObject: check the argument: customGenelist! It must be NULL or a character vector.")
    
    ## inheriting
    object <- new(
        Class = "iSCMEBObj",
        seuList = seuList,
        seulist = NULL,
        AdjList = NULL, 
        parameterList= list(),
        resList = NULL,
        project = project
    )
    seuList <- object@seuList 
    if(verbose)
        message("Filter spots and features from Raw count data...")
    seuList <- lapply(seuList, filter_spot, premin.features)
    seuList <- pbapply::pblapply(seuList, filter_gene, premin.spots)
    if(verbose) message(" \n ")
    
    if(is.null(customGenelist)){
        if(tolower(selectGenesMethod)=='spark-x'){
            seuList <- pbapply::pblapply(seuList, DR.SC::FindSVGs, nfeatures=gene.number, 
                                                                     num_core=numCores_sparkx, verbose=verbose)
            spaFeatureList <- lapply(seuList, DR.SC::topSVGs, ntop = gene.number)
        }else if(tolower(selectGenesMethod)=='hvgs'){
            seuList <- pbapply::pblapply(seuList ,FindVariableFeatures, nfeatures=gene.number, verbose=verbose)
            getHVGs <- function(seu){
                assay <- DefaultAssay(seu)
                seu[[assay]]@var.features
            }
            spaFeatureList <- lapply(seuList, getHVGs)
        }else{
            stop("CreateiSCMEBObject: check the argument: selectGenesMethod! It only support 'SPARK-X' and 'HVGs' to select genes now. You can provide self-selected genes using customGenelist argument.")
        }
        
        spaFeatureList <- lapply(spaFeatureList, function(x) x[!is.na(x)])
        if(any(sapply(spaFeatureList, length)< gene.number)){
            gene.number_old <- gene.number
            gene.number <- min(sapply(spaFeatureList, length))
            warning(paste0("Number of genes in one of sample is less than ", gene.number_old, ", so set minimum number of SVGs as gene.number=", gene.number) )
        }
        if(verbose)
            message("Select common top variable genes    for multiple samples...")
        
        genelist <- selectIntFeatures(seuList, spaFeatureList=spaFeatureList, IntFeatures=gene.number)
    }else{
        genelist <- customGenelist
        geneNames <- Reduce(intersect,(lapply(seuList, row.names))) # intersection of    genes from each sample
        if(any(!(customGenelist %in% geneNames)))
            stop("CreateiSCMEBObject: check the argument: customGenelist! It contains the gene not in seuList.")
    }
    
    seulist <- lapply(seuList, function(x) x[genelist, ])
    if(verbose)
         message("Filter spots and features from SVGs(HVGs) count data...")
    seulist <- lapply(seulist, filter_spot, postmin.features)
    seulist <- pbapply::pblapply(seulist, filter_gene, postmin.spots)
    seulist <- lapply(seulist, NormalizeData, verbose=verbose)
    object@seulist <- seulist
    
    if(!rawData.preserve){
        object@seuList <- NULL
    }
    return(object)
}