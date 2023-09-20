

######
#internal methods
######
#---get GRNs
#' @export
setGeneric("GRNs", function(x, ...) standardGeneric("GRNs"))

#' @export
setMethod("GRNs", "SingleCellExperiment", function(x){x@metadata$metaTF$GRNs})
#' @export
setGeneric("GRNs<-", function(x, ...) standardGeneric("GRNs<-"))
#' @export
setReplaceMethod("GRNs", c("SingleCellExperiment"), function(x, ..., value){
  x@metadata$metaTF$GRNs <- SimpleList(value)
  x
})

#---get GRN
#' @export
setGeneric("GRN", function(x, i, ...) standardGeneric("GRN"))

#' @export
setMethod("GRN", c("SingleCellExperiment","missing"), function(x, i, ...){
  if(length(GRNames(x))==0){
    stop("No available GRNs in your object! Please run inferGRNs() first!")
  }
  x@metadata$metaTF$GRNs[[1]]
  })
#' @export
setMethod("GRN", c("SingleCellExperiment","numeric"), function(x, i, ...){
  if(length(GRNames(x))<i){
    stop("No available GRNs in your object!")
  }
  x@metadata$metaTF$GRNs[[i]]
})
#' @export
setMethod("GRN", c("SingleCellExperiment","character"), function(x, i, ...){
  if(length(GRNames(x))==0){
    stop("No available GRNs in your object!")
  }
  x@metadata$metaTF$GRNs[[i]]
})

#' @export
setGeneric("GRN<-", function(x, i, ..., value) standardGeneric("GRN<-"))
#' @export
setReplaceMethod("GRN", c("SingleCellExperiment","missing"), function(x, i, ..., value){
  x@metadata$metaTF$GRNs[[1]] <- value
  x
})
#' @export
setReplaceMethod("GRN", c("SingleCellExperiment","character"), function(x, i,..., value){
  x@metadata$metaTF$GRNs[[i]] <- value
  x
})

#get names
#' @export
setGeneric("GRNames", function(x, ...) standardGeneric("GRNames"))

#' @export
setMethod("GRNames", "SingleCellExperiment", function(x){names(x@metadata$metaTF$GRNs)})

######
#---get all regulons
######

#' @export
setGeneric("Regulons", function(x, ...) standardGeneric("Regulons"))
#' @export
setMethod("Regulons", "SingleCellExperiment", function(x, ...){
  x@metadata$metaTF$regulons
  })

######
# assign regulons
######
#' @export
setGeneric("Regulons<-", function(x, i, ..., value) standardGeneric("Regulons<-"))

#' @export
setReplaceMethod("Regulons", "SingleCellExperiment", function(x, ..., value){
  if(length(RegulonNames(x))>0){
    i <- length(RegulonNames(x))+1
  }else{
    i <- 1L
  }
  if(!(is.list(value) || is(value,"SimpleList"))){
    stop("The value you input is not a list!")
  }
  x@metadata$metaTF$regulons <- SimpleList(value)
  x
})

######
#---get/set the regulon
######
#===get
#' @export
setGeneric("Regulon", function(x, i, ...) standardGeneric("Regulon"))
#' @export
setMethod("Regulon", c("SingleCellExperiment","missing"), function(x,i){
  x@metadata$metaTF$regulons[[1]]
})
#' @export
setMethod("Regulon", c("SingleCellExperiment","numeric"), function(x,i){
  regulons <- x@metadata$metaTF$regulons
  if(i > length(regulons)){
    warning("Only ",length(regulons)," found in your object, but you input: ",i,". using the first one!")
    i <- 1
  }
  regulons[[i]]
})
#' @export
setMethod("Regulon", c("SingleCellExperiment","character"), function(x,i){
  x@metadata$metaTF$regulons[[i]]
})

#===set
#' @export
setGeneric("Regulon<-", function(x, i, ..., value) standardGeneric("Regulon<-"))
#' @export
setReplaceMethod("Regulon", c("SingleCellExperiment","missing"), function(x, i, ..., value){
  x@metadata$metaTF$regulons[['regulon']] <- SimpleList(value)
  x
})
#' @export
setReplaceMethod("Regulon", c("SingleCellExperiment","numeric"), function(x, i, ..., value){
  x@metadata$metaTF$regulons[[i]] <- SimpleList(value)
  x
})
#' @export
setReplaceMethod("Regulon", c("SingleCellExperiment","character"), function(x, i, ..., value){
  x@metadata$metaTF$regulons[[i]] <- SimpleList(value)
  x
})

######
#get regulon names
######
#' @export
setGeneric("RegulonNames", function(x, ...) standardGeneric("RegulonNames"))
#' @export
setMethod("RegulonNames", "SingleCellExperiment", function(x){names(x@metadata$metaTF$regulons)})


######
#differential results
#####
#' @export
setGeneric("diffRegulon", function(x, ...) standardGeneric("diffRegulon"))
#' @export
setMethod("diffRegulon", c("SingleCellExperiment"), function(x){
  diff_regulon <- altExp(x)@metadata$metaTF$diffRegulon
  if(is.null(diff_regulon)) stop("No differential regulons found in your object, ",
                                 "please run cellTypeSpecific() first!")
  diff_regulon
  })
#' @export
setGeneric("diffRegulon<-", function(x, value) standardGeneric("diffRegulon<-"))
#' @export
setMethod("diffRegulon<-", "SingleCellExperiment", function(x, value){
  altExp(x)@metadata$metaTF[['diffRegulon']] <- value
  x
  })


######
#regulon pathway results
#####
#' @export
setGeneric("RegulonPathway", function(x, ...) standardGeneric("RegulonPathway"))
#' @export
setMethod("RegulonPathway", c("SingleCellExperiment"), function(x){
  altExp(x)@metadata$metaTF$RegulonPathway
})
#' @export
setGeneric("RegulonPathway<-", function(x, value) standardGeneric("RegulonPathway<-"))
#' @export
setMethod("RegulonPathway<-", "SingleCellExperiment", function(x, value){
  altExp(x)@metadata$metaTF[['RegulonPathway']] <- value
  x
})

######
#relevel colData
#####
#---getter
#' @export
setGeneric("levelsColData", function(x, col, ...) standardGeneric("levelsColData"))
#' @export
setMethod("levelsColData", c("SingleCellExperiment", "missing"), function(x, col, ...){
  if(is.factor(colData(x)[[1]])){
    levels(colData(x)[[1]])
  }else{
    levels(as.factor(colData(x)[[1]]))
  }
})
#' @export
setMethod("levelsColData", c("SingleCellExperiment", "numeric"), function(x, col, ...){
  if(col > ncol(colData(x))){
    warning("Only ",ncol(colData(x))," columns in your colData, but you input ",col,
            ". Showing the levels of 1-st column.")
    col <- 1
  }
  if(is.factor(colData(x)[[col]])){
    levels(colData(x)[[col]])
  }else{
    levels(as.factor(colData(x)[[col]]))
  }
})
#' @export
setMethod("levelsColData", c("SingleCellExperiment", "character"), function(x, col, ...){
  if(!(col %in% colnames(colData(x)))){
    warning("The ",col ," column is not found in your colData, check the column names of colData(x)!")
    col <- 1
  }
  if(is.factor(colData(x)[[col]])){
    levels(colData(x)[[col]])
  }else{
    levels(as.factor(colData(x)[[col]]))
  }
})
#---setter
#' @export
setGeneric("levelsColData<-", function(x, col, ..., value) standardGeneric("levelsColData<-"))
#' @export
setMethod("levelsColData<-", c("SingleCellExperiment", "missing"), function(x, col, ..., value){
  if(all(value %in% levelsColData(x,col))){
    colData(x)[[1]] <- factor(colData(x)[[1]], levels = value)
  }else{
    stop("The levels in 1-st column is not all matched with your input!")
  }
  x
})
#' @export
setMethod("levelsColData<-", c("SingleCellExperiment", "numeric"), function(x, col, ..., value){
  if(all(value %in% levelsColData(x,col))){
    colData(x)[[col]] <- factor(colData(x)[[col]], levels = value)
  }else{
    stop("The levels in ",col," column is not all matched with your input!")
  }
  x
})
#' @export
setMethod("levelsColData<-", c("SingleCellExperiment", "character"), function(x, col, ..., value){
  if(all(value %in% levelsColData(x,col))){
    colData(x)[[col]] <- factor(colData(x)[[col]], levels = value)
  }else{
    stop("The levels in ",col," column is not all matched with your input!")
  }
  x
})










