setGeneric("inferGRNs", function(x, ...) standardGeneric("inferGRNs"))

#' Inferring Gene Regulatory Network (GRN) from Single-cell Expression Data
#'
#' @param x An object used for GRN inferring . We support object with the
#' following types: \code{\link{matrix}} and an object of
#' \code{\link{SingleCellExperiment}}.
#' @param use_assay Assay used for GRN inferring, can be \code{numeric} or
#' \code{character} with assay names. Default:1,  the first assay.
#' @param method Methods used for weighted GRN inferring. Current version support
#' four types of methods: \code{\link{PUIC}}, \code{\link{GENIE3}},
#' \code{\link{PCOR}} and \code{\link{SINCERITIES}}. Default: PUIC.
#' @param use_regulator Regulators used for weighted GRN inferring, eg. TFs.
#' Default: NULL, all genes in rows.
#' @param use_target Targets used for weighted GRN inferring. Default: NULL, 
#' all genes in rows.
#' @param time_col A numeric \code{vector} or a one column \code{data.frame}
#' contains time points for each sample. If the input object x is a
#' \code{\link{SingleCellExperiment}} with `colData` information, the `time_col`
#' could be the column name with time information in `colData(x)`. The time points
#' must be numeric (eg. 0, 12, 24, 48) and at least 3 time points needed.
#' Only used for `sincerities` method. Default: NULL.
#' @param ncores Number of cores used for parallel calculation. Using
#' \code{\link{parallel::detectCores()}} to detect the available cores on your
#' computer. Default: 1.
#' @param ... Other parameters passed to \code{\link{inferWeightExpMat}}.
#'
#' @rdname inferGRNs
#' @export
setMethod("inferGRNs", "matrix", function(x,
                                          method = c("PUIC", "genie3", "pcor", "sincerities"),
                                          use_regulator = NULL,
                                          use_target = NULL,
                                          time_col = NULL,
                                          ncores=1,
                                          ...){
  exp_mat <- as.matrix(x)
  #---add time columns
  if(!is.null(time_col)){
    if(is.character(time_col)){
      stop("The 'time_col' can not be character when x is matrix!")
    }else if(is.numeric(time_col) && length(time_col)==ncol(x)){
      if(!is.null(names(time_col))){
        col_time <- data.frame(sample=colnames(x),time=time_col[colnames(x)])
      }else{
        col_time <- data.frame(sample=colnames(x),time=time_col)
      }
      row.names(col_time) <- colnames(x)
    }else if(is.data.frame(time_col) && nrow(time_col)==ncol(x)){
      if(length(intersect(colnames(x),row.names(time_col)))==ncol(x)){
        col_time <- data.frame(sample=colnames(x),time=time_col[colnames(x),1])
      }else{
        col_time <- data.frame(sample=colnames(x),time=time_col[,1])
      }
    }else{
      stop("Current version of inferGRNs() only support 3 types of data.")
    }
  }else{
    col_time=NULL
  }
  grn_mat <- inferWeightExpMat(exp_mat = exp_mat,
                               method = method,
                               use_regulator = use_regulator,
                               col_time = col_time, ...)
  grn_mat
})

#' @rdname inferGRNs
#' @export
setMethod("inferGRNs", "SingleCellExperiment", function(x,
                                                        use_assay = 1L,
                                                        method = c("PUIC","genie3", "pcor", "sincerities"),
                                                        use_regulator = NULL,
                                                        use_target = NULL,
                                                        time_col = NULL,
                                                        ncores=1,
                                                        ...){
  exp_mat <- as.matrix(assay(x = x, i = use_assay))
  #---add time columns
  if(!is.null(time_col)){
    if(is.character(time_col) && length(time_col)==1){
      if(time_col %in% colnames(colData(x))){
        col_time <- data.frame(sample=row.names(colData(x)),time=colData(x)[[time_col]])
        row.names(col_time) <- row.names(colData(x))
      }else{
        stop("The 'time_col' is not found in colData(x)!")
      }
    }else if(is.numeric(time_col) && length(time_col)==nrow(colData(x))){
      if(!is.null(names(time_col))){
        colData(x)$time <- time_col[row.names(colData(x))]
      }else{
        colData(x)$time <- time_col
      }
      col_time <- data.frame(sample=row.names(colData(x)),time=colData(x)[["time"]])
      row.names(col_time) <- row.names(colData(x))
    }else if(is.data.frame(time_col) && nrow(time_col)==nrow(colData(x))){
      if(length(intersect(row.names(colData(x)),row.names(time_col)))==nrow(colData(x))){
        colData(x)$time <- time_col[row.names(colData(x)),1]
      }else{
        colData(x)$time <- time_col[,1]
      }
      col_time <- data.frame(sample=row.names(colData(x)),time=colData(x)[["time"]])
      row.names(col_time) <- row.names(colData(x))
    }else{
      stop("Current version of inferGRNs() only support 3 types of data.")
    }
  }else{
    col_time <- NULL
  }
  #---runing the GRN infering
  grn_mat <- inferWeightExpMat(exp_mat = exp_mat,
                               method = method,
                               use_regulator = use_regulator,
                               col_time = col_time,
                               ncores=ncores, ...)
  GRN(x, i = method) <- S4Vectors::DataFrame(grn_mat, check.names = FALSE)
  x
})
###################################################################
#' Inferring weighted expression matrix
#'
#' @param exp_mat Matrix. Numeric matrix of gene expression data with rows were gene names and columns were samples.
#' @param method Methods used for weighted matrix calculation. Current version only
#' support four types of methods. The "genie3" for \code{\link{GENIE3}} method from package \code{\link{GENIE3}}.
#' The "pcor" for \code{\link{ppcor}} package. The "spearman", were common used method in \code{\link{cor()}} function.
#' @param use_regulator Regulators used for weighted matrix calculation, usually TFs. Default: NULL, stands for all
#' genes were used.
#' @param exp_cutoff Numeric. The number or percentage of cells expressed for each gene. Default: 1, at least one cell.
#' @param col_time A two columns data.frame. The first column must be sample names which should be identical
#' with sample names in \code{exp_mat}. The second column must be time points for each sample which should
#' be numeric (eg. 0, 12, 24, 48). If you do not have exactly time information for your cell types, you can
#' just set time points from zero (eg. 0, 1, 2, 3, 4, ...). Note: at least 3 time points needed.
#' @param ncores Number of cores used for parallel calculation. Using \code{\link{parallel::detectCores()}}
#' to detect the cores on your computer.Default: 1.
#' @param verbose Whether to show messages?
#' @param ... Other parameters passed to \pkg{GENIE3} \pkg{sincerities} and \pkg{sclink}.

#'
#' @return A matrix with weighted regulator-target correlation result.
#' @importFrom ppcor pcor
#' @importFrom GENIE3 GENIE3
#' @importFrom scLink sclink_cor sclink_norm
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#' @importFrom PUIC PUIC
#'
#' @export
inferWeightExpMat <- function(exp_mat,
                              method = c("pcor", "PUIC", "genie3", "sincerities"),
                              use_regulator = NULL,
                              exp_cutoff = 0,
                              col_time = NULL,
                              ncores = 1,
                              verbose = interactive(),
                              ...){
  method <- match.arg(method)
  if(verbose) message("Inferring GRNs using method: ", method)
  #------check matrix format
  exp_mat <- as.matrix(exp_mat)
  if(!is.numeric(exp_cutoff)) stop("The 'exp_cutoff' must be numeric between 0 and 1.")
  if(length(exp_cutoff)!=1) stop("The 'exp_cutoff' must be of length 1.")
  if(exp_cutoff == 0){
    exp_mat <- exp_mat[rowSums(exp_mat)>0,]
  }else if(exp_cutoff>0 && exp_cutoff<1){
    exp_mat <- exp_mat[apply(exp_mat,1,function(x){length(x[which(x>0)])>round(length(x)*exp_cutoff)}),]
  }else if(exp_cutoff>=1){
    exp_mat <- exp_mat[apply(exp_mat,1,function(x){length(x[which(x>0)])>round(exp_cutoff)}),]
  }else{
    stop("The 'exp_cutoff' should be an integer represent the number of cells ",
         "expressed or the ratio of expressed cells per gene.")
  }
  if(nrow(exp_mat)==0) stop("No gene names found in your input expression matrix!")
  if(!is.numeric(exp_mat)) stop("The value in your expression matrix is not numeric!")
  if(any(is.na(exp_mat))) stop("NA found in your expression matrix! Please check your input data!")
  #------check regulator
  if(is.null(use_regulator)){
    use_regulator <- row.names(exp_mat)
  }else{
    use_regulator <- intersect(use_regulator, row.names(exp_mat))
    if(length(use_regulator)==0) stop("No regulator detected while the 'use_regulator' is not NULL!")
  }
  #------check time_points
  if(is.null(col_time) && method=="sincerities"){
    stop("The 'col_time' can not be NULL when using 'sincerities' method!")
  }
  ##########
  #------start calculate GRNs
  if(method == "PUIC"){
    cor_mat <- PUIC(expMat = exp_mat, regulators = use_regulator, 
                          ncores = ncores, verbose = verbose, ...)
  }else if(method == "genie3"){
    #---GENIE3
    cor_mat <- suppressWarnings(GENIE3::GENIE3(exprMatrix=exp_mat, regulators=use_regulator, nCores=ncores,...))
  }else if(method == "sincerities"){
    #------sincerities
    sincerities_res <- run_sincerities(exp_mat = exp_mat, col_data = col_time, ...)
    cor_mat <- sincerities_res$adj_matrix
    cor_mat <- cor_mat[use_regulator,,drop=FALSE]
  }else if(method == "pcor"){
    #---PPCOR
    pcor_res <-suppressWarnings(ppcor::pcor(t(exp_mat), method = "spearman"))
    cor_mat <- pcor_res$estimate
    dimnames(cor_mat) <- list(row.names(exp_mat),row.names(exp_mat))
    cor_mat <- cor_mat[use_regulator,,drop=FALSE]
  }else{
    stop("Method should be included in one of the methods.")
  }
  cor_mat
}

######
#------calculating cor with normal method
######
.norm_methods <- function(matrix,
                          use_regulator = NULL,
                          method = c("spearman", "pearson", "kendall"),
                          ncores = 1){
  #method <- match.arg(method)
  if(is.null(use_regulator)){
    use_regulator <- row.names(matrix)
  }
  use_regulator <- intersect(use_regulator, row.names(matrix))
  #------pa
  cl <- parallel::makeCluster(ncores)
  cor_mat <- pbapply::pbsapply(use_regulator,function(gname, matrix){
    exp_i <- matrix[gname,]
    apply(matrix,1,function(x){cor(x,exp_i,method = "spearman")})
  }, matrix=matrix,cl=cl)
  parallel::stopCluster(cl)
  t(cor_mat)
}

######
#------sincerities
######

#' Calculating GRNs with SINCERITIES
#'
#' @param exp_mat A gene by sample expression matrix.
#' @param col_data A two columns data.frame. The first column must be sample names which should be identical
#' with sample names in \code{exp_mat}. The second column must be time points for each sample which should
#' be numeric (eg. 0, 12, 24, 48). If you do not have exactly time information for your cell types, you can
#' just set time points from zero (eg. 0, 1, 2, 3, 4, ...). Note: at least 3 time points needed.
#' @param ... Other parameters passed to \code{SINCERITIES}, \code{SINCERITIES_PLUS} and
#' \code{final_ranked_predictions}.
#'
#' @return A list containing the following information:
#' \itemize{
#' \item adj_matrix: \emph{m} by \emph{m} matrix containing the weights of
#' regulatory edges. The larger \code{adj_matrix[i,j]} indicates higher confidence
#' that the corresponding edge exists (i.e., gene \emph{i} regulating gene \emph{j}).
#' \item final_table: a \code{data.frame} containing the regulation relationships
#' between source genes and targets and its regulatory directions.}
#'
#' @author \code{SINCERITIES} were first created by Nan Papili Gao and R version implemented
#' by Ziyi Hua. Institute for Chemical and Bioengineering.ETH Zurich. E-mail: nanp@ethz.ch.\cr
#' Puwen Tan Integrated the main functions from \code{SINCERITIES} into \code{run_sincerities()}, which
#' was included in \pkg{metaTF} package.
#'
#' @export
run_sincerities <- function(exp_mat, col_data, ...){
  #------check exp_mat
  exp_mat <- as.matrix(exp_mat)
  if(!identical(sort(colnames(exp_mat)),sort(as.character(col_data[,1])))){
    stop("The samples in exp_mat is not identical in col_data!")
  }
  #------check col_data
  if(is.null(col_data)) stop("The 'col_data' can't be NILL!")
  if(!is.data.frame(col_data)) stop("The 'col_data' must be a data.frame!")
  if(!is.numeric(col_data[,2])) stop("The second column of 'col_data' must be numeric!")
  exp_mat <- t(exp_mat)
  gene_num <- ncol(exp_mat)
  gene_names <- colnames(exp_mat)
  col_data <- col_data[order(col_data[,2]),]
  time_line <- col_data[,2]
  time <- sort(unique(time_line))
  exp_mat <- exp_mat[as.character(col_data[,1]),]
  single_cell_data <- by(exp_mat,time_line,identity)
  DATA <- list(time=time,num_time_points=length(time),
               totDATA=exp_mat,timeline=time_line,numGENES=gene_num,
               genes=gene_names,singleCELLdata=single_cell_data)
  #------get gene regulatory network
  if(length(time)>=5){
    sincerities_res <- SINCERITIES(DATA = DATA, ...)
  }else if(length(time)<5 && length(time)>2){
    sincerities_res <- SINCERITIES_PLUS(DATA = DATA, ...)
  }else{
    stop("SINCERITIES method must have at least 3 time points.")
  }
  #---normalization
  adj_matrix <- sincerities_res$adj_matrix/max(sincerities_res$adj_matrix)
  final_table <- final_ranked_predictions(adj_matrix,DATA$genes,SIGN=1, ...)
  list(adj_matrix = adj_matrix, final_table = final_table)
}

