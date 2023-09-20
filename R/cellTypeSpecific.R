setGeneric("findCellTypeSignature", function(x, ...) standardGeneric("findCellTypeSignature"))
#' Identifying cell-type specific transcription factor regulons for indicated cluster
#'
#' @param x An object with gene expression data used for regulon activity evaluation.
#' Could be \code{matrix}, \code{SingleCellExperiment}.
#' @param pos_clust Indicated cluster used to calculate the cell-type specific score (CSS).
#' @param col_data A \code{data.frame} with cell cluster information. Each row
#' represent a cell and each column represent a cell metadata. The row names of
#' col_data must be identical with column names of x.
#' @param clust_name One column name in colData(x)(for \code{SingleCellExperiment})
#' or col_data (for \code{matrix}) specific the clusters used for cell-type
#' specific score calculation.
#' @param neg_clust Indicated cluster used as ccntrol set to calculate the
#' cell-type specific score (CSS) for. Default: NULL, all others cell clusters
#' except for pos_clust.
#' @param scale Whether to scale values for calculation. If the activate scores
#' have negative values, it must be scaled. Default: FALSE.
#' @param nperm Number of permutations used for calculating p-values against CSS.
#' Default: 1000.
#' @param ncores Number of threads used for parallel calculation. Default: 1.
#' @param verbose Whether to output message? Default: TRUE.
#' @param use_altExp Whether to use altExps in \code{SingleCellExperiment}.
#' Default: TRUE. The activate matrix stored in altExps(x). So when perform
#' cell-type specific analysis on regulon activate matrix, this parameter must be
#' TRUE. When performing expression level analysis, set this parameter to FALSE.
#' @param use_assay Setting the assays used for cell-type specific analysis. For
#' regulon analysis, it could be 'viper', 'aucell' or 'ssgsea'. For gene expression
#' level, it could be 'counts' or 'logcounts'. Default: the first assay.
#'
#' @return A \code{data.frame} with differential CSS for each TFRs.
#' @export
#'
#' @rdname findCellTypeSignature
#' @export
setMethod("findCellTypeSignature","matrix",function(x,
                                                    pos_clust,
                                                    col_data,
                                                    clust_name,
                                                    neg_clust = NULL,
                                                    scale = FALSE,
                                                    nperm = 1000,
                                                    ncores = 1,
                                                    verbose=interactive()){
  mat <- as.matrix(x)
  if(!is.numeric(mat)) stop("The x must be a numeric matrix!")
  if(!identical(sort(colnames(mat)),sort(row.names(col_data)))){
    stop("The colnames of x is not identical with rownames of col_data!")
  }
  if(!(clust_name %in% colnames(col_data))){
    stop("The 'clust_name' is not found in col_data")
  }
  if(!pos_clust %in% col_data[[clust_name]]) stop("The 'pos_clust' must be in ",
                                                  clust_name," column!")
  pos_cells <- row.names(col_data[col_data[[clust_name]] %in% pos_clust,])
  if(!is.null(neg_clust)){
    if(!neg_clust %in% col_data[[clust_name]]) stop("The 'neg_clust' must be in ",
                                                    clust_name," column!")
    neg_cells <- row.names(col_data[col_data[[clust_name]] %in% neg_clust,])
  }else{
    neg_cells <- NULL
  }
  css_res <- calcCSS(mat = mat,
                     pos_cells = pos_cells,
                     use_cluster = pos_clust,
                     neg_cells = neg_cells,
                     scale = scale,
                     nperm = nperm,
                     ncores = ncores,
                     verbose =verbose)
  css_res
})

#' @rdname findCellTypeSignature
#' @export
setMethod("findCellTypeSignature", "SingleCellExperiment", function(x,
                                                                    pos_clust,
                                                                    clust_name,
                                                                    neg_clust = NULL,
                                                                    use_altExp = TRUE,
                                                                    use_assay=1L,
                                                                    scale = FALSE,
                                                                    expression_per = 0.1,
                                                                    nperm = 1000,
                                                                    ncores = 1,
                                                                    verbose=interactive()){
  if(use_altExp){
    mat <- as.matrix(assay(altExp(x=x,e=use_assay)))
  }else{
    mat <- as.matrix(assay(x=x,i=use_assay))
  }
  col_data <- as.data.frame(SummarizedExperiment::colData(x))
  if(!(clust_name %in% colnames(col_data))){
    stop("The 'clust_name' is not found in colData(x), please check the names ",
         "of colData() in your object!")
  }
  if(!pos_clust %in% col_data[[clust_name]]) stop("The 'pos_clust' must be in ",
                                                  clust_name," column!")
  pos_cells <- row.names(col_data[col_data[[clust_name]] %in% pos_clust,])
  if(!is.null(neg_clust)){
    if(!neg_clust %in% col_data[[clust_name]]) stop("The 'neg_clust' must be in ",
                                                    clust_name," column!")
    neg_cells <- row.names(col_data[col_data[[clust_name]] %in% neg_clust,])
  }else{
    neg_cells <- NULL
  }
  css_res <- calcCSS(mat = mat,
                     pos_cells = pos_cells,
                     use_cluster = pos_clust,
                     scale = scale,
                     nperm = nperm,
                     ncores = ncores,
                     verbose =verbose)
  #---pos cluster gene expression
  if(use_altExp){
    if(verbose) message("Summarizing the expression percentge...")
    percent_mat <- summarizeClusterValue(object = x,
                                         clust_name = clust_name,
                                         value = "percent",
                                         gene_name = as.character(css_res$regulon),
                                         use_assay = 1L)
    css_res[['exp_percent']] <- percent_mat[css_res$regulon,pos_clust]
    regulon_list <- sapply(Regulon(x),function(x){paste(x,collapse = ",")})
    css_res[['target']] <- regulon_list[css_res$regulon]
  }
  css_res
})

################################################
# find all cell type specific regulons
setGeneric("findAllCTSRegulons", function(x, ...) standardGeneric("findAllCTSRegulons"))

#' Calculating the significant enriched regulons for each cluster
#'
#' @param x An object with gene expression data used for regulon activity evaluation.
#' Could be \code{matrix}, \code{SingleCellExperiment}.
#' @param col_data A \code{data.frame} with cell cluster information. Each row
#' represent a cell and each column represent a cell metadata. The row names of
#' col_data must be identical with column names of x,
#' @param clust_name One column name in colData(x)(for \code{SingleCellExperiment})
#' or col_data (for \code{matrix}) specific the clusters used for cell-type
#' specific score calculation.
#' @param use_altExp Whether to use altExps in \code{SingleCellExperiment}.
#' Default: TRUE. The activate matrix stored in altExps(x). So when perform
#' cell-type specific analysis on regulon activate matrix, this parameter must be
#' TRUE. When performing expression level analysis, set this parameter to FALSE.
#' @param use_assay Setting the assays used for cell-type specific analysis. For
#' regulon analysis, it could be 'viper', 'aucell' or 'ssgsea'. For gene expression
#' level, it could be 'counts' or 'logcounts'. Default: the first assay.
#' @param fdr The FDR cutoff used for filtering significant regulons. Default:0.05.
#' @param css The CSS cutoff used for filtering significant regulons. Default: 0.3.
#' @param expression_per The expression percentage cutoff for filtering
#' significant regulons. Default:0.1.
#' @param scale Whether to scale values for calculation. If the activate scores
#' have negative values, it must be scaled. Default: FALSE.
#' @param nperm Number of permutations used for calculating p-values against CSS.
#' Default: 1000.
#' @param ncores Number of threads used for parallel calculation. Default: 1.
#' @param verbose Whether to output message? Default: TRUE.
#'
#' @return A \code{data.frame} with differential CSS for each TFRs.
#'
#' @import data.table
#' @export

#' @rdname findAllCTSRegulons
#' @export
setMethod("findAllCTSRegulons","matrix",function(x,
                                                 col_data,
                                                 clust_name,
                                                 fdr = 0.05,
                                                 css = 0.3,
                                                 scale=FALSE,
                                                 nperm = 1000,
                                                 ncores = 1,
                                                 verbose=interactive()){
  #---check data
  mat <- as.matrix(x)
  if(!is.numeric(mat)) stop("The x must be a numeric matrix!")
  if(!identical(sort(colnames(mat)),sort(row.names(col_data)))){
    stop("The colnames of x is not identical with rownames of col_data!")
  }
  if(!clust_name %in% colnames(col_data)) stop(clust_name," not found in col_data!")
  #---get list
  pos_clusters <- unique(as.character(col_data[[clust_name]]))
  pos_cell_list <- lapply(pos_clusters,function(x,col_data,clust_name){
    row.names(col_data[col_data[[clust_name]]==x,])
  },col_data=col_data,clust_name=clust_name)
  names(pos_cell_list) <- pos_clusters
  pos_mat_list <- replicate(n = length(pos_clusters), expr = {
    mat}, simplify = FALSE)
  names(pos_mat_list) <- pos_clusters
  all_res <- purrr::pmap(.l = list(mat=pos_mat_list,
                                   pos_cells=pos_cell_list,
                                   use_cluster=pos_clusters),
                         .f = calcCSS, nperm=nperm, scale=scale,
                         ncores=ncores, verbose = verbose)
  all_res <- data.table::rbindlist(all_res, use.names = TRUE, idcol = clust_name)
  all_res <- dplyr::tibble(all_res)
  all_res <- all_res %>% filter(padj <= fdr & CSS >= css)
  return(all_res)
})

#' @rdname findAllCTSRegulons
#' @export
setMethod("findAllCTSRegulons","SingleCellExperiment",function(x,
                                                               clust_name,
                                                               use_altExp = TRUE,
                                                               use_assay = 1L,
                                                               fdr = 0.05,
                                                               css = 0.3,
                                                               expression_per = 0.1,
                                                               scale=FALSE,
                                                               nperm=1000,
                                                               ncores=1,
                                                               verbose=interactive()){
  if(use_altExp){
    mat <- as.matrix(assay(altExp(x=x,e=use_assay)))
  }else{
    mat <- as.matrix(assay(x=x,i=use_assay))
  }
  regulon_list <- as.list(Regulon(x))
  col_data <- SummarizedExperiment::colData(x)
  if(length(clust_name)!=1) stop("The 'clust_name' must be of length 1!")
  if(!(clust_name %in% colnames(col_data))){
    stop("The 'clust_name' must be one of the colnames of colData(x)!")
  }
  pos_clusters <- unique(as.character(col_data[[clust_name]]))
  pos_cell_list <- lapply(pos_clusters,function(x,col_data,clust_name){
    row.names(col_data[col_data[[clust_name]]==x,])
  },col_data=col_data,clust_name=clust_name)
  names(pos_cell_list) <- pos_clusters
  pos_mat_list <- replicate(n = length(pos_clusters), expr = {
    mat}, simplify = FALSE)
  names(pos_mat_list) <- pos_clusters
  all_res <- purrr::pmap(.l = list(mat=pos_mat_list,
                                   pos_cells=pos_cell_list,
                                   use_cluster=pos_clusters),
                         .f = calcCSS, nperm=nperm, scale=scale,
                         ncores=ncores, verbose = verbose)
  all_res <- data.table::rbindlist(all_res, use.names = TRUE, idcol = clust_name)
  all_res <- dplyr::tibble(all_res)
  if(use_altExp){
    if(verbose) message("Summarizing the expression percentge...")
    percent_mat <- summarizeClusterValue(object = x,
                                         clust_name = clust_name,
                                         value = "percent",
                                         gene_name = unique(as.character(all_res$regulon)),
                                         use_assay = 1L)
    all_res$exp_percent <- apply(all_res,1,function(x, percent_mat){
      percent_mat[x['regulon'],x[clust_name]]},percent_mat=percent_mat)
    all_res$targets <- regulon_list[all_res$regulon]
    all_res <- all_res %>% filter(exp_percent >= expression_per)
  }
  diff_res <- all_res %>% filter(padj <= fdr & CSS>=css)
  diff_res
})

################################################
#---calculate cell-type specific score
################################################
#' Identifying cell-type specific transcription factor regulons
#'
#' @param mat A numeric \code{matrix} contains transcription factor regulons activate score.
#' @param pos_cells Cells from the cluster(s) to calculate the cell-type specific
#' score. It must be included in column names of mat.
#' @param neg_cells Cells from the control cluster(s) to calculate the cell-type
#' specific score. It must be included in column names of mat. Default: NULL,
#' all cells expect for pos_cells.
#' @param use_cluster The cluster used for cell-type specific calculation.
#' Default: NULL.
#' @param scale Whether to scale values for calculation. If the activate scores
#' have negative values, it must be scaled. Default: FALSE.
#' @param nperm Number of permutations used for calculating p-values against CSS.
#' Default:1000.
#' @param ncores Number of threads used for parallel calculation. Default: 1.
#' @param verbose Whether to output message? Default: TRUE.
#'
#' @importFrom purrr map2 map
#' @importFrom parallel makeCluster stopCluster
#' @importFrom pbapply pblapply
#'
#' @return A \code{data.frame} with differential CSS for each TFRs.
#'
calcCSS <- function(mat, pos_cells, use_cluster=NULL, neg_cells=NULL, scale = FALSE,
                    nperm=1000, ncores=1, verbose = interactive()){
  if(!is.null(use_cluster) && verbose){
    message(use_cluster,": ",nrow(mat)," regulator(s).")
  }
  mat <- as.matrix(mat)
  if(is.null(neg_cells)){
    neg_cells <- setdiff(colnames(mat),pos_cells)
  }else{
    neg_cells <- intersect(colnames(mat),neg_cells)
  }
  pos_cells <- intersect(colnames(mat),pos_cells)
  if(length(intersect(pos_cells, neg_cells))>0){
    warning("The pos_cells and neg_cells have overlaps!")
  }
  mat <- mat[, c(pos_cells, neg_cells), drop = FALSE]
  if(nrow(mat)==0 || ncol(mat)==0) stop("No activate score found in matrix!")
  #---normalization
  if(scale || min(mat)<0) mat <- t(apply(mat,1,function(x){(x-min(x))/(max(x)-min(x))}))
  #-----calculate average activate score
  pos_score <- apply(mat[,pos_cells,drop=FALSE],1,mean)
  neg_score <- apply(mat[,neg_cells,drop=FALSE],1,mean)
  ave_diff <- pos_score - neg_score
  #---get CSS list
  mat_list <- lapply(X = row.names(mat),FUN = function(x, mat){
    mat[x,]/sum(mat[x,])},mat=mat)
  names(mat_list) <- row.names(mat)
  const_freq <- rep(0,ncol(mat))
  names(const_freq) <- colnames(mat)
  const_freq[pos_cells] <- 1/length(pos_cells)
  #--par
  cl <- mkCluster(ncores = ncores)
  parallel::clusterExport(cl = cl, varlist = c(".getCSS",".calcJSsp",".H"), envir = environment())
  CSS <- pbapply::pblapply(X = mat_list, FUN = .getCSS, const_freq=const_freq, cl = cl)
  parallel::stopCluster(cl = cl)
  CSS <- as.data.frame(do.call(rbind,CSS))
  CSS$padj <- p.adjust(CSS$pval, method = "BH")
  CSS <- data.frame(regulon = row.names(CSS),
                    pos_score = pos_score[row.names(CSS)],
                    neg_score = neg_score[row.names(CSS)],
                    ave_diff = ave_diff[row.names(CSS)],
                    CSS, check.rows = TRUE)
  CSS[order(CSS$padj,CSS$CSS,decreasing = c(FALSE, TRUE),method = "radix"),]
}
#---calculating css score
.getCSS <- function(pos_freq, const_freq, nperm=1000){
  regulon_css <- .calcJSsp(pos_freq, const_freq)
  #---permutation
  nperm_freq <- replicate(n = nperm, expr = {
    sample(x = unname(pos_freq), size = length(pos_freq),replace = FALSE)},
    simplify = FALSE)
  nperm_css <- purrr::map2(.x = nperm_freq,.y = list(const_freq),.f = .calcJSsp)
  #------
  pval <- length(nperm_css[which(nperm_css>regulon_css)])/nperm
  c(CSS = regulon_css, pval = pval)
}

#---calculate cell-type specific score using JS divergence
.calcJSsp <- function(p1, p2, unit = c("log","log2","log10")){
  unit <- match.arg(unit)
  #---calculate JSD
  JSD <-.H(0.5*(p1+p2), unit = unit) - 0.5*(.H(p1, unit = unit)+.H(p2, unit = unit))
  1-sqrt(JSD)
}

#---calculate entropy from probability
.H <- function(freqs, unit = c("log", "log2", "log10")){
  unit = match.arg(unit)
  #freqs = freqs/sum(freqs)
  H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
  if (unit == "log2") H = H/log(2)
  if (unit == "log10") H = H/log(10)
  H
}

#---calculating the average value for genes in each cluster

#' Summarizing the values for genes in each cluster
#'
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param clust_name One column name in colData(x)(for \code{SingleCellExperiment})
#' clusters used for summarizing the value.
#' @param gene_name Gene names used for summarizing. Default: NULL. All genes in
#' object.
#' @param value The values calculated for each cell cluster. Default: mean.
#' @param scale Whether to scale the value?
#' @param use_assay Setting the assays with the expression data. It could be
#' 'counts', 'logcounts' or 'TPM'. Default: the first assay.
#'
#' @return A matrix with average values for each cluster.
#' @export
#'
#' @examples
summarizeClusterValue <- function(object,
                                  clust_name,
                                  value = c("mean","median","percent"),
                                  gene_name = NULL,
                                  scale=FALSE,
                                  use_assay = 1L){
  value <- match.arg(value)
  if(!is(object,"SingleCellExperiment")){
    stop("The object must be a SingleCellExperiment object!")
  }
  col_data <- SummarizedExperiment::colData(object)
  if(length(clust_name)!=1) stop("The 'clust_name' must be of length 1!")
  if(!(clust_name %in% colnames(col_data))){
    stop("The 'clust_name' must be one of the colnames of colData(object)!")
  }
  mat <- as.matrix(assay(x=object,i=use_assay))
  if(!is.null(gene_name)){
    gene_name_overlap <- intersect(gene_name,row.names(object))
    if(length(gene_name_overlap)==0) stop("All your input genes not found in your object!")
    gene_name_left <- setdiff(gene_name, row.names(object))
    if(length(gene_name_left)> 0){
      warning("The following genes were not found in your object: \n",gene_name_left)
    }
    gene_name <- gene_name_overlap
  }else{
    gene_name <- row.names(object)
  }
  mat <- mat[gene_name,,drop=FALSE]
  if(scale){
    mat <- t(apply(mat,1,function(x){(x-min(x))/(max(x)-min(x))}))
  }
  if(value=="mean"){
    mean_mat <- t(pbapply::pbapply(mat,1,function(x,col_data,clust_name){
      tapply(X = x, INDEX = as.character(col_data[names(x),clust_name]),FUN = mean)
    },col_data=col_data,clust_name=clust_name))
  }else if(value=="median"){
    mean_mat <- t(pbapply::pbapply(mat,1,function(x,col_data,clust_name){
      tapply(X = x, INDEX = as.character(col_data[names(x),clust_name]),FUN = median)
    },col_data=col_data,clust_name=clust_name))
  }else if(value=="percent"){
    mean_mat <- t(pbapply::pbapply(mat,1,function(x,col_data,clust_name){
      tapply(X = x, INDEX = as.character(col_data[names(x),clust_name]),
             FUN = function(y){
               length(y[y>0])/length(y)
             })
    },col_data=col_data,clust_name=clust_name))
  }
  mean_mat
}

#' Summarizing the values for regulons in each cluster
#'
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param clust_name One column name in colData(x)(for \code{SingleCellExperiment})
#' specific the clusters used for mean value calculation.
#' @param regulon_name Gene names used for mean value calculation. Default: NULL.
#' All genes in object.
#' @param value The values calculated for each cell cluster. Default: mean.
#' @param scale Whether to scale the value?
#' @param use_assay Setting the alternative assays with the regulon activate
#' score matrix. It could be 'viper', 'aucell' or 'ssgsea'. Default: the first assay.
#'
#' @return A matrix with average values for each cluster.
#' @export
#'
#' @examples
summarizeClusterValueRegulon <- function(object,
                                         clust_name,
                                         value = c("mean","median","percent"),
                                         regulon_name = NULL,
                                         scale=FALSE,
                                         use_assay = 1L){
  value <- match.arg(value)
  atfr <- altExp(x = object, e = use_assay)
  mean_act_mat <- summarizeClusterValue(object = atfr,
                                        clust_name = clust_name,
                                        gene_name = regulon_name,
                                        value = value,
                                        scale=scale,
                                        use_assay = "atfr")
  mean_act_mat
}







