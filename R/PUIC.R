#' Inferring Gene Regulatory Network (GRN) using partial information decomposition and context (PUIC) method
#'
#' The PUIC package implements the method in paper
#' (\href{https://linkinghub.elsevier.com/retrieve/pii/S2405471217303861}{2017 Chan et al.})
#' with few modifications. Besides, The code used for calculating the \code{Multivariate
#' Information} (MI) and \code{specific.information}, discretizing gene expression data
#' were  partly adopted from package \code{\link{Informeasure}}.
#'
#' @param expMat A \code{matrix} storing the gene expression data. Rows corresponding
#' to features (eg. gene symbols) and columns corresponding to samples (cells).
#' Raw read counts, UMIs, TPMs or logNormalized counts were supported.
#' @param regulators The regulator genes used for GRN inferring (eg, transcription
#' factors). At least two regulators required. Default: NULL, using all the genes
#' in gene expression matrix. Default: NULL, all genes.
#' @param targets The target genes used for GRN inferring. Default: NULL, using
#' all the genes in gene expression matrix. Default: NULL, all genes.
#' @param logScale Whether to log-scale the input data? Default: FALSE.
#' @param ncores Number of cores used for parallel calculation. The running time
#' heavily depend on the number of regulators and targets, for genes over 5000,
#' we strongly suggest you to use multi-cores. Default: 1.
#' @param diag Numeric. The weight in the diagonal of output square matrix. Default: 1.
#' @param verbose Whether to print message while running. Default: TRUE.
#'
#' @return A matrix with weighted value between regulators and targets.
#'
#' @importFrom purrr map map2 pmap transpose
#' @importFrom parallel makeCluster clusterExport clusterEvalQ stopCluster
#' @importFrom pbapply pbapply pblapply
#' @importFrom methods as is
#' @importFrom stats ecdf
#' @export
#'
#' @examples
#' data(expMat)
#' PUIC_res <- PUIC(expMat)
#' head(PUIC_res)
#'
#'
#'
PUIC <- function(expMat, regulators=NULL, targets=NULL, logScale=FALSE,
                 ncores=1, diag=c("auto","zero","one"), verbose = interactive()){
  if(verbose) message("[1] Filtering and Normalizing data...")
  .checkPUICArgs(expMat=expMat, regulators=regulators, targets=targets,
                 logScale=logScale, ncores=ncores)
  diag <- match.arg(diag)
  #---check expMat
  expMat <- as.matrix(expMat)
  if(logScale){
    expMat <- log2(expMat+1)
  }
  #---check regulators
  if(is.null(regulators)){
    regulators <- row.names(expMat)
  }
  #---check targets
  if(is.null(targets)){
    targets <- row.names(expMat)
  }
  ######calculating
  #1.discretize Gene
  if(verbose) message("[2] Discretizing expression matrix into bins...")
  expMat <- as.list(as.data.frame(t(expMat),check.names=FALSE))
  discret_list <- pbapply::pblapply(X=expMat, FUN = .discretizeGene)
  regulators <- intersect(regulators, names(discret_list))
  targets <- intersect(targets, names(discret_list))
  #2. Calculating the proportional unique contribution (PUC) score matrix
  if(verbose) message("[3] Calculating proportional unique contribution(PUC) ",
                      "matrix using ",length(regulators)," regulators and ",
                      length(targets)," targets...\n(This step may taken dozens ",
                      "of hours, please be patient)")
  #---parallel calculation
  if(.Platform$OS.type=="unix"){
    cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores), type = "FORK")
  }else{
    cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores),type = "PSOCK")
  }
  parallel::clusterEvalQ(cl,library(purrr))
  parallel::clusterEvalQ(cl,library(dplyr))
  parallel::clusterEvalQ(cl,library(pbapply))
  PUC_list <- pbapply::pblapply(X = targets, FUN = .getPUC,
                                regulators = regulators,
                                discret_list = discret_list,
                                cl=cl)
  parallel::stopCluster(cl)
  names(PUC_list) <- targets
  PUC_mat <- do.call(cbind, PUC_list)
  #3. Summaring Uxy and calculating the weight matrix
  if(verbose) message("[4] Summaring Uxy and calculating the weight matrix...")
  out_res <- .FUxy(Uxy_mat = PUC_mat)
  if(diag=="auto"){
    return(out_res)
  }else if(diag=="zero"){
    diag(out_res) <- 0
  }else if(diag=="one"){
    diag(out_res) <- 1
  }
  out_res
}


#' Inferring gene regulons from Gene Regulatory Network
#'
#' Inferring gene regulons from gene regulatory network using weighted square matrix from PUIC output.
#'
#' @param weightMat A square matrix whose element should be positive values.
#' @param methods The name of the network inference algorithm. Default: aracne.
#' @param cutoff Set the cutoff of regulatory networks. Default: NULL, means
#' inferring cutoff by itself.
#'
#' @return A data.frame or a list with regulatory networks.
#' @importFrom minet aracne clr mrnet mrnetb
#' @importFrom reshape2 melt
#' @importFrom stats sd ks.test rnorm qnorm density dnorm
#' @importFrom utils installed.packages capture.output
#' @importFrom parallel makeCluster stopCluster
#' @importFrom pbapply pblapply
#'
#' @export
#'
#' @examples
#' data(expMat)
#' PUIC_res <- PUIC(expMat)
#' PUIC_net <- matToNet(PUIC_res)
matToNet <- function(weightMat,
                     methods = c("aracne", "clr", "mrnet", "mrnetb"),
                     cutoff = NULL){
  methods <- match.arg(methods)
  mat <- as.matrix(weightMat)
  if(nrow(mat)!=ncol(mat)) stop("The input matrix must be square matrix!")
  if(min(mat)<0) stop("The input matrix can not have negative values!")
  if(!is.null(cutoff)){
    if(!is.numeric(cutoff) || !(cutoff>=0 && cutoff<=1) || length(cutoff)!=1){
      stop("The cutoff must be a numeric value between 0 and 1!")
    }
  }
  message("[1] Generating Gene regulatory networks using: ", methods, ".")
  if(methods=="aracne"){
    net <- minet::aracne(mat)
  }else if(methods=="clr"){
    net <- minet::clr(mat)
  }else if(methods=="mrnet"){
    net <- minet::mrnet(mat)
  }else if(methods=="mrnetb"){
    net <- minet::mrnetb(mat)
  }
  out_tab <- reshape2::melt(data = net)
  colnames(out_tab) <- c("regulator","target","weight")
  if(is.null(cutoff)){
    out_res <- out_tab
  }else{
    out_res <- out_tab[out_tab$weight>cutoff,]
  }
  row.names(out_res) <- NULL
  out_res[order(out_res$regulator,out_res$weight,out_res$target,
                decreasing = c(FALSE,TRUE,FALSE), method = "radix"),]
}

######
#---1. check input arguments
######
.checkPUICArgs <- function(expMat, regulators, targets, logScale, ncores){
  #---check expMat
  if(!is.matrix(expMat) && !is.array(expMat) && !is(expMat,"dgCMatrix")){
    stop("The expMat must be a two-dimensional matrix where the row corresponds",
         " to a gene and each column corresponds to a sample")
  }
  if (length(dim(expMat)) != 2) {
    stop("The expMat must be a two-dimensional matrix where the row corresponds",
         " to a gene and each column corresponds to a sample")
  }
  if (is.null(rownames(expMat))) {
    stop("The expMat must contain the names of the genes as rownames.")
  }
  expMat <- as.matrix(expMat)
  if(!is.numeric(expMat)) stop("The expMat contain non-numeric values.")
  if(any(is.na(expMat))) stop("The expMat contain NA.")
  gene_names <- unique(rownames(expMat))
  #---check regulators in expMat
  regulators <- unique(regulators)
  if(!is.null(regulators)){
    if(!is.character(regulators)) stop("The regulators must characters!")
    if(length(intersect(gene_names, regulators))<2){
      stop("At least two regulators requiered, please check the names of the ",
           "expMat and the gene ids (names) in your regulators!")
    }else if(length(intersect(gene_names, regulators)) < length(regulators)){
      ratio <- length(setdiff(regulators,gene_names))/length(regulators)
      warning(length(setdiff(regulators,gene_names))," out of ",length(regulators),
              " (",round(ratio*100,2),"%) regulators not found in your expMat!")
    }
  }
  #---check targets in expMat
  targets <- unique(targets)
  if(!is.null(targets)){
    if(!is.character(targets)) stop("The targets must characters!")
    if(length(intersect(gene_names, targets))<1){
      stop("At least one targets requiered, please check the names of the ",
           "expMat and the gene ids (names) in your targets!")
    }else if(length(intersect(gene_names, targets)) < length(targets)){
      ratio <- length(setdiff(targets,gene_names))/length(targets)
      warning(length(setdiff(targets,gene_names))," out of ",length(targets),
              " (",round(ratio*100,2),"%) targets not found in your expMat!")
    }
  }
  if(!is.logical(logScale)) stop("The logScale must be TRUE(T) or FALSE(F)!")
  if(!is.numeric(ncores) || ncores<1) stop("The ncores must be positive numbers!")
}
######
#---2. discrete one gene expression into equal length items
######
.discretizeGene <- function(x, method = c("uniform_width","uniform_frequency")){
  method <- match.arg(method)
  if(method=="uniform_width"){
    numBins <- floor(sqrt(length(x)))
    r <- range(x)
    b <- seq(from = r[1], to = r[2], length.out = (numBins + 1))
    X <- cut(x, breaks = b , include.lowest = TRUE)
  }else if("uniform_frequency"){
    numBins <- floor(sqrt(length(x)))
    nrepl <- floor(length(x)/numBins)
    nplus <- sample(seq_len(numBins), length(x) - nrepl * numBins)
    nrep <- rep(nrepl, numBins)
    nrep[nplus] <- nrepl + 1
    X <- x[order(x)] <- rep(seq.int(numBins), nrep)
  }
  return(X)
}
######
#---3. get pairwise probability from  frequency table
######
# x discrete factor
# y discrete factor
.freqTable <- function(x, y, method=c("ML", "Jeffreys", "Laplace", "SG", "minimax")){
  method <- match.arg(method)
  if(is.list(x) || is.list(y)){
    x <- unlist(x)
    y <- unlist(y)
  }
  tab <- table(x,y)
  dimnames(tab) <- NULL
  if(method=="ML"){
    probs <- tab/sum(tab)
  }else if(method=="Jeffreys"){
    probs <- (tab + 0.5)/sum(tab + 0.5)
  }else if(method=="Laplace"){
    probs <- (tab + 1)/sum(tab + 1)
  }else if(method=="SG"){
    probs <- (tab + length(tab))/sum(tab + length(tab))
  }else if(method=="minimax"){
    probs <- (tab + sqrt(sum(tab))/length(tab))/sum(tab + sqrt(sum(tab))/length(tab))
  }
  probs
}

######
#---3.get all regulator and one target pair PUC scores
######
.getPUC <- function(target, discret_list, regulators = names(discret_list)){
  if(length(regulators)<2) stop("At least two regulators needed!")
  if(length(target)!=1) stop("Only support one target each time!")
  if(is.null(names(discret_list))) stop("The names of discret_list must be gene symbols!")
  Z <- discret_list[target]
  Xi <- discret_list[regulators]
  p_XiZ <- purrr::map2(.x = Xi, .y = Z, .f = .freqTable)
  p_Xi <- purrr::map(.x = p_XiZ, .f = rowSums)
  p_Z <- purrr::map(.x = p_XiZ, .f = colSums)
  I_XiZ <- purrr::pmap(.l=list(freqs=p_XiZ, p_i=p_Xi, p_z=p_Z), .f = .MI)
  Ispec_XiZ <- purrr::pmap(.l=list(p_iz=p_XiZ, p_i=p_Xi, p_z=p_Z), .f = .specific.information)
  U_xy <- .puc_per_target(I_XiZ = I_XiZ, Ispec_XiZ = Ispec_XiZ, p_Z = p_Z)
  return(U_xy)
}

#' I_XiZ, Ispec_XiZ and p_Z is list with equally length, theoretically, p_Z is the same value.
.puc_per_target <- function(I_XiZ, Ispec_XiZ, p_Z){
  gene_names <- names(Ispec_XiZ)
  Ispec_XiZ_bins <- lapply(purrr::transpose(Ispec_XiZ),unlist)
  Ispec_XiZ_bins_pos <- which(sapply(Ispec_XiZ_bins,sum)>0)
  p_Z_bins <- lapply(purrr::transpose(p_Z),unlist)
  p_Z_bins_pos <- which(sapply(p_Z_bins,sum)>0)
  if(!identical(Ispec_XiZ_bins_pos, p_Z_bins_pos)){
    warning("The bins with values may not match between Ispec and p_Z!")
  }
  overlap_bins_pos <- union(Ispec_XiZ_bins_pos, p_Z_bins_pos)
  Ispec_XiZ_bins_filter <- Ispec_XiZ_bins[overlap_bins_pos]
  p_Z_bins_filter <- p_Z_bins[overlap_bins_pos]
  #---calculate sum redundance of Xi Z for all Yi
  per_bin_redundancy <- purrr::map2(.x = Ispec_XiZ_bins_filter,
                                    .y = p_Z_bins_filter,
                                    .f = function(x, y){
                                      if(is.null(names(x)) || is.null(names(y))){
                                        stop("The Ispec and p_Z in each bin must have the same gene names!")
                                      }
                                      x_sort <- sort(x)
                                      n_gene <- length(x_sort)
                                      x_cumsum <- cumsum(x_sort)
                                      x_over_num <- n_gene-order(x_sort)-1
                                      names(x_over_num) <- names(x_sort)
                                      x_over_sum <- x_sort*x_over_num
                                      x_value <- x_cumsum+x_over_sum
                                      y_value <- y[names(x_value)]
                                      out <- x_value*y_value
                                      out[names(x)]
                                    })
  #---check gene name
  check_gnames <- sapply(per_bin_redundancy,function(x,gname){
    identical(names(x),gname)},gname=gene_names)
  if(!all(check_gnames)){
    per_bin_redundancy <- lapply(per_bin_redundancy,function(x,gnames){
      x[gnames]},gname=gene_names)
  }
  per_gene_redundancy <- sapply(purrr::transpose(per_bin_redundancy),
                                function(x){sum(unlist(x))})
  per_gene_puc <- length(gene_names)-1-per_gene_redundancy/unlist(I_XiZ)[names(per_gene_redundancy)]
  return(per_gene_puc)
}

#---calculation of specific.information for each regulator-target pair
.specific.information <- function(p_iz, p_i, p_z, unit = c("log", "log2", "log10")){
  unit <- match.arg(unit)
  p_i_z <- t(t(p_iz)/p_z)  ##p(i|z)
  p_z_i <- t(p_iz/p_i) ##p(z|i)
  tmp <- t((log(1/p_z)) - (log(1/p_z_i))) * p_i_z
  if (unit == "log2")  tmp <- t((log(1/p_z, 2))  - (log(1/p_z_i, 2))) * p_i_z  # change from log to log2 scale
  if (unit == "log10") tmp <- t((log(1/p_z, 10)) - (log(1/p_z_i, 10))) * p_i_z # change from log to log10 scale
  colSums(tmp, na.rm = TRUE)
}

#---calculation of mutual information
.MI <- function(freqs, p_i, p_z, unit = c("log", "log2", "log10")){
  if(length(dim(freqs))!=2) stop("The dims of freqs must be 2!")
  unit <- match.arg(unit)
  MI <- .H(p_i, unit = unit) + .H(p_z, unit = unit) - .H(freqs, unit = unit)
  MI
}

.H <- function(freqs, unit = c("log", "log2", "log10")){
  unit = match.arg(unit)
  #freqs = freqs/sum(freqs)
  H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
  if (unit == "log2") H = H/log(2)
  if (unit == "log10") H = H/log(10)
  H
}

######
#---4. get the cumulative distribution score of PUC
######
.FUxy <- function(Uxy_mat){
  mat <- as.matrix(Uxy_mat)
  if(nrow(mat)==ncol(mat)){
    mat <- mat + t(mat)
  }
  F_x <- t(apply(mat,1,function(x){stats::ecdf(x)(x)}))
  colnames(F_x) <- colnames(mat)
  F_y <- apply(mat,2,function(x){stats::ecdf(x)(x)})
  row.names(F_y) <- row.names(mat)
  F_xy <- (F_x+F_y)*0.5
  return(F_xy)
}

