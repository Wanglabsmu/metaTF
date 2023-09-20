#' Integrating weighted regulator-target regulatory networks
#'
#' This method integrate weighted regulatory networks from different resources
#' and filter the top ranked regulatory-target pairs for further analysis.
#'
#' @param grnNets List. A list of weighted regulator-targets networks from
#' different resources, for example: motif prediction, Chip-seq and Knockout
#' experiment. For every network in this list, three columns corresponding to
#' regulator, target and weight must be offered. If one of your network do not
#' have weights, you can set 1 to all the regulator-target pairs in 3-rd column.
#' @param weight Numeric vector. Positive numbers of the same length with
#' \code{grnNets} specific the weight from each resource. default: 1.
#' @param norm_weight Logical. Whether to normalize the weight to [0,1] of
#' regulator-target pairs from each resource. Default: FALSE.
#' @param keep_neg_regu Logical. Whether to keep negative regulatory
#' relationships between regulator and targets. Default: FALSE.
#' @param join_type Character. The merge types among different \code{data.frame},
#' corresponding to \code{\link{full_join}}, \code{\link{inner_join}},
#' \code{\link{left_join}} and \code{\link{right_join}} methods in \pkg{dplyr}.
#' Default: full.
#' @param return_sig Logical. Only return top weighted regulatory networks?
#' Default: FALSE.
#' @param max_target_num Numeric. Set the maximum number for secondly filtering.
#' Number over this value will filter for multiple times. Default: 5000.
#' @param min_target_num Numeric. Set the minimum cutoff numbers of targets from
#' each regulator for filtering. Number under this value will not filter.
#' Default: 100.
#' @param regulon_list Logical. Return a \code{list} of regulons or a
#' \code{data.frame} with networks. Default: \code{data.frame}.
#' @param beta0 Numeric. The sparsity parameter.
#' @param maxInter Maximum number of iteration.
#' @param ncores Number of cores use. Default: 1.
#' @param verbose Whether to show information?
#'
#' @return A \code{data.frame} with integrated weighted regulator-target
#' regulatory networks.
#'
#' @importFrom dplyr full_join inner_join left_join right_join
#' @importFrom purrr reduce pmap_dbl
#'
#' @export
#'
#' @examples
integraNets <- function(grnNets,
                        weight = rep(1,length(grnNets)),
                        norm_weight = FALSE,
                        keep_neg_regu = FALSE,
                        join_type = c("full","inner","left","right"),
                        return_sig = FALSE,
                        max_target_num = 5000L,
                        min_target_num = 100L,
                        regulon_list=FALSE,
                        beta0 = 0,
                        maxInter=1,
                        ncores=1,
                        verbose = interactive()){
  #------check arguments
  if(!is.list(grnNets)) stop("The parameter grnNets must be a list!")
  if(length(grnNets)<2) stop("There need more than two (included) networks for integration!")
  if(verbose) message("[1] Checking the input list...")
  grnNets <- .check_df(grnNets)
  net_names <- names(grnNets)
  if(is.null(net_names)){
    net_names <- paste("Net",1:length(grnNets),sep = "")
  }
  names(grnNets) <- net_names
  #---merge all tables
  if(verbose) message("[2] Mergeing the weighted data.frame ...")
  join_type <- match.arg(join_type)
  if(join_type=="full"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::full_join,by=c("tf","target"))
  }else if(join_type=="inner"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::inner_join,by=c("tf","target"))
  }else if(join_type=="left"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::left_join,by=c("tf","target"))
  }else if(join_type=="right"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::right_join,by=c("tf","target"))
  }
  grn_tbl[is.na(grn_tbl)] <- 0
  #---weight direction
  if(keep_neg_regu){
    if(verbose) message("Getting regulatory direction ...")
    regu_direc <- apply(grn_tbl[,-c(1:2),drop=FALSE],1,function(x){
      prod(x)/abs(prod(x))
    })
  }
  #--- check weights and normalization
  if(!is.numeric(weight)) stop("The weight must be numeric!")
  if(length(weight)!=length(grnNets)) stop("The weight must be the same length with grnNets!")
  names(weight) <- net_names
  if(verbose) message("[3] Normalizing the merged weight...")
  grnNets <- pbapply::pblapply(net_names, function(x, grnNets, weight, norm_weight){
    tbl_one <- grnNets[[x]]
    colnames(tbl_one) <- c("regulator","target", x)
    if(norm_weight){
      tbl_one[[x]] <- (tbl_one[[x]]-min(tbl_one[[x]]))/(max(tbl_one[[x]])-min(tbl_one[[x]]))
    }
    tbl_one[[x]] <- tbl_one[[x]]*weight[x]
    tbl_one
  }, grnNets=grnNets, weight=weight, norm_weight=norm_weight)
  names(grnNets) <- net_names
  if(verbose) message("[4] Calculating the merged weight...")
  merge_weight <- purrr::pmap_dbl(.l = grn_tbl[,-c(1:2)],
                                  .f = function(..., beta0=0){
                                    x <- as.numeric(unlist(list(...)))
                                    1/(1+exp(-1*sum(beta0,abs(x))))
                                    }, beta0=beta0)
  final_tbl <- dplyr::bind_cols(grn_tbl[,1:2],weight=merge_weight) %>% unique(.)
  final_tbl$weight <- (final_tbl$weight-min(final_tbl$weight))/(max(final_tbl$weight)-min(final_tbl$weight))
  merged_net <- final_tbl[order(final_tbl$tf,final_tbl$weight,final_tbl$target,
                                decreasing=c(FALSE, TRUE, FALSE), method = "radix"),]
  row.names(merged_net) <- NULL
  #---get the cutoff
  if(return_sig){
    if(verbose) message("[5] Inferring weight cutoff and filtering each networks...")
    merged_net <- filterNets(df = merged_net,
                             minSize = min_target_num,
                             maxSize = max_target_num,
                             maxInter = maxInter,
                             ncores = ncores)

  }
  #---get regulon
  if(regulon_list){
    merged_net <- split(merged_net$target,merged_net$tf)
  }
  merged_net
}

#---check the GRNs in your data
.check_df <- function(df_list){
  df_class <- sapply(df_list, is.data.frame)
  if(!all(df_class)) stop("Not all the data stored in you list is data.frame!")
  df_cols <- sapply(df_list, ncol)
  if(any(df_cols < 2)) stop("There must be at least 2 columns in your data.frame!")
  if(any(df_cols > 3)){
    warning("There are more than 3 columns in your data.frame, only using the first three columns!")
  }
  df_list <- lapply(df_list,function(x){
    if(ncol(x)>3){
      x <- x[,1:3]
      if(!is.numeric(x[[3]])) warning("The third column is not numeric")
      x[[3]] <- as.numeric(x[[3]])
    }else if(ncol(x)==2){
      x[[3]] <- rep(1,nrow(x))
    }
    colnames(x) <- c("tf","target","weight")
    x
  })
  df_list
}


#---calculate the cutoff of score using method from AUCell package with some modification
#' Filtering the weighted regulatory network
#'
#' Filtering the regulatory network by evaluating the threshold of each regulator
#' based on the distribution of weight.
#'
#' @param df A data.frame with only one regulator
#' @param minSize Minimum number of targets to keep. Default: 100.
#' @param maxSize Maximum number of targets to filter. Default: 5000.
#' @param maxInter Maximum number of iteration. At least 1. Default: 5.
#' @param ncores Number of cores to use. Default: 1.
#'
#' @importFrom pbapply pblapply
#' @importFrom parallel makeCluster clusterExport stopCluster
#'
#' @return A weighted \code{data.frame} with regulator and target information.
#' @export
#'
#' @examples
#' net <- readRDS(system.file("./inst/extdata/mm_pantissue_network.rds",package = "metaTF"))
#' net <- head(net,100000)
#' net <- filterNets(df=net)
#'
filterNets <- function(df, minSize=100L, maxSize=5000L, maxInter=1L, ncores=1L){
  if(!is(df,"data.frame")) stop("The df must be data.frame!")
  if(ncol(df)<3) stop("The df must at least 3 columns!")
  if(!is.numeric(maxInter) || maxInter<0) stop("The 'maxInter' must be positive numbers.!")
  df <- as.data.frame(df)
  df_list <- split(df, as.character(df[[1]]))
  if(.Platform$OS.type=="unix"){
    cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores), type = "FORK")
  }else{
    cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores), type = "PSOCK")
  }
  parallel::clusterEvalQ(cl = cl, library(metaTF))
  parallel::clusterExport(cl = cl, varlist = c(".getWeightThreshold", ".getAUCThreshold"), envir = environment())
  filter_res <- pbapply::pblapply(X = df_list, FUN = .getWeightThreshold,
                                  minSize=minSize, maxSize=maxSize,
                                  maxInter=maxInter, cl = cl)
  parallel::stopCluster(cl)
  df <- do.call(rbind, filter_res)
  row.names(df) <- NULL
  df
}

.getWeightThreshold <- function(df, minSize=100L, maxSize=5000L, maxInter=5L, pb = NULL){
  #---filter regulons
  if(!is.null(pb)) pb$tick()
  df <- unique(df)
  if(nrow(df) <= minSize){
    df_final <- df
  }else{
    inters <- 0
    while(nrow(df) > minSize){
      auc <- df[["weight"]]
      names(auc) <- df[['target']]
      auc_thres <- .getAUCThreshold(auc = auc)
      if(length(auc[auc>=auc_thres]) < minSize){
        break;
      }else{
        df <- df[df$weight>=auc_thres,]
      }
      inters <- inters+1
      if(inters > maxInter) break;
    }
    df_final <- df
  }
  df_final <- df_final[order(df_final$tf,df_final$weight,df_final$target,
                             decreasing = c(FALSE,TRUE,FALSE),method = "radix"),]
  row.names(df_final) <- NULL
  if(nrow(df_final) > maxSize){
    df_final <- head(df_final,maxSize)
  }
  return(df_final)
}

#----calculate AUC
.getAUCThreshold <- function(auc, smallestPopPercent=.1,
                             densAdjust=2, thrP=0.01, nBreaks=100){
  gSetName <- names(auc)
  nCells <- length(auc)
  skipGlobal <- TRUE
  skipRed <- FALSE
  skipSmallDens <- FALSE
  commentMsg <- ""
  aucThrs <- c()

  notPopPercent <- 1 - smallestPopPercent
  if(sum(auc==0) > (nCells*notPopPercent))
  {
    skipGlobal <- FALSE
    commentMsg <- paste(commentMsg,
                        round((sum(auc==0)/nCells)*100),
                        "% (more than ", notPopPercent,"%) of AUC are zero. ", sep="")
  }

  meanAUC <- mean(auc)
  sdAUC <- sd(auc)
  maybeNormalDistr <- !suppressWarnings(
    ks.test(auc, rnorm(max(100,length(auc)),mean=meanAUC, sd = sdAUC),
            alternative = "less")$p.value < .01)
  if(maybeNormalDistr){
    commentMsg <- paste0(commentMsg,
                         "The AUC might follow a normal distribution (random gene-set?). ")
    skipGlobal <- FALSE

    # aucThrs["outlierOfGlobal"] <- meanAUC + 2*sdAUC
    aucThrs["outlierOfGlobal"] <- qnorm(1-(thrP/nCells), mean=meanAUC, sd=sdAUC)
  }

  #V6
  histogram <- hist(c(0, auc/max(auc)), breaks=100, plot=FALSE)$count
  if((sum(histogram[1:5]) / sum(histogram)) >= notPopPercent*.75) {
    skipGlobal <- FALSE
    skipRed <- TRUE
    skipSmallDens <- TRUE
  }
  if((sum(histogram[1:10]) / sum(histogram)) >= notPopPercent*.50) {
    skipSmallDens <- TRUE
    skipGlobal <- FALSE
    # skipRed <- TRUE ?
    aucThrs["tenPercentOfMax"] <- max(auc)*.10
  }
  # print(skipRed)

  densCurve <- density(auc, adjust=densAdjust, cut=0)
  maximumsDens <- NULL
  inflPoints <- diff(sign(diff(densCurve$y)))
  maximumsDens <- which(inflPoints==-2)
  globalMax <- maximumsDens[which.max(densCurve$y[maximumsDens])]
  minimumDens <- which(inflPoints==2)
  smallMin <- NULL
  if(!skipSmallDens)
    smallMin <- data.table::last(minimumDens[which(minimumDens < globalMax)]) #1prev to max
  minimumDens <- c(smallMin,
                   minimumDens[which(minimumDens > globalMax)]) # all after maximum

  # Density-based threshold (V4):
  # First minimum after the biggest maximum   (adjust=2)
  densTrh <- NULL
  if(length(minimumDens)>0) # && (!skipMinimumDens))
  {
    densTrh <- densCurve$x[min(minimumDens)]
    # Commented on V6
    # Only keep if it is a real inflextion point
    # (i.e. next max at least 5% of the global max)
    if(length(maximumsDens)>0)
    {
      nextMaxs <- maximumsDens[which(densCurve$x[maximumsDens] > densTrh)]
      if((max(densCurve$y[nextMaxs])/max(densCurve$y))<.05)
      {
        densTrh <- NULL
        # print(gSetName)
      }
    }
  }

  ## TO DO: Check special cases with many zeroes
  auc <- sort(auc)
  distrs <- list()
  distrs[["Global_k1"]] <- list(mu=c(meanAUC, NA), sigma=c(sdAUC, NA), x=auc)


  if("mixtools" %in% rownames(installed.packages()))
  {
    na <- capture.output(distrs[["k2"]] <-
                           tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=2, verb=FALSE),
                                    # With fast, if there are many zeroes, it fails quite often
                                    error = function(e) {
                                      return(NULL)
                                    }))

    na <- capture.output(distrs[["k3"]] <-
                           tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=3, verb=FALSE),
                                    error = function(e) {
                                      return(NULL)
                                    }))

    if(is.null(distrs[["k2"]]) && is.null(distrs[["k3"]]))
    {
      if(sum(auc==0)<(nCells*notPopPercent*.5))
        skipGlobal <- FALSE    # only if not too many zeroes??
      # qpois(1-(thrP/nCells), 1, log = FALSE)
      # plot(sort(rpois(auc, lambda=var(auc)), decreasing=TRUE))

      # commented V6 why was it here??
      # qPop <- quantile(auc, 1-smallestPopPercent)
      # if(sum(auc<qPop) >0)
      #   distrs[["k2"]] <- list(mu=c(mean(auc[auc<qPop]), NA),
      #                          sigma=c(sd(auc[auc<qPop]), NA),
      #                          lambda=c(1,NA), x=auc)
    }
    # if(!skipGlobal) print(gSetName) Warning?

    if(!is.null(distrs[["k2"]]))
    {
      compL <- which.min(distrs[["k2"]][["mu"]])
      compR <- which.max(distrs[["k2"]][["mu"]])
      ### Check distributions
      # Second distribution is "taller" than first one
      height1 <- .4/distrs[["k2"]][["sigma"]][compL]*
        distrs[["k2"]][["lambda"]][compL]
      height2 <- .4/distrs[["k2"]][["sigma"]][compR]*
        distrs[["k2"]][["lambda"]][compR]
      taller <- height1 < height2
      # Use global distr:
      # Mean of the global distr is included within the SD of the first
      # & Both means are included within the mean+SD of the Global distribution
      globalInclInFirst <-
        (distrs[["Global_k1"]]$mu[1] <
           (distrs[["k2"]][["mu"]][compL]+(1.5*distrs[["k2"]][["sigma"]][compL])))
      includedInGlobal <-
        ((distrs[["k2"]][["mu"]][compL] >
            (distrs[["Global_k1"]]$mu[1]-distrs[["Global_k1"]]$sigma[1])) &&
           (distrs[["k2"]][["mu"]][compR] <
              (distrs[["Global_k1"]]$mu[1]+distrs[["Global_k1"]]$sigma[1])))
      if(taller || (globalInclInFirst && includedInGlobal))
      {
        skipGlobal <- FALSE

        if(globalInclInFirst && includedInGlobal)
          commentMsg <- paste(commentMsg,
                              "The global distribution overlaps the partial distributions. ")
        if(taller && !includedInGlobal)
          commentMsg <- paste(commentMsg, "The right distribution is taller. ")
      }
    }
  }else{
    warning("Package 'mixtools' is not available to calculate the sub-distributions.")
  }

  glProb <- 1-(thrP/nCells + smallestPopPercent)   ## CORRECT?!?!
  aucThrs["Global_k1"] <- qnorm(glProb,# qnorm(1-(thrP/nCells),
                                mean=distrs[["Global_k1"]][["mu"]][1],
                                sd=distrs[["Global_k1"]][["sigma"]][1])
  if(!is.null(distrs[["k2"]]))
  {
    k2_L <- which.min(distrs[["k2"]][["mu"]]) # (sometimes the indexes are shifted)
    aucThrs["L_k2"] <- qnorm(1-(thrP/nCells),
                             mean=distrs[["k2"]][["mu"]][k2_L],
                             sd=distrs[["k2"]][["sigma"]][k2_L])
  }

  if(!is.null(distrs[["k3"]]))
  {
    k3_R <- which.max(distrs[["k3"]][["mu"]]) # R: right distribution
    k3_R_threshold <- qnorm(thrP,
                            mean=distrs[["k3"]][["mu"]][k3_R],
                            sd=distrs[["k3"]][["sigma"]][k3_R])
    if(k3_R_threshold > 0) aucThrs["R_k3"] <- k3_R_threshold
  }

  if(!is.null(densTrh))
  {
    aucThrs["minimumDens"] <- densTrh
  }

  aucThr <- aucThrs
  if(skipGlobal)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "Global_k1")]
  # TO DO: Decide when to merge with GLOBAL

  if(skipRed)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "L_k2")]
  # TO DO: Decide when to merge with GLOBAL

  aucThr <- aucThr[which.max(aucThr)] # to keep name
  if((length(aucThr)>0) && (names(aucThr) == "minimumDens"))
  {
    maximumsDens <- maximumsDens[which(densCurve$y[maximumsDens]>1)]
    if(length(maximumsDens) > 2)
    {
      tmp <- cbind(minimumDens[seq_len(length(maximumsDens)-1)],
                   maximumsDens[-1])
      FCs <- densCurve$y[tmp[,2]]/densCurve$y[tmp[,1]]
      FCs <- na.omit(FCs)
      if(any(FCs > 1.5))
        warning(gSetName,
                ":\tCheck the AUC histogram. ",
                "'minimumDens' was selected as the best threshold, ",
                "but there might be several distributions in the AUC.")
    }
  }

  if("minimumDens" %in% names(aucThrs))
    aucThr <- aucThrs["minimumDens"]
  if(length(aucThr)==0)
    aucThr <-  aucThrs[which.max(aucThrs)]
  aucThr
}


