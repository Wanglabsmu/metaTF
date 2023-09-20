#' SINCERITIES
#'
#' SINCERITIES is a novel computational method for inferring
#' gene regulatory network (GRN) from time-stamped cross-sectional
#' single-cell expression data. SINCERITIES_PLUS function is an extension of
#' the default version of SINCERITIES, designed for single cell datasets with
#' fewer than five time points.
#'
#' @param DATA a list containing the following information:\cr
#' \itemize{
#' \item \code{singleCELLdata}: list of length \emph{n}, where \emph{n} is the
#' number of capture time points. \code{DATA$singleCELLdata[[k]]} is a \emph{s_k}
#' by \emph{m} \code{matrix}(\code{data.frame}) containing observed expression
#' levels of \emph{m} genes in \emph{s_k} single cells.
#' \item \code{totDATA}: \emph{S} by \emph{m} matrix, where \emph{S} is the
#' total number of single cells (i.e., \emph{S=s_1+s_2+...+s_n} where \emph{n}
#' the number of capture time points) and \emph{m} is the number of genes.
#' \item \code{time}: vector of length \emph{n} containing the cell capture
#' time points or time-stamps.}
#' @param distance This parameter selects the distribution distance. 1 stands
#' for KS (Kolmogorov-Smirnov), 2 stands for CM (Cramer-von Mises), 3 stands
#' for AD (Anderson-Darling), 4 stands for Mean Expression Difference. Default: 1.
#' @param method This parameter selects the regularization regression strategy.
#' 1 stands for RIDGE, 2 stands for ELASTIC-NET with automatic detection of
#' optimal alpha parameter, 3 stands for LASSO, 4 stands for ELASTIC-NET with
#' manual selection of alpha parameter. Default: 1.
#' @param noDIAG This parameter selects whether the auto-regulatory edge is
#' inferred. 0 stands for GRN contains no auto-regulatory edge; 1 stands for
#' GRN contain auto-regulatory edge. Default: 0
#' @param SIGN This parameter selects whether the sign / mode of the gene
#' regulations is inferred. 0 stands for unsigned GRN; 1 stands for signed GRN.
#' Default: 1.\cr
#' \code{SINCERITIES} uses partial correlation analysis where a positive (negative)
#' correlation is taken as an indication of activation (repression).
#' @param CV_nfolds Defines a partition of the data into \code{CV_n_folds} disjoint
#' subsets for the cross validation. DEFAULT: 5.
#'
#' @return A list containing the following information:
#' \itemize{
#' \item adj_matrix: \emph{m} by \emph{m} matrix containing the weights of
#' regulatory edges. The larger \code{adj_matrix[i,j]} indicates higher confidence
#' that the corresponding edge exists (i.e., gene \emph{i} regulating gene \emph{j}).
#' \item DISTANCE_matrix: \emph{n-1} by m matrix containing the (normalized)
#' distribution distance (DD) computed during the network inference.}
#' @author \code{SINCERITIES} were created by Nan Papili Gao and R version implemented
#' by Ziyi Hua. Institute for Chemical and Bioengineering.ETH Zurich.
#' E-mail: nanp@ethz.ch.\cr
#' Puwen Tan Integrated the \code{SINCERITIES} into \pkg{metaTF}.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom ppcor pcor
#' @importFrom cvTools cvFolds
#' @importFrom kSamples ad.test
#'
#'
#------when time points > 5
#' @rdname SINCERITIES
#' @export
SINCERITIES <- function(DATA, distance=1, method=1, noDIAG=0, SIGN=1){
  if(!distance%in%c(1,2,3,4)){
    stop("Choose distance metric with 1,2 and 3: 1 for Kolmogorov-Smirnov(KM), 2 for Cramer-von Mises(CM),
         3 for Anderson-Darling(AD)")
  }
  if(!method%in%c(1,2,3,4)){
    stop("Choose regularization regression strategy with 1,2,3 and 4: 1 for RIGDE, 2 for ELASTIC NET with automatic
         detection of optimal alpha parameter, 3 for LASSO, 4 for ELASTIC NET with manual selection of alpha parameter")
  }
  if(!noDIAG%in%c(0,1)){
    stop("noDIAG should be either 0 or 1")
  }
  if(!SIGN%in%c(0,1)){
    stop("SIGN should be either 0 or 1")
  }
  #Initialization
  single_cell_data <- DATA$singleCELLdata
  time <- DATA$time
  numGENES <- DATA$numGENES
  num_time_points <- length(time)
  gene_names <- DATA$genes
  #Catch error
  if(num_time_points<5){
    stop('** DATA with a number of time points < 5. Please run SINCERITIES_CROSS_VALIDATION function **')
  }

  #Distribution Distance
  DISTANCE_matrix <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
  row_names_dist_mat <- paste(time[1:(num_time_points-1)],time[2:num_time_points],sep="_to_")
  dimnames(DISTANCE_matrix) <- list(row_names_dist_mat, gene_names)
  #add names to DISTANCE_matrix(by tpw, 2021-4-10)
  #h <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
  totalDATA <- single_cell_data[[1]]

  #cmtest2 <- dget("SINCERITIES functions/cmtest2.R")

  for (ti in 1:(num_time_points-1)) {
    totalDATA <- rbind(totalDATA,single_cell_data[[ti+1]])
    data_ti <- t(single_cell_data[[ti]])
    data_ti_plus1 <- t(single_cell_data[[ti+1]])

    for (gi in 1:numGENES) {
      p1 <- data_ti[gi,]
      p2 <- data_ti_plus1[gi,]
      if(distance==1){
        test.stat <- ks.test(p1,p2)
        DISTANCE_matrix[ti,gi] <- test.stat$statistic
      }else if(distance==2){
        DISTANCE_matrix[ti,gi] <- cmtest2(p1,p2)$CM_limiting_stat
      }else if(distance==3){
        test.stat <- kSamples::ad.test(p1,p2)
        DISTANCE_matrix[ti,gi] <- test.stat$ad[2,1]
      }else if(distance==4){
        DISTANCE_matrix[ti,gi] <- abs(mean(p1)-mean(p2))
      }
    }
  }

  #normalization
  deltaT <- replicate(dim(DISTANCE_matrix)[2],time[2:length(time)]-time[1:(length(time)-1)])
  DISTANCE_matrix_normed <- DISTANCE_matrix/deltaT

  #Generate Y and X_matrix for glmnet
  if(method==1){
    alphas <- 0
  }else if(method==2){
    alphas <- seq(0,1,0.1)
  }else if(method==3){
    alphas <- 1
  }else{
    input <- readline(' *** Please input manually the alpha values (between 0 and 1) separated by comma: ')
    alphas <- as.numeric(unlist(strsplit(input,',')))
  }
  DISTANCE_matrix <- DISTANCE_matrix_normed
  X_matrix <- DISTANCE_matrix[1:(num_time_points-2),]

  #LOOCV settings
  nfold <- dim(X_matrix)[1]
  foldid <- 1:nfold
  keep <- TRUE
  pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)
  dimnames(pred_lambda_min) <- list(gene_names, gene_names)
  #set gene names for adj matrix.
  lambda_res <- vector()
  alpha_res <- vector()

  for (gi in 1:numGENES) {

    lambda <-  vector()
    cvERROR <-  vector()
    beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))

    for (test in 1:length(alphas)) {
      Y_vector <- DISTANCE_matrix[2:(num_time_points-1),gi]
      if(noDIAG==1){
        CV_results <- glmnet::cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
                                        keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
      }else{
        CV_results <- glmnet::cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
                                        keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
      }
      lambda[test] <- CV_results$lambda.min
      cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
      #---new version do not support coef.cv.glmnet but replaced with coef()
      #---coef(CV_results, s = "lambda.min") (by Puwen Tan, 2021-4-9)
      coef.cv.glmnet=function(object,s=c("lambda.1se","lambda.min"),...){
        if(is.numeric(s))lambda=s
        else
          if(is.character(s)){
            s=match.arg(s)
            lambda=object[[s]]
          }
        else stop("Invalid form for s")
        coef(object$glmnet.fit,s=lambda,...)
      }
      #---end coef.cv.glmnet
      #coef.CV_results <- coef.cv.glmnet(CV_results,s='lambda.min')
      #---here we replace the coef.cv.glmnet with coef. (by Puwen Tan, 2021-4-9)
      coef.CV_results <- coef(CV_results,s='lambda.min')
      beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    }
    minIdx <- max(which(cvERROR==min(cvERROR)))
    lambda_res[gi] <- lambda[minIdx]
    alpha_res[gi] <- alphas[minIdx]
    pred_lambda_min[,gi] <- beta[,minIdx]

  }

  if(SIGN==1){
    parcorr_matrix <- ppcor::pcor(DATA$totDATA,method = 'spearman')$estimate
    pred_lambda_min <- pred_lambda_min*sign(parcorr_matrix)
  }

  result <- list(DISTANCE_matrix=DISTANCE_matrix,adj_matrix=pred_lambda_min)
  return(result)
}

#------SINCERITIES_PLUS
#' @rdname SINCERITIES
#' @export
SINCERITIES_PLUS <- function(DATA, noDIAG = 0, SIGN = 1, CV_nfolds = 5){
  #Initialization
  single_cell_data <- DATA$singleCELLdata
  time <- DATA$time
  numGENES <- DATA$numGENES
  num_time_points <- length(time)
  gene_names <- DATA$genes
  #Catch error
  if(num_time_points<3){
    stop('** The data must contain at least 3 time points **')
  }

  #K-fold CV
  library(cvTools)
  K <- CV_nfolds
  indices <- list()
  for(cv in 1:length(single_cell_data)){
    N <-  dim(single_cell_data[[cv]])[1]
    indices[[cv]] <- cvFolds(N,K = K)
  }

  error_CV <- array(0,dim = c(K,numGENES,100))

  for(cross in 1:K){
    data4test <- list()
    data4train <- list()
    for(cv in 1:length(single_cell_data)){
      test <- indices[[cv]]$subsets[indices[[cv]]$which==cross]
      train <- indices[[cv]]$subsets[indices[[cv]]$which!=cross]
      data4test[[cv]] <- single_cell_data[[cv]][test,]
      data4train[[cv]] <- single_cell_data[[cv]][train,]
    }

    #Distribution Distance
    DISTANCE_matrix_train <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
    row_names_dist_mat <- paste(time[1:(num_time_points-1)],time[2:num_time_points],sep="_to_")
    dimnames(DISTANCE_matrix_train) <- list(row_names_dist_mat, gene_names)
    #names added by tpw(2021-4-10)
    totalDATA <- data4train[[1]]

    for(ti in 1:(num_time_points-1)){
      totalDATA <- rbind(totalDATA,data4train[[ti+1]])
      data_ti <- t(data4train[[ti]])
      data_ti_plus1 <- t(data4train[[ti+1]])
      for(gi in 1:numGENES){
        p1 <- data_ti[gi,]
        p2 <- data_ti_plus1[gi,]
        test.stat <- ks.test(p1,p2)
        DISTANCE_matrix_train[ti,gi] <- test.stat$statistic
      }
    }

    #Normalization
    deltaT <- replicate(dim(DISTANCE_matrix_train)[2],time[2:length(time)]-time[1:(length(time)-1)])
    DISTANCE_matrix_train_normed <- DISTANCE_matrix_train/deltaT
    X_matrix <- DISTANCE_matrix_train_normed[1:(num_time_points-2),]

    #Distribution Distance Test
    DISTANCE_matrix_test <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
    totalDATA <- data4test[[1]]

    for(ti in 1:(num_time_points-1)){
      totalDATA <- rbind(totalDATA,data4test[[ti+1]])
      data_ti <- t(data4test[[ti]])
      data_ti_plus1 <- t(data4test[[ti+1]])
      for(gi in 1:numGENES){
        p1 <- data_ti[gi,]
        p2 <- data_ti_plus1[gi,]
        test.stat <- ks.test(p1,p2)
        DISTANCE_matrix_test[ti,gi] <- test.stat$statistic
      }
    }

    #Normalization
    deltaT <- replicate(dim(DISTANCE_matrix_test)[2],time[2:length(time)]-time[1:(length(time)-1)])
    DISTANCE_matrix_test_normed <- DISTANCE_matrix_test/deltaT
    X_matrix_test <- DISTANCE_matrix_test_normed[1:(num_time_points-2),]

    #Generate Y and X_matrix for glmnet
    alphas <- 0 #Ridge Regression
    lambdas <- 10^seq(-2,2,length.out = 100)

    for(gi in 1:numGENES){
      Y_vector <- DISTANCE_matrix_train_normed[2:(num_time_points-1),gi]
      if(noDIAG==1){
        CV_results <- glmnet::glmnet(X_matrix,Y_vector,alpha = alphas,exclude = gi,lambda = lambdas,
                                     lower.limits = 0, upper.limits = Inf)
      }else{
        CV_results <- glmnet::glmnet(X_matrix,Y_vector,alpha = alphas, lambda = lambdas,
                                     lower.limits = 0, upper.limits = Inf)
      }
      Y_vector_test <- DISTANCE_matrix_test_normed[2:(num_time_points-1),gi]

      for(lambdacount in 1:length(CV_results$lambda)){
        beta_lambda <- as.matrix(CV_results$beta)[,lambdacount]
        error_CV[cross,gi,lambdacount] <- sum((Y_vector_test - X_matrix_test%*%beta_lambda)^2)
      }
    }
  }

  mean_error_CV <- apply(error_CV,c(2,3),mean)
  standard_error_mean <- apply(error_CV,c(2,3),sd)/sqrt(K)

  #Lambda min
  min_mean_error_CV <- apply(mean_error_CV,1,min)
  idx_lambda_min <- apply(mean_error_CV,1,which.min)
  lambda_min <- CV_results$lambda[idx_lambda_min]

  #Lambda 1SE
  idx_lambda_1SE <- vector(length = numGENES)
  for(gi in 1:numGENES){
    min_plu_1SE <- mean_error_CV[gi,idx_lambda_min[gi]]+standard_error_mean[gi,idx_lambda_min[gi]]
    idx_lambda_1SE[gi] <- which(mean_error_CV[gi,]<=min_plu_1SE)[1]
  }
  lambda_1SE <- CV_results$lambda[idx_lambda_1SE]

  #SINCERITIES_final <- dget("SINCERITIES functions/SINCERITIES_final.R")
  pred_lambda <- SINCERITIES_final(single_cell_data,time,numGENES,num_time_points,lambdas,alphas,idx_lambda_min,noDIAG)
  dimnames(pred_lambda) <- list(gene_names,gene_names)
  if(SIGN==1){
    parcorr_matrix <- ppcor::pcor(DATA$totDATA,method = 'spearman')$estimate
    pred_lambda <- pred_lambda*sign(parcorr_matrix)
  }

  result <- list(DISTANCE_matrix=DISTANCE_matrix_train_normed,adj_matrix=pred_lambda)
  return(result)
}

#------cmtest2
cmtest2 <- function(x1,x2,alpha=0.05){
  if(alpha<0||alpha>1) stop("alpha should be between 0 and 1")

  sorted <- sort(c(x1,x2))
  binEdges <- c(c(min(sorted)-1,sorted),max(sorted)+1)
  binCounts1 <- c(hist(x1,breaks = binEdges, right = FALSE, plot = FALSE)$counts,0)
  binCounts2 <- c(hist(x2,breaks = binEdges, right = FALSE, plot = FALSE)$counts,0)
  sumCounts1 <- cumsum(binCounts1)/sum(binCounts1)
  sumCounts2 <- cumsum(binCounts2)/sum(binCounts2)
  sampleCDF1 <- sumCounts1[1:(length(sumCounts1)-1)]
  sampleCDF2 <- sumCounts2[1:(length(sumCounts2)-1)]
  N1 <- length(x1)
  N2 <- length(x2)
  N <- N1+N2

  #Compute test statistic of interest
  CMstatistic <- N1*N2/N^2*sum((sampleCDF1 - sampleCDF2)^2)

  #table of the limiting distribution
  z=c(0.00000,0.02480,0.02878,0.03177,0.03430,0.03656,0.03865,0.04061,0.04247,0.04427,0.04601,0.04772,0.04939,0.05103,0.05265,0.05426,0.05586,0.05746,0.05904,0.06063,0.06222,0.06381,0.06541,0.06702,0.06863,
      0.07025,0.07189,0.07354,0.07521,0.07690,0.07860,0.08032,0.08206,0.08383,0.08562,0.08744,0.08928,0.09115,0.09306,0.09499,0.09696,0.09896,0.10100,0.10308,0.10520,0.10736,0.10956,0.11182,0.11412,0.11647,
      0.11888,0.12134,0.12387,0.12646,0.12911,0.13183,0.13463,0.13751,0.14046,0.14350,0.14663,0.14986,0.15319,0.15663,0.16018,0.16385,0.16765,0.17159,0.17568,0.17992,0.18433,0.18892,0.19371,0.19870,0.20392,
      0.20939,0.21512,0.22114,0.22748,0.23417,0.24124,0.24874,0.25670,0.26520,0.27429,0.28406,0.29460,0.30603,0.31849,0.33217,0.34730,0.36421,0.38331,0.40520,0.43077,0.46136,0.49929,0.54885,0.61981,0.74346,
      1.16786)
  Pz=c(seq(0,0.99,0.01),0.999)

  #compute parameters of the statistic's distribution
  T_mean <- 1/6+1/6/N;
  T_var <- 1/45*(N+1)/N^2 * ( 4*N1*N2*N-3*(N1^2+N2^2)-2*N1*N2 ) / (4*N1*N2)
  #translate the T statistic into the limiting distribution
  CM_limiting_stat <- ( CMstatistic - T_mean ) / sqrt(45*T_var) + 1/6;
  #interpolate
  if(CM_limiting_stat>z[length(z)]) pValue <- 1
  else if(CM_limiting_stat<z[1]) pValue <- 0
  else pValue <- approx(z,Pz,CM_limiting_stat)$y

  H <- 1*(alpha>(1-pValue))

  return(list(H=H,CM_limiting_stat=CM_limiting_stat,pValue=pValue))
}

#------SINCERITIES_final
SINCERITIES_final <- function(single_cell_data,
                              time,
                              numGENES,
                              num_time_points,
                              lambdas,
                              alphas,
                              idx_lambda_min,
                              noDIAG){
  #Distribution Distance
  DISTANCE_matrix <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
  totalDATA <- single_cell_data[[1]]

  for(ti in 1:(num_time_points-1)){
    totalDATA <- rbind(totalDATA,single_cell_data[[ti+1]])
    data_ti <- t(single_cell_data[[ti]])
    data_ti_plus1 <- t(single_cell_data[[ti+1]])
    for(gi in 1:numGENES){
      p1 <- data_ti[gi,]
      p2 <- data_ti_plus1[gi,]
      test.stat <- ks.test(p1,p2)
      DISTANCE_matrix[ti,gi] <- test.stat$statistic
    }
  }

  #Normalization
  deltaT <- replicate(dim(DISTANCE_matrix)[2],time[2:length(time)]-time[1:(length(time)-1)])
  DISTANCE_matrix_normed <- DISTANCE_matrix/deltaT
  X_matrix <- DISTANCE_matrix_normed[1:(num_time_points-2),]

  #glmnet
  pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)
  for(gi in 1:numGENES){
    Y_vector <- DISTANCE_matrix_normed[2:(num_time_points-1),gi]
    if(noDIAG==1){
      CV_results <- glmnet::glmnet(X_matrix,Y_vector,alpha = alphas,exclude = gi,lambda = lambdas,
                                   lower.limits = 0, upper.limits = Inf)
    }else{
      CV_results <- glmnet::glmnet(X_matrix,Y_vector,alpha = alphas, lambda = lambdas,
                                   lower.limits = 0, upper.limits = Inf)
    }
    pred_lambda_min[,gi] <- as.matrix(CV_results$beta)[,idx_lambda_min[gi]]
  }

  return(pred_lambda_min)
}

#------final_ranked_predictions
final_ranked_predictions <- function(connectivityMATRIX,
                                     genes,
                                     SIGN=0,
                                     norm_interaction=FALSE,
                                     fileNAME=NULL,
                                     saveFile=FALSE){
  if(is.null(fileNAME)){
    fileNAME <- 'GRNprediction.txt'
  }

  numGENES <- length(genes)
  interactions <- as.vector(connectivityMATRIX)

  edges <- vector(length = length(interactions))
  if(SIGN==1){
    edges[which(interactions<0)] <- "repression"
    edges[which(interactions>0)] <- "activation"
    edges[which(interactions==0)] <- "no regulation"
  }else{
    edges[which(interactions==0)] <- "no regulation"
    edges[which(interactions!=0)] <- "activation/repression"
  }
  if(norm_interaction){
    interactions <- abs(interactions)
  }

  targetGENES <- as.vector(replicate(numGENES,genes))
  sourceGENES <- as.vector(t(replicate(numGENES,genes)))

  df <- data.frame(sourceGENES,targetGENES,interactions,edges)
  colnames(df) <- c('SourceGENES','TargetGENES','Interaction','Edges')
  df <- df[order(-df$Interaction),]
  row.names(df) <- 1:dim(df)[1]
  if(saveFile) write.csv(df,file = fileNAME,row.names = FALSE,quote = FALSE)

  return(df)
}
