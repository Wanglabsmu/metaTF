#' Estimating the jaccard similarity between two gene lists
#'
#' This method estimate the similarities between two gene lists by using
#' Jaccard/Tanimoto similarity coefficients.
#'
#' @param tfr_list A list of transcription factor regulons. Names of the list 
#' were TF genes, while the value of each list is its target genes.
#' @param pathway_list A list of pathway gene sets. Names of the list could be 
#' pathway names, while the value of each pathway is genes.
#' @param method Method used compute p-values in \pkg{jaccard} package. The asymptotic
#' is the fastest, while the bootstrap is better but is a little bit slow. The exact
#' solution is slow for a moderately large vector. Default: asymptotic.
#' @param adjust Adjust methods for p-value. Default: "BH".
#' @param fdr The cutoff for significantly overlap terms. Default: 0.05.
#' @param collapse_set Collapse enriched terms with high similarities.
#' @param min_overlap The minimum enriched genes in pathway list. Default: 5.
#' @param max_cluster The maximum clusters used for merge enriched pathways. 
#' Default: auto, \code{sqrt(nrow(enrich_tab))}.
#' @param nperm Number of permutations. Used only for bootstrap method. 
#' Default: 1000.
#' @param ncores Number of cores used for parallel calculation. Default: 1.
#' @param verbose Whether to output information while running?
#'
#' @return a list of similarity test \code{data.frame} or a \code{matrix} of
#' normalized similarity scores.
#' @importFrom jaccard jaccard.test jaccard
#' @importFrom pbapply pblapply
#' @importFrom  dplyr tibble
#' @import data.table
#' @importFrom factoextra eclust
#'
#' @export
calcuModuleSimilar <- function(tfr_list,
                               pathway_list,
                               method = c("asymptotic", "bootstrap", "mca", "exact"),
                               adjust = c("BH", "bonferroni", "fdr"),
                               fdr=0.05,
                               collapse_set=TRUE,
                               min_overlap = 5,
                               max_cluster="auto",
                               nperm=1000,
                               ncores=1,
                               verbose=interactive()){
  method <- match.arg(method)
  adjust <- match.arg(adjust)
  if(!is.list(tfr_list) || !is.list(pathway_list)){
    stop("The 'tfr_list' and the 'pathway_list' must be list!")
  }
  set.seed(10086)
  #------check args
  total_genes <- union(unique(unlist(tfr_list)),unique(unlist(pathway_list)))
  tfr_name <- names(tfr_list)
  pathway_name <- names(pathway_list)
  if(is.null(tfr_name) || is.null(pathway_name)){
    stop("The names of 'tfr_list' and 'pathway_list' must not be NULL!")
  }
  if(verbose) message("[1] Perform similarity analysis using ",length(tfr_name),
                      " regulons with ", length(pathway_name),
                      " pathways. Using method: ", method)
  #---generate binary matrix for jaccard similarity analysis
  binary_mat <- matrix(0L,nrow=length(total_genes),ncol=length(c(tfr_name,pathway_name)),
                       dimnames = list(total_genes, c(tfr_name,pathway_name)))
  for(i in tfr_name){
    genes <- tfr_list[[i]]
    binary_mat[genes,i] <- 1L
  }
  for(j in pathway_name){
    genes <- pathway_list[[j]]
    binary_mat[genes,j] <- 1L
  }
  #---jaccard test
  #------par
  cl <- mkCluster(ncores = ncores)
  parallel::clusterEvalQ(cl,library(pbapply))
  parallel::clusterEvalQ(cl,library(jaccard))
  jaccard_list <- pbapply::pblapply(tfr_name,function(x,
                                                      pathway_name,
                                                      binary_mat,
                                                      method = c("bootstrap","asymptotic","mca","exact"),
                                                      adjust = c("BH","bonferroni","fdr"),
                                                      nperm=1000){
    method <- match.arg(method)
    adjust <- match.arg(adjust)
    jaccard_res <- apply(binary_mat[,pathway_name,drop=FALSE],2,function(p){
      if(method=="bootstrap"){
        test_res <- suppressMessages(
          jaccard::jaccard.test(x=p,
                                y= binary_mat[names(p),x,drop=TRUE],
                                method = method, B=nperm,seed=10086))
      }else{
        test_res <- jaccard::jaccard.test(x=p,
                                          y= binary_mat[names(p),x,drop=TRUE],
                                          method = method)
      }
      jaccard_similarity <- jaccard::jaccard(x = p, y = binary_mat[names(p),x,drop=TRUE], center = FALSE)
      c(similarity=jaccard_similarity,pval=test_res$pvalue)
    })
    jaccard_res <- as.data.frame(t(jaccard_res))
    jaccard_res$pathway <- row.names(jaccard_res)
    jaccard_res$padj <- p.adjust(jaccard_res$pval, method = adjust)
    jaccard_res[,c("pathway","similarity","pval","padj")]
  }, pathway_name=pathway_name, binary_mat=binary_mat, method=method, nperm=nperm, cl=cl)
  names(jaccard_list) <- tfr_name
  parallel::stopCluster(cl)
  #---collapse similar modules
  if(collapse_set){
    if(verbose) message("[2] Mergeing similar modules...")
    pb <- myPbs(total = length(jaccard_list))
    jaccard_list <- purrr::map2(
      .x = jaccard_list,
      .y = tfr_list[tfr_name],
      .f = function(tbl, tfr, pathway_list, min_overlap, fdr, max_cluster){
        pb$tick()
        pathway_list <- pathway_list[as.character(tbl$pathway)]
        tbl$overlap_targets <- lapply(pathway_list,function(x,tfr){
          intersect(x,tfr)},tfr=tfr)
        tbl <- tbl[sapply(tbl$overlap_targets,length)>=min_overlap,]
        tbl <- tbl[tbl$padj < fdr, ]
        tbl <- collapse_module(enrich_tab = tbl, 
                               key_ids = "pathway",
                               module_name = "overlap_targets",
                               group_by = "padj",
                               fdr = fdr, 
                               min_cluster = 5,
                               max_cluster=max_cluster,
                               verbose=FALSE)
        tbl
      }, pathway_list=pathway_list, min_overlap=min_overlap, fdr=fdr, max_cluster=max_cluster)
  }
  jaccard_list <- data.table::rbindlist(l = jaccard_list, fill=TRUE, idcol = "regulon")
  dplyr::tibble(jaccard_list)
}

#' Merging Enriched pathways based on term similarity
#' 
#' The \code{collapse_module} function is designed to clustering redundant terms 
#' into several modules based on the Jaccard similarity of enriched gene. The
#' term with minimum adjusted p-value in each module used to represent this 
#' module.
#' 
#' @param enrich_tab The enrichment result table. It must be a \code{\link{data.table}}
#' or \code{\link{tibble}}.
#' @param key_ids The key column to merge, usually enriched terms, it must be unique.
#' @param module_name The enriched genes of each term in key_ids column. It must 
#' be a list. So the input 'enrich_tab' could only be a \code{\link{data.table}}
#' or \code{\link{tibble}}.
#' @param group_by The value used to represent the enriched term in each module.
#' Default: padj.
#' @param fdr The cutoff used filtering significantly enriched terms. Default: 1, 
#' no filtering.
#' @param similarity_method The method used to calculate the similarity between
#' enriched terms.
#' @param min_cluster The minimum number passed to cluster terms, terms small 
#' than this value will not be clustered. Default: 5. 
#' @param max_cluster The maximum number of cluster after clustering. Default: 
#' "auto", \code{sqrt(nrow(enrich_tab))}.
#' @param verbose Whether to output information while running?
#' @param ... Other parameters passed to \code{\link{eclust()}}
#'
#' @export
#' @import data.table
#' @importFrom factoextra eclust
#' 
collapse_module <- function(enrich_tab,
                            key_ids = "pathway",
                            module_name = "overlap_targets",
                            group_by = c("padj","similarity","pval"),
                            fdr = 1L,
                            similarity_method = c("Jaccard","Ochiai"),
                            min_cluster = 5,
                            max_cluster = "auto",
                            verbose = interactive(), ...){
  #------check args
  similarity_method <- match.arg(similarity_method)
  group_by <- match.arg(group_by)
  enrich_tab <- data.table::as.data.table(enrich_tab)
  fgsea_out <- enrich_tab[padj < fdr,]
  if(any(duplicated(fgsea_out[[key_ids]]))){
    stop("The 'key_ids' can not have duplicated rows in 'fgsea_out'.")
  }
  #---remove duplicated values in fgsea output.
  fgsea_out <- fgsea_out[!duplicated(fgsea_out[[module_name]]),]
  if(max_cluster=="auto"){
    max_cluster <- ceiling(sqrt(nrow(fgsea_out)))
  }else if(is.numeric(max_cluster)){
    max_cluster <- floor(min(nrow(fgsea_out),max_cluster))
  }else{
    stop("The max_cluster should be set to 'auto' or a numeric!")
  }
  if(nrow(fgsea_out)==0) {
    return(fgsea_out[,`:=`(group=integer(),firstInGroup=integer())])
  }else if(nrow(fgsea_out) <= min_cluster){
    fgsea_out[, group:=rep(1L,nrow(fgsea_out))]
    if(group_by == "similarity"){
      fgsea_out[, firstInGroup:=ifelse(similarity == max(abs(similarity)),1,0), by=group]
    }else if(group_by == "padj"){
      fgsea_out[, firstInGroup:=ifelse(padj == min(padj),1,0), by=group]
    }else if(group_by == "pval"){
      fgsea_out[, firstInGroup:=ifelse(pval == min(pval),1,0), by=group]
    }
  }else{
    #------generating gene list
    leading_edge_list <- lapply(fgsea_out[[module_name]],function(x){x[!is.na(x)]})
    names(leading_edge_list) <- fgsea_out[[key_ids]]
    leading_edge_num <- sapply(leading_edge_list,length)
    leading_edge_list <- leading_edge_list[which(leading_edge_num>0)]
    if(any(which(leading_edge_num==0))) warning("Removing gene sets without valid genes.")
    #------calculating the similarity score
    if(similarity_method == "Jaccard"){
      leading_edge_simi_mat <- sapply(leading_edge_list,function(x){
        sapply(leading_edge_list,function(y){length(intersect(x,y))/length(union(x,y))})})
    }else if(similarity_method == "Ochiai"){
      leading_edge_simi_mat <- sapply(leading_edge_list,function(x){
        sapply(leading_edge_list,function(y){length(intersect(x,y))/sqrt(length(x)*length(y))})})
    }
    k_max <- min(nrow(leading_edge_simi_mat), max_cluster+1)-1
    leading_edge_eclust <- factoextra::eclust(x = leading_edge_simi_mat,
                                              FUNcluster = "hclust", k=NULL,
                                              graph=FALSE, k.max=k_max,
                                              stand=FALSE, verbose=verbose, ...)
    fgsea_out[, group:=leading_edge_eclust$cluster[fgsea_out[[key_ids]]]]
    if(group_by == "similarity"){
      fgsea_out[, firstInGroup:=ifelse(similarity == max(abs(similarity)),1,0), by=group]
    }else if(group_by == "padj"){
      fgsea_out[, firstInGroup:=ifelse(padj == min(padj),1,0), by=group]
    }else if(group_by == "pval"){
      fgsea_out[, firstInGroup:=ifelse(pval == min(pval),1,0), by=group]
    }
  }
  fgsea_out
}




