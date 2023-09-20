
#' @importFrom ROCR prediction performance
#'
#' @export
evalPerform <- function(active_score,
                        col_data,
                        pos_cluster,
                        neg_cluster=NULL,
                        key_col="cluster",
                        type=c('KO','OE'),...){
  type <- match.arg(type)
  #---check cell names
  if(!identical(sort(names(active_score)),sort(row.names(col_data)))){
    stop("The names of 'active_score' is not identical with rownames of col_data!")
  }
  pos_cells <- row.names(col_data[col_data[[key_col]] %in% pos_cluster,])
  if(length(pos_cells)==0) stop("No cells found in your input cluster: ", pos_cluster)
  if(is.null(neg_cluster)){
    neg_cells <- setdiff(names(active_score), pos_cells)
  }else{
    neg_cells <- row.names(col_data[col_data[[key_col]] %in% neg_cluster,])
  }
  if(length(intersect(pos_cells, neg_cells))!=0) stop("No overlap cells allowed between ",
                                                      "positive and negative clusters!")
  all_cells <- c(pos_cells,neg_cells)
  #---
  pos_label <- rep(1,length(pos_cells))
  names(pos_label) <- pos_cells
  neg_label <- rep(-1,length(neg_cells))
  names(neg_label) <- neg_cells
  all_label <- c(pos_label, neg_label)[all_cells]
  if(type=='KO'){
    active_score <- -1*active_score
  }
  list(predict=active_score[all_cells],label=all_label)
}






