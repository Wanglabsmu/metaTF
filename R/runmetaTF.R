
#' Runing ATFR in one function
#'
#' @param expression_data expression matrix
#' @param col_data sample annotation matrix
#' @param run_steps Run from this step
#' @param organism species
#' @param use_regulon regulon source
#' @param grn_method method used for regulon calculation
#' @param activity_method method for regulon activity
#' @param ncores number of cores
#' @param file_out path to output files.
#' @param ... others parameters
#'
#' @return
#' @export
#'
#' @examples
runATFR <- function(expression_data,
                    col_data=NULL,
                    run_steps=c(0,1,2,3,4),
                    organism = c("human","mouse"),
                    use_regulon = c("dorothea","cisTarget"),
                    grn_method = c("genie3", "sincerities", "sclink", "pcor", "spearman"),
                    activity_method = c("viper","aucell","gsva","ssgsea","zscore", "plage"),
                    ncores=1,
                    file_out="./", ...){
  st <- Sys.time()
  run_steps <- match.arg(run_steps)
  use_regulon <- match.arg(use_regulon)
  organism <- match.arg(organism)
  grn_method <- match.arg(grn_method)
  activity_method <- match.arg(activity_method)
  if(use_regulon == "dorothea"){
    if(organism=="human"){
      regulons <- dfToList(df = metaTF::dorothea_regulons$hs_regulons, tf_col = "tf",target_col = "target")
    }else if(organism=="mouse"){
      regulons <- dfToList(df = metaTF::dorothea_regulons$mm_regulons, tf_col = "tf",target_col = "target")
    }
  }else if(use_regulon == "cisTarget"){

  }
  #---output files
  file_out <- file.path(file_out,"output")
  dir.create(file_out, recursive = TRUE)

  #---step by step
  init_out <- file.path(file_out,"init")
  dir.create(init_out,recursive = TRUE)

  message("[1] Loading data...")
  if(run_steps==1) atfr <- readRDS(file.path(init_out,"01_init_object.rds"))
  atfr <- metaTF(exp_data = expression_data, col_data = col_data, regulons = regulons, ...)
  saveRDS(atfr,file.path(init_out,"01_init_object.rds"))

  message("[2] inferring GRNs ...")
  if(run_steps==2) atfr <- readRDS(file.path(init_out,"02_inferGRNs.rds"))
  atfr <- inferGRNs(x = atfr, method=grn_method, ncores=ncores, ...)
  saveRDS(atfr,file.path(init_out,"02_inferGRNs.rds"))

  message("[3] Filtering regulons ...")
  if(run_steps==3) atfr <- readRDS(file.path(init_out,"03_filterRegulons.rds"))
  atfr <- filterRegulons(x=atfr, ncores=ncores, ...)
  saveRDS(atfr,file.path(init_out,"03_filterRegulons.rds"))

  message("[4] regulonActivity ...")
  if(run_steps==4) atfr <- readRDS(file.path(init_out,"04_regulonActivity.rds"))
  atfr <- regulonActivity(atfr, method=activity_method,ncores=ncores,...)
  saveRDS(atfr,file.path(init_out,"04_regulonActivity.rds"))

  active_out <- file.path(file_out,"active_mat")
  dir.create(active_out,recursive = TRUE)
  out_mat <- assay(altExp(atfr,"atfr"))
  write.table(out_mat,file.path(active_out,paste(use_regulon,grn_method,activity_method,"active_mat.txt",sep="_")),sep="\t",quote=F)
  saveRDS(atfr,file.path(active_out,paste(use_regulon,grn_method,activity_method,"obj.rds",sep="_")))

  time_diff <- Sys.time()-st
  time_diff <- as.data.frame(paste(use_regulon,grn_method,activity_method,time_diff,sep="\t"))
  write.table(time_diff,path(active_out,"time_cost.txt"),sep="\t",quote=F)
  Sys.time()-st
}
