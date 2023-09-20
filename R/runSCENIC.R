setGeneric("runSCENIC", function(x, ...) standardGeneric("runSCENIC"))
#'
#' Running simplified version of SCENIC
#' 
#' This method aims to obtain the TF regulatory active matrix by running the 
#' RcisTarget and AUCell from SCENIC.
#' 
#' @param x The expression matrix with normalized values where rows were
#' genes and columns were cells. It can be a matrix or a 
#' \code{\link{SingleCellExperiment}} object with such expression matrix.
#' @param grn_tab The weighted TF regulatory matrix where rows were TF genes and 
#' columns were target genes.
#' @param motif_ranking The motif ranking object including \code{.feather} files 
#' import from \pkg{RcisTarget} database using \code{\link{importRankings()}} 
#' method.
#' @param motif_annotation The motif annotation database containing the 
#' annotations of the motif to transcription factors.
#' @param use_regulator Character. The TFs used for regulon analysis.
#' @param top_fraction The percentage of top ranked targets used for initial 
#' regulon analysis. Default: 0.05.
#' @param nesThreshold The normalized enrichment score for motif enrichment 
#' result. Default: 0, all enriched result.
#' @param minSize The minimum number of genes in one TF regulon.
#' @param maxSize The maximum number of genes in one TF regulon.
#' @param ncores Number of cores to use. Default: 1.
#' @param verbose Print message while runing?
#' 
#' @return A matrix or object contains TF regulon activity in each single-cell.
#' 
#' @rdname runSCENIC
#' @export
setMethod("runSCENIC", "matrix", function(x,
                                          motif_ranking,
                                          motif_annotation,
                                          use_regulator=NULL,
                                          grn_tab=NULL,
                                          top_fraction=0.05,
                                          nesThreshold=0L,
                                          minSize = 5L,
                                          maxSize = 500L,
                                          ncores=1,
                                          verbose=interactive()){
  exprMatr <- as.matrix(x)
  #1. get gene regulatory network with genie3
  if(!is.numeric(exprMatr) || !is.matrix(exprMatr)) stop("The exprMatr must be numeric matrix!")
  if(is.null(use_regulator)){
    use_regulator <- unique(motif_annotation[['TF']])
  }
  if(is.null(grn_tab)){
    if(verbose) message("No Gene Regulatory Network input, Runing GENIE3...")
    grn_tab <- GENIE3::GENIE3(exprMatrix = exprMatr, regulators = use_regulator, nCores = ncores)
  }else if(!is.numeric(grn_tab) || !is.matrix(grn_tab)){
    stop("The grn_tab must be numeric matrix!")
  }
  grn_tab <- grn_tab[intersect(row.names(grn_tab),use_regulator),,drop=FALSE]
  if(nrow(grn_tab)==0) stop("No regulators found in your GRN table!")
  active_mat <- runSCENICMethod(exprMatr=exprMatr,
                                grn_tab=grn_tab,
                                motif_ranking=motif_ranking,
                                motif_annotation=motif_annotation,
                                top_fraction=top_fraction,
                                nesThreshold=nesThreshold,
                                minSize= minSize,
                                maxSize=maxSize,
                                verbose=verbose)
  active_mat
})

#' @rdname runSCENIC
#' @export
setMethod("runSCENIC", "SingleCellExperiment", function(x,
                                                        motif_ranking,
                                                        motif_annotation,
                                                        use_assay=1L,
                                                        use_grn="genie3",
                                                        use_regulator=NULL,
                                                        top_fraction=0.05,
                                                        nesThreshold=0,
                                                        minSize=5L,
                                                        maxSize=500L,
                                                        verbose=interactive()){
  exprMatr <- assay(x,use_assay)
  grn_tab <- GRN(x,use_grn)
  if(is.null(grn_tab)) stop("No GRNs found in your object x, please run GENIE3",
                            " using inferGRNs() first!")
  if(is.null(use_regulator)){
    use_regulator <- unique(motif_annotation[['TF']])
  }
  grn_tab <- grn_tab[intersect(row.names(grn_tab),use_regulator),,drop=FALSE]
  if(nrow(grn_tab)==0) stop("No regulators found in your GRN table!")
  active_mat <- runSCENICMethod(exprMatr=exprMatr,
                                grn_tab=grn_tab,
                                motif_ranking=motif_ranking,
                                motif_annotation=motif_annotation,
                                top_fraction=top_fraction,
                                nesThreshold=nesThreshold,
                                minSize= minSize,
                                maxSize=maxSize,
                                verbose=verbose)
  sce <- SingleCellExperiment::SingleCellExperiment(as.matrix(active_mat))
  SummarizedExperiment::assayNames(sce) <- "auc"
  SummarizedExperiment::colData(sce) <- SummarizedExperiment::colData(x)
  SingleCellExperiment::altExp(x, e = "scenic") <- sce
  x
})

#' Running simplified version of SCENIC
#' 
#' This method aims to obtain the TF regulatory active matrix by running the 
#' RcisTarget and AUCell from SCENIC.
#'
#' @param exprMatr The expression matrix with normalized values where rows were
#' genes and columns were cells.
#' @param grn_tab The weighted TF regulatory matrix where rows were TF genes and 
#' columns were target genes.
#' @param motif_ranking The motif ranking object including '.feather' files 
#' import from \pkg{RcisTarget} database using \code{\link{importRankings()}} 
#' method.
#' @param motif_annotation The motif annotation database containing the 
#' annotations of the motif to transcription factors.
#' @param top_fraction The percentage of top ranked targets used for initial 
#' regulon analysis. Default: 0.05.
#' @param nesThreshold The normalized enrichment score for motif enrichment result. 
#' Default: 0, all enriched result.
#' @param minSize The minimum number of genes in one TF regulon.
#' @param maxSize The maximum number of genes in one TF regulon.
#' @param verbose Print message while runing?
#'
#' @return A matrix contains TF regulon activity in each single-cell.
#' 
#' 
#' 
#' @importFrom dplyr group_by top_n `%>%`
#' @importFrom GENIE3 GENIE3 getLinkList
#' @importFrom RcisTarget cisTarget 
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @importFrom stringr str_replace_all str_split
#'
#' @examples
#' @export
runSCENICMethod <- function(exprMatr,
                            grn_tab,
                            motif_ranking,
                            motif_annotation,
                            top_fraction=0.05,
                            nesThreshold=0L,
                            minSize = 5L,
                            maxSize = 500L,
                            verbose=interactive()){
  #1.get grn cis-target
  grn_tab <- as.matrix(grn_tab)
  top_num <- round(top_fraction*ncol(grn_tab))
  grn_frm <- GENIE3::getLinkList(weightMatrix = grn_tab, threshold = 0)
  grn_frm <- grn_frm %>% dplyr::group_by(regulatoryGene) %>% dplyr::top_n(n = top_num, wt = weight)
  grn_list <- split(as.character(grn_frm$targetGene),as.character(grn_frm$regulatoryGene))
  if(verbose) message("[1] Runing RcisTarget...")
  suppressMessages(suppressWarnings(
    grn_cistarget <- RcisTarget::cisTarget(geneSets = grn_list, 
                                         motifRankings = motif_ranking, 
                                         motifAnnot = motif_annotation, 
                                         highlightTFs = names(grn_list),
                                         nesThreshold = nesThreshold)
    ))
  grn_cistarget <- grn_cistarget[TFinDB %in% c("**","*"),][enrichedGenes!="",]
  if(nrow(grn_cistarget)==0) stop("The cistarget result is NULL! Please check ",
                                  "the species used between your data and annotation!")
  #---filter tf genes
  filterTFlists <- function(str){
    str <- stringr::str_replace_all(str, "\\s\\(directAnnotation\\)\\.",";")
    str <- stringr::str_replace_all(str, "\\s\\(inferredBy_Orthology\\)\\.",";")
    str <- stringr::str_replace_all(str, "\\s\\(inferredBy_MotifSimilarity\\)\\.",";")
    str <- stringr::str_replace_all(str, "\\s\\(inferredBy_MotifSimilarity_n_Orthology\\)\\.",";")
    str <- stringr::str_replace_all(str, "\\;\\s?+$","")
    str <- stringr::str_split(str,pattern = "\\;\\s")
    str
  }
  grn_cistarget$TF_highConf <- filterTFlists(grn_cistarget$TF_highConf)
  grn_cistarget$TF_lowConf <- filterTFlists(grn_cistarget$TF_lowConf)
  grn_cistarget <- grn_cistarget[apply(grn_cistarget,1,function(x){
    x['geneSet'] %in% unlist(x['TF_highConf']) || x['geneSet'] %in% unlist(x['TF_lowConf'])
  }),]
  grn_cistarget$enrichedGenes <- stringr::str_split(grn_cistarget$enrichedGenes, pattern =";")
  grn_regulon <- split(grn_cistarget$enrichedGenes,as.character(grn_cistarget$geneSet))
  grn_regulon <- lapply(grn_regulon, function(x){unique(unlist(x))})
  grn_regulon_num <- sapply(grn_regulon, length)
  grn_regulon <- grn_regulon[names(grn_regulon_num[which(grn_regulon_num>=minSize)])]
  grn_regulon <- grn_regulon[names(grn_regulon_num[which(grn_regulon_num<=maxSize)])]
  #3. runing AUCell
  if(verbose) message("[2] Runing AUCell...")
  suppressMessages(suppressWarnings(
  cells_rankings <- AUCell::AUCell_buildRankings(exprMat = exprMatr,
                                                 nCores=1, 
                                                 plotStats=FALSE)
  ))
  suppressMessages(suppressWarnings(
  cells_AUC <- AUCell::AUCell_calcAUC(geneSets=grn_regulon, 
                                      rankings = cells_rankings,
                                      nCores=1,
                                      normAUC = TRUE)
  ))
  auc_value <- AUCell::getAUC(object = cells_AUC)
  auc_value
}

######################
# 4. Obtaining TFRs from cisTarget
######################

#' extract TF regulons from Rcistarget database
#'
#' @param motifRnking The \code{\link{Rcistarget}} motif ranking file
#' @param motif_anno The motif annotation file
#' @param cutoff The cutoff of motif enrichment.
#' @param use_TFs User defined TFs.
#' @param motif_level The levels of motif used for filtering.
#' @param levels The level of motif used.
#'
#' @return
#'
#' @examples
#' @export
tfrFromRcisTarget <- function(motifRnking,
                              motif_anno,
                              cutoff=0.01,
                              use_TFs = NULL,
                              motif_level=c("high", "low", "custom"),
                              levels = c("directAnnotation","inferredBy_Orthology",
                                         "inferredBy_MotifSimilarity",
                                         "inferredBy_MotifSimilarity_n_Orthology")){
  motif_level <- match.arg(motif_level)
  if(motif_level=="high"){
    levels <- c("directAnnotation","inferredBy_Orthology")
  }else if(motif_level=="low"){
    levels <- c("directAnnotation","inferredBy_Orthology","inferredBy_MotifSimilarity",
                "inferredBy_MotifSimilarity_n_Orthology")
  }else{
    levels <- match.arg(levels,several.ok = TRUE)
  }
  ranking_tbl <- RcisTarget::getRanking(object = motifRnking)
  motif_anno_tbl <- dplyr::as_tibble(motif_anno) %>% dplyr::filter(annotationSource %in% levels)
  if(!is.null(use_TFs)){
    use_TFs <- intersect(use_TFs,motif_anno_tbl$TF)
    if(length(use_TFs)==0) stop("The TFs you input is not found!")
    motif_anno_tbl <- motif_anno_tbl %>% dplyr::filter(TF %in% use_TFs)
  }
  overlap_motif_id <- intersect(ranking_tbl$features, motif_anno_tbl$motif)
  if(length(overlap_motif_id)==0) stop("The motifRnking and motif_anno is not matching!")
  ranking_tbl <- ranking_tbl %>% dplyr::filter(features %in% overlap_motif_id)
  motif_anno_tbl <- motif_anno_tbl %>% dplyr::filter(motif %in% overlap_motif_id)
  
  #---get top ranked targets for each motif
  ranking_motif <- as.matrix(ranking_tbl[,-1])
  row.names(ranking_motif) <- ranking_tbl[['features']]
  message("Extracting motif list from motif ranking...")
  ranking_motif_list <- pbapply::pbapply(ranking_motif,1,function(x, cutoff){
    list(names(head(sort(x),round(length(x)*cutoff))))
  },cutoff=cutoff)
  ranking_motif_list <- lapply(ranking_motif_list,unlist)
  #---get enriched motif ids for each TF
  tf_motif_list <- split(motif_anno_tbl$motif,motif_anno_tbl$TF)
  message("Mergeing TF regulon list...")
  tf_regulon_list <- pbapply::pblapply(tf_motif_list,function(x,motif_list){
    unique(unlist(motif_list[intersect(x,names(motif_list))]))
  },motif_list=ranking_motif_list)
  tf_regulon_list
}
