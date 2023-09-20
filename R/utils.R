######################
# 1. transform data.frame to gene list
######################
#' Transform data.frame to gene list
#'
#' @param df A data.frame with TF and it's target information.
#' @param tf_col The column names of TFs.
#' @param target_col The column names of targets.
#'
#' @return A list of Transcription factor regulons. The names of list is TFs
#' @export
#'
#' @examples
#' data(dorothea_regulons)
#' mouse_regulons <- dfToList(df=dorothea_regulons$mm_regulons,
#' tf_col="tf", target_col="target")
#'
dfToList <- function(df, tf_col="tf", target_col="target"){
  if(!tf_col %in% colnames(df)) stop("The 'tf_col' is not in df!")
  if(!target_col %in% colnames(df)) stop("The 'target' is not in df!")
  regulons <- split(as.character(df[[target_col]]),as.character(df[[tf_col]]))
  lapply(regulons, unique)
}

######################
# 2. Get homologene from local database
######################
#' Get homologene from local database
#'
#' @param input_gene Input genes, which can be gene symbols or entrez gene ids.
#' @param input_tax Input gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}.
#' @param output_tax Output gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}
#' @param output_type Output types of genes. Default: gene_symbol.
#' @param db Database used in our homologene transformation.
#'
#' @return A data.frame include homologene gene information
#' @export
#'
#' @examples
getHomologGene <- function(input_gene=NULL,
                           input_tax="human",
                           output_tax="mouse",
                           output_type=c("gene_symbol","gene_id","both"),
                           db=metaTF::homologene_data){
  require(dplyr)
  output_type <- match.arg(output_type)
  #---names
  short_name <- c("human","mouse","rat","zebrafish","chicken")
  names(short_name) <- c("Homo sapiens","Mus musculus","Rattus norvegicus","Danio rerio","Gallus gallus")
  tax_tbl <- browserTaxonomy()
  tax_tbl$short_name <- short_name[tax_tbl$Taxonomy_Name]
  tax_tbl[is.na(tax_tbl)] <- ""
  #------check input id
  input_tax_id <- .check_tax_id(input_tax=input_tax, tax_tbl = tax_tbl)
  output_tax_id <- .check_tax_id(input_tax=output_tax, tax_tbl = tax_tbl)
  #------check input gene
  if(is.null(input_gene)) input_gene <- db$Gene_Symbol
  input_gene <- unique(input_gene)
  #---input
  input_tax_tab <- db %>%
    dplyr::filter(Taxonomy_ID %in% input_tax_id & (Gene_Symbol %in% input_gene | Gene_ID %in% input_gene)) %>%
    dplyr::select(HID, Gene_ID,Gene_Symbol)
  names(input_tax_tab)[2:3] <- paste(c("Gene_ID","Gene_Symbol"),input_tax,sep="_")
  #---output
  output_tax_tab <- db %>%
    dplyr::filter(Taxonomy_ID %in% output_tax_id & HID %in% input_tax_tab$HID) %>%
    dplyr::select(HID, Gene_ID, Gene_Symbol)
  names(output_tax_tab)[2:3] <- paste(c("Gene_ID","Gene_Symbol"),output_tax,sep="_")

  output_tab <- dplyr::inner_join(input_tax_tab,output_tax_tab,by="HID",keep=FALSE,na_matches="never") %>%
    dplyr::select(2,3,4,5)
  #------arrage output
  if(output_type=="gene_id"){
    output_tab <- output_tab[,grep("Gene_ID",colnames(output_tab))] %>% dplyr::distinct()
  }else if(output_type=="gene_symbol"){
    output_tab <- output_tab[,grep("Gene_Symbol",colnames(output_tab))] %>% dplyr::distinct()
  }
  output_tab$sortBy <- factor(output_tab[,1], levels = input_gene)
  output_tab <- dplyr::arrange(output_tab, sortBy)
  output_tab$sortBy <- NULL
  output_tab
}

######
#---check output
.check_tax_id <- function(input_tax, tax_tbl=browserTaxonomy()){
  if(is.numeric(input_tax)){
    if(!(input_tax %in% tax_tbl$Taxonomy_ID)) stop("No 'Taxonomy_ID' found for your input: ", input_tax)
  }else if(is.character(input_tax)){
    if(input_tax %in% tax_tbl$Taxonomy_Name){
      input_tax <- tax_tbl[tax_tbl$Taxonomy_Name==input_tax,"Taxonomy_ID"]
    }else if(input_tax %in% tax_tbl$short_name){
      input_tax <- tax_tbl[tax_tbl$short_name==input_tax,"Taxonomy_ID"]
    }else{
      if(as.numeric(input_tax) %in% tax_tbl$Taxonomy_ID){
        input_tax <- as.integer(input_tax)
      }else{
        stop("No 'Taxonomy_ID' found for your input: ", input_tax)
      }
    }
  }else{
    stop("The input_tax is a ", class(input_tax), ", it must be intiger or character!")
  }
  input_tax
}

#---browser
#' @export
browserTaxonomy <- function(db = metaTF::homologene_data){
  require(dplyr)
  as.data.frame(db %>% dplyr::select(Taxonomy_ID, Taxonomy_Name) %>%
                  dplyr::distinct(.))
}

#---homolog GRNs

#' Get gene regulatory network homologous
#'
#' @param df The input gene regulatory network used for homologous transformation.
#' please use \code{\link{browserTaxonomy()}}.
#' @param output_tax Output gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}.
#' @param output_type Output types of genes. Default: gene_symbol.
#' @param db Database used in our homologene transformation.
#'
#' @return A gene regulatory network of homologous genes.
#' @export
#'
#' @examples
#' hs_pantissue <- readRDS(system.file("./extdata/hs_pantissue_network.rds",package = "metaTF"))
#' mm_pantissue <- getHomologGRN(df = hs_pantissue, output_tax = "mouse")
getHomologGRN <- function(df, input_tax = "human",
                          output_tax = "mouse",
                          output_type = c("gene_symbol", "gene_id"),
                          db = metaTF::homologene_data){
  output_type <- match.arg(output_type)
  if(ncol(df)!=3) stop("The df must be a three column weighted network!")
  colnames(df) <- c("tf","target","weight")
  homolog_gene <- getHomologGene(input_gene = unique(c(as.character(df[[1]]), as.character(df[[2]]))), 
                                 output_tax = output_tax, output_type = output_type, db=db)
  colnames(homolog_gene) <- c("source_id","target_id")
  homolog_gene_net <- df %>% inner_join(x = ., y = homolog_gene, by=c("tf"="source_id")) %>% 
    inner_join(x = ., y = homolog_gene, by=c("target"="source_id")) %>% select(c(4,5,3))
  colnames(homolog_gene_net) <- c("tf","target","weight")
  homolog_gene_net
}

######
#---get homolog gene list
#' Title
#'
#' @param input_gene Input gene list
#' @param input_tax Input gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}
#' @param output_tax Output gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}
#' @param isListNameGene Whether the gene names is gene symbols?
#' @param ... Others parameters passed to \code{\link{getHomologGene})
#'
#' @return A gene list with homolog transformed list.
#' @export
#'
#' @examples
mapHomologList <- function(input_gene,
                           input_tax="human",
                           output_tax="mouse",
                           isListNameGene=FALSE, ...){
  input_gene <- as.list(input_gene)
  if(isListNameGene){
    all_genes <- unique(c(names(input_gene),unlist(input_gene)))
  }else{
    all_genes <- unique(unlist(input_gene))
  }
  gene_tbl <- getHomologGene(input_gene = all_genes,
                             input_tax=input_tax,
                             output_tax=output_tax,...)
  colnames(gene_tbl) <- c("human","mouse")
  output_list <- lapply(input_gene,function(x, gene_tbl){
    out <- gene_tbl %>% dplyr::filter(human %in% x)
    out <- out[['mouse']]
    out[!is.na(out)]},gene_tbl)
  if(isListNameGene){
    gene_tbl_filter <- gene_tbl %>% dplyr::filter(human %in% names(output_list))
    gene_tbl_filter$list_data <- output_list[gene_tbl_filter[["human"]]]
    output_list <- gene_tbl_filter$list_data
    names(output_list) <- gene_tbl_filter[["mouse"]]
  }
  output_list
}

######################
# 3. Obtaining TFRs from DoRoThEA
######################

#' Get TFRs From DoRoThEA
#'
#' @param input The regulons \code{data.frame} from \pkg{dorothea} package
#' @param levels The levels of data used as regulons.
#' @param min_target_num Miniumn target gene bumber.
#' @param max_target_num Maxiumn target gene bumber.
#'
#' @export
tfrFromDoRoThEA <- function(input,
                            levels = c("A","B","C","D","E"),
                            min_target_num = 5,
                            max_target_num = 1000){
  levels <- match.arg(levels, several.ok = TRUE)
  input <- as.data.frame(input[input$confidence %in% levels,])
  input$tf <- as.factor(input$tf)
  regulons <- split(input$target, input$tf)
  regulon_gene_num <- sapply(regulons, function(x){
    length(x[!is.na(x)])>=min_target_num && length(x[!is.na(x)])<=max_target_num})
  regulons[regulon_gene_num]
}

######################
# 4. re-level the cell-types
######################
setGeneric("reLevel", function(x, ...) standardGeneric("reLevel"))

#' Reset the levels of column data from object
#' 
#' Reset the levels of column data from object
#' 
#' @param x The SingleCellExperiment object.
#' @param colable The column in colData(x) used for re-order.
#' @param level_order The new orders.
#'
#' @return An object
#' @rdname reLevel
#' @export
setMethod("reLevel", "SingleCellExperiment", function(x,
                                                      colable,
                                                      level_order){
  #---relevel itself
  col_data <- SummarizedExperiment::colData(x)
  if(colable %in% colnames(col_data)){
    col_data[[colable]] <- factor(col_data[[colable]], levels = level_order)
    SummarizedExperiment::colData(x) <- col_data
  }else{
    warning("The colable: ",colable,"is not found in your object!\nPlease use colnames(colData(x)) to check it!")
  }
  #---relevel alt_exps
  alt_exp_name <- SingleCellExperiment::altExpNames(x)
  if(length(alt_exp_name)>0){
    alt_exps <- altExps(x)
    alt_exps <- lapply(alt_exps, function(y, colable, level_order){
      reLevel(x = y, colable=colable, level_order=level_order)
    },colable=colable, level_order=level_order)
    SingleCellExperiment::altExps(x) <- alt_exps
  }
  return(x)
})

setMethod("reLevel", "Seurat", function(x,
                                        colable,
                                        level_order){
  #---relevel itself
  col_data <- x@meta.data
  if(colable %in% colnames(col_data)){
    col_data[[colable]] <- factor(col_data[[colable]], levels = level_order)
    x@meta.data <- col_data
  }else{
    warning("The colable: ",colable,"is not found in your object!\nPlease use colnames(colData(x)) to check it!")
  }
  return(x)
})


######################
#5. my progress bar
######################
#' @importFrom progress progress_bar
#' @export
myPbs <- function(total){
  progress::progress_bar$new(format = "  |:bar| :percent ~:eta", complete='+',
                             incomplete=' ', current='+', total = total,
                             clear = FALSE, width= 64)
}
#' @importFrom parallel makeCluster
#' @export
mkCluster <- function(ncores=1){
  if(.Platform$OS.type=="unix"){
    cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores), type = "FORK")
  }else{
    cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores),type = "PSOCK")
  }
  cl
}



