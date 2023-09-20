#Setter/Getter
setGeneric("RadvizEnrichPlot",function(x, ...) standardGeneric("RadvizEnrichPlot"))

#' RadvizEnrichPlot
#'
#'
#' @param x A \code{\link{data.frame}} or objects used for \code{\link{Radviz}} plot.
#' @param anchors Anchors used to plot.
#' @param color_by Groups used to label the color.
#' @param trim_ratio Whether to trim the top and bottom of each anchor which correspond to outliers.
#' @param optim_anchor Whether to optimize the anchor order?
#' @param recenter User defined start anchor.
#' @param text_color Text color, inherited from \code{\link{getRadviz}}.
#' @param text_size Text size, inherited from \code{\link{getRadviz}}.
#' @param return_rv Only return \code{\link{getRadviz}} object, not generate figures.
#' @param outline_circle Whether to add circle in figure?
#' @param plot_type Figure types to plot. Default: point.
#' @param point_size Label size of points in \code{\link{geom_point()}}.
#' @param point_shape Label shape of points in \code{\link{geom_point()}}. Default: 16.
#' @param alpha The point transparency. Default: 0.8.
#' @param ... Other parameters passed to \code{\link{getRadviz}}.
#' @rdname RadvizEnrichPlot
#' @export
setMethod("RadvizEnrichPlot","data.frame",function(x,
                                                   anchors = colnames(x),
                                                   color_by = NULL,
                                                   trim_ratio = 0.005,
                                                   norm_data = TRUE,
                                                   optim_anchor = TRUE,
                                                   recenter = NULL,
                                                   text_color = "orangered4",
                                                   text_size = 8,
                                                   return_rv = FALSE,
                                                   outline_circle = TRUE,
                                                   plot_type = c("point","density","hexagonal","buble"),
                                                   point_size = 0.5,
                                                   point_shape = 16,
                                                   alpha=1, ...){
  rv <- getRadviz(input_frm=x, anchors=anchors, trim_ratio=trim_ratio, 
                  norm_data=norm_data, optim_anchor=optim_anchor,
                  recenter=recenter, text_color=text_color, text_size=text_size)
  
  
  rv_plot <- do_RadvizEnrichPlot(rv = rv, color_by = color_by, 
                                 outline_circle=outline_circle,
                                 plot_type=plot_type, point_size=point_size,
                                 point_shape=point_shape, alpha=alpha)
  if(return_rv) rv_plot <- rv
  rv_plot
})


####
#RadvizFeaturePlot
setGeneric("RadvizFeaturePlot",function(x, ...) standardGeneric("RadvizFeaturePlot"))

#' RadvizFeaturePlot
#'
#'
#' @param x A \code{\link{data.frame}} or objects used for \code{\link{Radviz}} plot.
#' @param anchors Anchors used to plot.
#' @param use_feature Features used to split the screens.
#' @param return_rv Only return \code{\link{getRadviz}} object, not generate figures.
#' @param anchors Dimensions used to generate the map.
#' @param optim_anchor Whether to optimize the anchor order?
#' @param recenter User defined start anchor.
#' @param text_color Text color, inherited from \code{\link{do.radviz}}.
#' @param text_size Text size, inherited from \code{\link{do.radviz}}.
#' @param color_by Groups used to label the color.
#' @param outline_circle Whether to add circle in figure?
#' @param plot_type Figure types to plot.
#' @param text_size Text size.
#' @param label_size Label size of point in \code{\link{geom_point()}}.
#' @param label_shape Label shape of point in \code{\link{geom_point()}}.
#' @param ... Other parameters passed to \code{\link{getRadviz}}.
#' @rdname RadvizFeaturePlot
#' @export
setMethod("RadvizFeaturePlot","data.frame",function(x,
                                                    anchors = colnames(x),
                                                    use_feature,
                                                    return_rv = FALSE,
                                                    trim_ratio=0.005,
                                                    norm_data = TRUE,
                                                    optim_anchor = TRUE,
                                                    recenter = NULL,
                                                    text_color = "orangered4",
                                                    text_size = 8,
                                                    outline_circle=TRUE,
                                                    point_size=1,
                                                    point_alpha=0.8,
                                                    point_shape=16,
                                                    plot_type=c("point","hexagonal"),
                                                    colors=scale_color_gradient(low='grey80', high="dodgerblue4"),
                                                    ...){
  rv <- getRadviz(input_frm = x, anchors = anchors, trim_ratio = trim_ratio, 
                  norm_data = norm_data, optim_anchor = optim_anchor, 
                  recenter = recenter, text_color = text_color, text_size = text_size)
  rv_plot <- do_RadvizFeaturePlot(rv = rv, 
                                  use_feature = use_feature, 
                                  outline_circle = outline_circle, 
                                  point_size = point_size,
                                  point_alpha = point_alpha, 
                                  point_shape = point_shape,
                                  plot_type=plot_type,
                                  colors = colors,
                                  text_size = text_size)
  if(return_rv) rv_plot <- rv
  rv_plot
})

#######################
# getRadviz object
#######################
#' getRadviz object
#'
#'
#' @param input_frm Matrix used to generate the Radviz object.
#' @param anchors Dimensions used to generate the map.
#' @param trim_ratio Whether to trim the top and bottom of each anchor which correspond to outliers.
#' @param optim_anchor Whether to optimize the anchor order?
#' @param recenter User defined start anchor.
#' @param text_color Text color, inherited from \code{\link{do.radviz}}.
#' @param text_size Text size, inherited from \code{\link{do.radviz}}.
#' @param ... Other parameters passed to \code{\link{do.radviz}}.
#' @return getRadviz object
#' 
#' @export
#' @importFrom Radviz do.L make.S cosine do.optimRadviz get.optim recenter do.radviz
getRadviz <- function(input_frm,
                      anchors=colnames(input_frm),
                      trim_ratio=0.005,
                      norm_data = TRUE,
                      optim_anchor = TRUE,
                      recenter = NULL,
                      text_color = "orangered4",
                      text_size = 8,
                      verbose=interactive(), ...){
  anchors <- intersect(anchors,colnames(input_frm))
  if(!all(anchors %in% colnames(input_frm))){
    not_in <- setdiff(anchors, colnames(input_frm))
    stop("Anchors [",not_in,"] not found from columns of input_frm!")
  }
  exp_mat <- as.matrix(input_frm[,anchors])
  if(!is.numeric(exp_mat)){
    stop("The 'anchors' columns in 'input_frm' must be numeric!")
  }
  #---normalization
  trans <- function(vec, cutoff=trim_ratio){
    Radviz::do.L(vec,fun = function(x) quantile(x,c(cutoff,1-cutoff)))
  }
  if(norm_data){
    exp_mat <- t(apply(exp_mat,1,trans,cutoff=trim_ratio))
    colnames(exp_mat) <- anchors
    other_cols <-setdiff(colnames(input_frm), anchors)
    if(length(other_cols)>0){
      input_frm <- data.frame(exp_mat,
                              input_frm[,other_cols,drop=FALSE],
                              check.names=F)
    }else{
      input_frm <- exp_mat
    }
  }
  #---generate springs
  mat.S <- Radviz::make.S(anchors)
  if(optim_anchor){
    mat.sim <- Radviz::cosine(exp_mat)
    optim.mat <- Radviz::do.optimRadviz(mat.S,mat.sim,iter=10,n=100)
    mat.S <- Radviz::make.S(Radviz::get.optim(optim.mat))
    if(!is.null(recenter) && length(recenter)==1 && (recenter %in% anchors)){
      mat.S <- Radviz::recenter(mat.S, newc = recenter)
    }
  }
  Radviz::do.radviz(x = input_frm, 
                    springs=mat.S, 
                    trans = trans,
                    label.color = text_color,
                    label.size = text_size)
}

################################
#------feature plot
################################
#' do_RadvizEnrichPlot
#'
#' @param rv A \code{\link{getRadviz}} object generate from \code{\link{do.radviz}}.
#' @param color_by Groups used to label the color.
#' @param outline_circle Whether to add circle in figure?
#' @param plot_type Figure types to plot.
#' @param text_size Text size.
#' @param point_size Label size of point in \code{\link{geom_point()}}.
#' @param point_shape Label shape of point in \code{\link{geom_point()}}.
#' @param ... Other input parameters.
#'
#' @return a ggplot object
#'
#' @export
do_RadvizEnrichPlot <- function(rv,
                                color_by = NULL,
                                outline_circle=TRUE,
                                plot_type=c("point","density","hexagonal","buble"),
                                point_size=0.5,
                                point_shape=16,
                                alpha=1,
                                ...){
  plot_type <- match.arg(plot_type)
  rv_data_frm <- as.data.frame(rv$proj$data)
  if(outline_circle) rv <- addRadvizCircle(rv)
  if(!is.null(color_by)){
    color_by <- intersect(color_by,colnames(rv_data_frm))
    if(length(color_by)!=1) stop("The 'color_by' shold only have one element!")
  }
  p <- rv$proj
  if(plot_type=="point"){
    if(is.character(point_shape) && point_shape %in% colnames(rv_data_frm)){
      p <- p + geom_point(aes_string(color=color_by, shape=point_shape),
                          size=point_size,alpha=alpha)
    }else if(is.numeric(point_shape)){
      p <- p + geom_point(aes_string(color=color_by),shape=point_shape,
                          size=point_size, alpha=alpha)
    }
  }else if(plot_type=="density"){
    p <- Radviz::smoothRadviz(rv)+geom_point(shape='.',alpha=alpha)
  }else if(plot_type=="hexagonal"){
    p <- Radviz::hexplot(rv)
  }else if(plot_type=="buble"){
    p <- Radviz::bubbleRadviz(rv,group=color_by)
  }
  p
}

################################
#------feature plot
################################
#' do_RadvizFeaturePlot
#'
#'
#' @param rv A \code{\link{getRadviz}} object generate from \code{\link{do.radviz}}
#' @param use_feature Features used to split the screens.
#' @param outline_circle Whether to add circle in figure?
#' @param plot_type Figure types to plot.
#' @param colors Colors for gradient.
#' @param text_size Text size.
#'
#' @return a ggplot object
#'
#' @export
do_RadvizFeaturePlot <- function(rv,
                                 use_feature=NULL,
                                 outline_circle=TRUE,
                                 point_size=1,
                                 point_alpha=0.8,
                                 point_shape=16,
                                 plot_type=c("point","hexagonal"),
                                 colors=scale_color_gradient(low='grey80',high="dodgerblue4"),
                                 text_size=8){
  plot_type <- match.arg(plot_type)
  rv_data_frm <- as.data.frame(rv$proj$data)
  use.feature <- intersect(use_feature, colnames(rv_data_frm))
  if(length(use.feature)==0) stop("The 'use.features' in your input is not found!")
  feature_mat <- as.matrix(rv_data_frm[,use.feature])
  if(!is.numeric(feature_mat)) stop("All the values in your features should be numeric!")

  #---start plot
  if(outline_circle) rv <- addRadvizCircle(rv)
  p <- rv$proj
  plots <- lapply(use.feature,function(fea){
    if(plot_type=="point"){
      geoms <- geom_point(aes_string(color=fea), shape=point_shape, size=point_size, alpha=point_alpha)
    }else if(plot_type=="hexagonal"){
      geoms <- stat_summary_hex(aes_string(z=fea), fun=median)
    }
    p+geoms+colors+labs(title=fea)+
      theme(text = element_text(size=text_size, hjust = 1,color="black"),
            plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5))
  })
  patchwork::wrap_plots(plots, ncol = ceiling(sqrt(length(use_feature))))
}

#------add out.liner for circle
addRadvizCircle <- function(rv, color="black"){
  circle <- annotate("path", x=cos(seq(0,2*pi,length.out=100)),y=sin(seq(0,2*pi,length.out=100)),color=color)
  rv$proj <- rv$proj+circle
  rv
}




