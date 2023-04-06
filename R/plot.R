#' Plot EVPPI estimates
#' 
#' Plot EVPPI estimates as simple dot or curve plots.
#'
#' These plotting functions are intended for quick interactive exploration of
#' EVPPI results, so they deliberately have limited options for customising them.
#' 
#' For publication quality graphics, it is advised to use `ggplot2` by hand
#' on the data returned by `evppi`.   Examine the code for `plot_evppi_dots`
#' and `plot_evppi_curves` to see how these plots might be constructed. 
#'
#' @param x Object returned from \code{\link{evppi}}.
#'
#' @param type 
#'
#' \code{"dots"} for a dot plot of the EVPPI by parameter.  If \code{x} includes multiple
#' willingness-to-pay values for the same parameter, these are shown as multiple dots.
#'
#' \code{"curves"} for a plot of EVPPI against willingness-to-pay, with different
#' parameters distinguished as different curves.  This is only applicable if there
#' are multiple willingness-to-pay values included in \code{x}. 
#'
#' @param order For dot plots, order the plot with highest EVPPI values at the top.
#' 
#' @param top A positive integer. If specified, for example as `top=5` then only 
#' five parameters are included in the plot, those with the top five maximum EVPPI 
#' values by parameter. 
#' 
#' @param ... Other arguments (currently unused).
#'
#' @return A `ggplot2` object. 
#'
#' @export
plot.evppi <- function(x, type=NULL, order=FALSE, top=NULL, ...){
  lower <- upper <- pars <- NULL
  if (is.null(type)) {
    if (is.null(x$k) || (length(unique(x$k)) == 1))
      type <- "dots"
    else type <- "curves"
  }
  if (!is.null(top)){
    if (!is.numeric(top) || top<0) stop("`top` must be a positive number")
    max_by_par <- sapply(split(x, x$pars), function(y){max(y$evppi)})
    top_evppis <- rev(sort(max_by_par))[1:top]
    top_pars <- names(top_evppis)
    x <- x[x$par %in% top_pars, ,drop=FALSE]
  }
  if (type=="dots") { 
    plot_evppi_dots(x=x, order=order, ...)
  } else if (type=="curves"){
    if (is.null(x$k) || (length(unique(x$k))==1))
      stop("`curves` plots are only applicable with multiple willingness-to-pay values")
    plot_evppi_curves(x=x, ...)
  }
}

plot_evppi_dots <- function(x, order=FALSE, ...){
  pars <- evppi <- lower <- upper <- NULL
  if (order) x$pars <- reorder(x$pars, x$evppi)
  g <- ggplot2::ggplot(x, ggplot2::aes(y=pars, x=evppi)) +
    ggplot2::geom_point() + 
    ggplot2::ylab("") + 
    ggplot2::xlab("EVPPI") 
  if (!is.null(x$se)) {
    x$lower <- pmax(0, x$evppi - 2*x$se)
    x$upper <- pmax(0, x$evppi + 2*x$se)
    g <- g + ggplot2::geom_linerange(aes(xmin=lower, xmax=upper))
  } 
  g
}

plot_evppi_curves <- function(x, top=NULL, ...){
  k <- evppi <- pars <- NULL
  ggplot2::ggplot(x, ggplot2::aes(x=k, y=evppi, group=pars, col=pars)) + 
    ggplot2::geom_line() + 
    ggplot2::geom_point() +
    ggplot2::ylab("EVPPI") + 
    ggplot2::xlab("Willingness-to-pay") +
    ggplot2::labs(col="Parameters")
}
