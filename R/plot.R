#' Plot EVPPI estimates as a simple dot plot
#'
#' This is a simple `ggplot2` plot and returns a `ggplot2` object, which
#' can be manipulated.  Or examine the source for this function to see how
#' to construct a similar plot that can be customised. 
#'
#' @param x Object returned from \code{\link{evppi}}. 
#'
#' @param order Order the plot with highest EVPPI values at the top.
#'
#' @param wtp Not implemented yet. Suggestions welcome for how. Multiple curves on same plot, for top few
#' or selected?
#'
#' @param ... Other arguments (currently unused).
#'
#' @return A `ggplot2` object. 
#'
#' @export
plot.evppi <- function(x, order=FALSE, wtp=NULL, ...){
  lower <- upper <- pars <- NULL
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
