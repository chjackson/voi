##' Check the fit of a regression model used to estimate EVPPI or EVSI
##'
##' Produces diagnostic plots and summaries of regression models used to estimate EVPPI or EVSI, 
##' mainly in order to check that the residuals have mean zero.
##'
##' @details For VoI estimation, the key thing we are looking for is that the residuals
##' have mean zero, hence that the mean of the model output is represented well by the
##' regression function of the model input parameters.  It should not matter if the 
##' variance of the residuals is non-constant, or non-normally distributed.
##'
##' Models produced with `method="gam"` are summarised using \code{\link{gam.check}}.
##'
##' Models produced `method="earth"` are summarised using \code{\link{plot.earth}}.
##' 
##' For any regression model, if `fitted()` and `residuals()` methods are defined for those models,
##' then a histogram of the residuals and a scatterplot of residuals against fitted values is produced.
##'
##' @param x Output from \code{\link{evppi}} or \code{\link{evsi}}. The argument \code{check=TRUE}
##' must have been used when calling \code{evppi} or \code{evsi}, to allow the regression model
##' objects from \code{gam} or \code{earth} to be preserved.  (This is not done by
##' default, since these objects can be large.).   \code{attr(x, "models")} contains these objects.
##' 
##' @param pars Parameter (or parameter group) whose EVPPI calculation is to be checked.
##' This should be in the \code{pars} component of the object returned by \code{\link{evppi}}.
##' Only relevant if \code{x} is the result of an \code{\link{evppi}} calculation.  By default,
##' the first calculation shown in \code{x} is checked.
##' 
##' @param n Sample size whose EVSI calculation is to be checked. 
##' This should be in the \code{n} component of the object returned by \code{\link{evsi}}.
##' Only relevant if \code{x} is the result of an \code{\link{evsi}} calculation.
##'
##' @param comparison Only relevant if there are more than two treatments in the decision model.
##' Different regression models are then used for the comparisons of different treatments
##' with the baseline treatment. 
##' \code{comparison} is an integer identifying which of these models is checked. 
##'
##' @param outcome \code{"costs"} or \code{"effects"}.  Only relevant if `outputs` was
##' in cost-effectiveness format when
##' calling \code{evppi} or \code{evsi}, hence different regressions are used for costs and
##' effects.  By default, \code{outcome="costs"} is used, so that the regression
##' for costs is checked.
##'
##' @param plot If \code{FALSE}, only numerical statistics are returned, and a plot is not made.
##'
##' @return Where possible, an appropriate statistic is returned that allows the regression
##' model to be compared with other regression models implemented using the same \code{method}
##' but with different assumptions.   For \code{method="gam"},
##' this is Akaike's information criterion (AIC).    
##' For \code{method="earth"}, this is the generalised cross-validation statistic
##' \code{gcv}.    Currently not implemented for other methods. 
##'
##' @examples
##' pars <- c("p_side_effects_t1", "p_side_effects_t2")
##' evtest <- evppi(chemo_nb, chemo_pars, pars=pars, check=TRUE)
##' evtest
##' check_regression(evtest)
##' 
##' ## with no interaction term 
##' evtest2 <- evppi(chemo_nb, chemo_pars, pars=pars, 
##'                 gam_formula="s(p_side_effects_t1)+s(p_side_effects_t2)",
##'                 check=TRUE)
##' evtest2
##' check_regression(evtest2)
##'
##' ## doesn't make much difference to the estimate
##' ## fit is OK in either case
##' 
##' @export
check_regression <- function(x, pars=NULL, n=NULL, comparison=1, outcome="costs", plot=TRUE){
  if (inherits(x, "evppi")) {
    if (is.null(pars)) pars <- x$pars[1]
    if (!(pars %in% x$pars)) stop(sprintf("parameter `%s` not found", pars))
    method <- attr(x, "methods")[match(pars, unique(x$pars))]
  }
  else if (inherits(x, "evsi")){
    if (is.null(n)) n <- x$n[1]
    else if (!(n %in% x$n)) stop(sprintf("sample size `%s` not found", n))
    pars <- as.character(n)
    method <- attr(x,"method")
  }
  else stop("`x` should be an object returned by evppi() or evsi()")
  if (method %in% npreg_methods){
    cea <- (attr(x, "outputs") == "cea")
    mods <- attr(x, "models")
    if (is.null(mods)) 
      stop("evppi() or evsi() should be run with `check=TRUE` to enable regression checks")
    if (inherits(x, "evppi")){
    } else if (inherits(x, "evsi")) {
    }
    ncomp <- if (cea) length(mods[[1]][[1]]) else length(mods[[1]])
    if (!(comparison %in% 1:ncomp)) stop(sprintf("`comparison` should be a positive integer <= %s", ncomp))
    if (cea){
      if (!(outcome %in% c("costs","effects"))) stop("`outcome` should be \"costs\" or \"effects\"")
      outcome <- if (outcome=="costs") "c" else "e"
      mod <- mods[[pars]][[outcome]][[comparison]]
    } else { 
      mod <- mods[[pars]][[comparison]]
    }
    if (plot) {
        check_plot_default(mod)
    }
  } else {
    message("`check_reg` is only applicable when method=\"gam\", \"earth\", \"gp\" or \"inla\"")
  }
  check_stats_fn <- get(sprintf("check_stats_%s", method))
  check_stats_fn(mod)
}

## Should work for any model for which fitted() and residuals() works.

check_plot_default <- function(mod){
  fit <- as.numeric(fitted(mod))
  res <- as.numeric(residuals(mod))
  if (!is.numeric(fit))
    warning("fitted() does not work on regression model object, so can't produce diagnostic plots")
  if (!is.numeric(res))
    warning("residuals() does not work on regression model object, so can't produce diagnostic plots")
  dat <- data.frame(fit = fit, res = res)
  bw <- 2 * IQR(dat$res) / length(dat$res)^(1/3)
  p1 <- ggplot2::ggplot(dat, aes(x=res)) +
    ggplot2::geom_histogram(binwidth=bw) +
    xlab("Residuals") + ylab("Frequency")
  p2 <- ggplot2::ggplot(dat, aes(x=fit, y=res)) +
    ggplot2::geom_point() +
    xlab("Fitted values") + ylab("Residuals")
  gridExtra::grid.arrange(p1, p2, nrow=1)
}
