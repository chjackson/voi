
npreg_methods <- c("gam", "gp", "inla", "earth")

evppi_npreg <- function(outputs, ...){
    UseMethod("evppi_npreg", outputs)
}

evppi_npreg.nb <- function(outputs, inputs, pars, method, se, B, verbose, ...){
    if (verbose) message("Fitting nonparametric regression") 
    fit <- fitted_npreg(outputs, inputs=inputs, pars=pars, method=method, se=se, B=B, 
                        verbose=verbose, ...)
    res <- data.frame(evppi = calc_evppi(fit))
    if (se){
        if (verbose) message("Calculating replicate EVPPIs from simulation")
        evppi_rep <- numeric(B)
        for (i in 1:B){
            ## This bit could be faster.  Is there a vectorised way?  Naive apply doesn't help
            evppi_rep[i] <- calc_evppi(attr(fit, "rep")[i,,])   
        }
        res$se <-  sd(evppi_rep)
    }
    attr(res, "models") <- attr(fit, "models")
    res
}

evppi_npreg.cea <- function(outputs, inputs, pars, method, se, B, verbose, ...){
    wtp <- outputs$k
    nwtp <- length(wtp)
    res <- resse <- numeric(nwtp)
    if (verbose) message("Fitting nonparametric regression for costs") 
    cfit <- fitted_npreg(outputs$c, inputs=inputs, pars=pars, method=method, se=se, B=B, verbose=verbose, ...)
    if (verbose) message("Fitting nonparametric regression for effects") 
    efit <- fitted_npreg(outputs$e, inputs=inputs, pars=pars, method=method, se=se, B=B, verbose=verbose, ...)
    for (i in 1:nwtp){
        evppi_rep <- numeric(B)
        inbfit <- efit*wtp[i] - cfit
        if (se){
            inbrep <- attr(efit, "rep")*wtp[i] - attr(cfit, "rep")
            if (verbose) message("Calculating replicate EVPPIs from simulation")
            for (j in 1:B){
                evppi_rep[j] <- calc_evppi(inbrep[j,,])   
            }
            resse[i] <- sd(evppi_rep)
        }
        res[i] <- calc_evppi(inbfit)
    }
    res <- data.frame(k=wtp, evppi=res)
    attr(res, "models") <- list(c=attr(cfit, "models"), e=attr(efit, "models"))
    if (se) res$se <- resse
    res
}

fitted_npreg <- function(nb, inputs, pars, method, se=FALSE, B, verbose, ...){
    nopt <- ncol(nb)
    nsim <- nrow(nb)
    ## Transforming to incremental net benefit allows us to do one fewer regression
    ## Assume it doesn't matters which option is the baseline for this purpose
    inb <- nb[, -1, drop=FALSE] - nb[,1]
    fitted <- matrix(0, nrow=nsim, ncol=nopt)
    if (se) 
        fitted_rep <- array(0, dim=c(B, nsim, nopt))
    models <- vector(nopt-1, mode="list")
    for (i in 1:(nopt-1)){
        if (verbose) message(sprintf("Decision option %s",i+1)) 
        fit <- fitted_npreg_fn(method)(y=inb[,i], inputs=inputs, pars=pars, verbose=verbose, ...) 
        fitted[,i+1] <- as.vector(fit)
        if (se){
            fitted_rep[,,i+1] <- fitted_npreg_rep_call(method, attr(fit,"model"), B, verbose)
        }
        models[[i]] <- attr(fit, "model")
    }
    if (se) attr(fitted, "rep") <- fitted_rep
    attr(fitted, "models") <- models
    fitted
}

fitted_npreg_fn <- function(method){
    switch(method, 
           gam = fitted_gam,
           gp = fitted_gp,
           inla = fitted_inla,
           earth = fitted_earth)
}

fitted_npreg_rep_call <- function(method, model, B, verbose=FALSE){
    if (verbose) message("Simulating parameters to calculate standard errors")
    if (method=="gam") {
        frep <- fitted_rep_gam(model, B)
    } else stop("Standard errors only currently available for GAM method")
    frep
}

calc_evppi <- function(fit) {
    ## NAs are removed in BCEA here. Shouldn't we make users investigate them and remove by hand if they know what they are doing?  At least warn users if there are NAs in their samples
    mean(apply(fit, 1, max)) - max(colMeans(fit))
}

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
##' Models produced with `method="gam"` are checked using \code{\link{gam.check}}
##'
##' Models produced `method="earth"` are checked using \code{\link{plot.earth}}
##' 
##' Models produced with `method="gp"` or `method="inla"` are checked using a histogram
##' of the residuals, and plots of the residuals and response against the fitted values. 
##'
##' @param x Output from \code{\link{evppi}} or \code{\link{evsi}}. The argument \code{check=TRUE}
##' must have been used when calling \code{evppi} or \code{evsi}, to allow the regression model
##' objects from \code{gam} or \code{earth} to be preserved.  (This is not done by
##' default, since these objects can be large.).   \code{attr(x, "models")} cont
##' 
##' @param pars Parameter (or parameter group) whose EVPPI calculation is to be checked.
##' This should be in the \code{pars} component of the object returned by \code{\link{evppi}}.
##' Only relevant if \code{x} is the result of an \code{\link{evppi}} calculation.
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
##' @return Where possible, an appropriate statistic is returned that allows the regression
##' model to be compared with other regression models implemented using the same \code{method} 
##' but with different assumptions.   For \code{method=="gam"},
##' this is Akaike's information criterion (AIC).    
##' For \code{method=="earth"}, this is the generalised cross-validation statistic
##' \code{gcv}.    Currently not implemented for other methods. 
##'
##' @examples # TODO and refer in vignette
##' 
##' @export
check_regression <- function(x, pars=NULL, n=NULL, comparison=1, outcome="costs"){
  ## TODO check for evsi vs evppi here.
  ## then treat n in evsi like pars in evppi 
  if (inherits(x, "evppi")) {
    if (is.null(pars)) pars <- x$pars[1]
    if (!(pars %in% x$pars)) stop(sprintf("parameter `%s` not found", pars))
    method <- attr(x, "methods")[match(pars, x$pars)]
  }
  else if (inherits(x, "evsi")){
    if (is.null(n)) pars <- as.character(x$n[1])
    if (!(n %in% x$n)) stop(sprintf("sample size `%s` not found", n))
    method <- attr(x,"method")
  }
  else stop("`x` should be an object returned by evppi() or evsi()")
  if (method %in% npreg_methods){
    cea <- (attr(x, "outputs") == "cea")
    mods <- attr(x, "models")
    if (is.null(mods)) 
      stop("evppi() or evsi() should be run with `check=TRUE` to enable regression checks")
    check_fn <- if (method=="gam") mgcv::gam.check else if (method %in% c("gp","inla")) gp.check else plot
    if (inherits(x, "evppi")){
    } else if (inherits(x, "evsi")) {
    }
    ncomp <- if (cea) length(mods[[1]][[1]]) else length(mods[[1]])
    if (!(comparison %in% 1:ncomp)) stop(sprintf("`comparison` should be a positive integer <= %s", ncomp))
    if (cea){
      if (!(outcome %in% c("costs","effects"))) stop("`outcome` should be \"costs\" or \"effects\"")
      mod <- mods[[pars]][[outcome]][[comparison]]
    } else { 
      mod <- mods[[pars]][[comparison]]
    }
    check_fn(mod)
  } else {
    message("`check_reg` is only applicable when method=\"gam\", \"earth\", \"gp\" or \"inla\"")
  }
  check_regression_stats(mod, method)
}

check_regression_stats <- function(mod, method){
  if (method=="gam") 
    list(AIC = stats::AIC(mod)) 
  else if (method=="earth")
    list(gcv = mod$gcv)
  else invisible()
}