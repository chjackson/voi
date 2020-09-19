##' Calculate the expected value of partial perfect information from a decision model
##'
##' Calculate the expected value of partial perfect information from a decision model
##'
##' @param outputs This could take one of two forms
##'
##' "net benefit" form: a matrix or data frame of samples from the uncertainty distribution of the expected net benefit.  The number of rows should equal the number of samples, and the number of columns should equal the number of decision options.
##'
##' "cost-effectiveness analysis" form: a list with the following named components:
##'
##' \code{"c"}: a matrix or data frame of samples from the distribution of costs.  There should be one column for each decision option.
##' 
##' \code{"e"}: a matrix or data frame of samples from the distribution of effects, likewise.
##'
##' \code{"k"}: a vector of willingness-to-pay values.
##'
##' Objects of class \code{"bcea"}, as created by the \pkg{BCEA} package, are in this "cost-effectiveness analysis" format, therefore they may be supplied as the \code{outputs} argument.
##'
##' If \code{outputs} is a matrix or data frame it is assumed to be of "net benefit" form.  Otherwise if it is a list, it is assumed to be of "cost effectiveness analysis" form.
##'
##' @param inputs Matrix or data frame of samples from the uncertainty distribution of the input parameters of the decision model.   The number of columns should equal the number of parameters, and the columns should be named.    This should have the same number of rows as there are samples in \code{outputs}, and each row of the samples in \code{outputs} should give the model output evaluated at the corresponding parameters.
##'
##' @param pars A character vector giving the parameters of interest, for which
##'   a single EVPPI calculation is required.  If the vector has multiple
##'   element, then the joint expected value of perfect information on all these
##'   parameters together is calculated.
##'
##'   Alternatively, \code{pars} may be a list.  Multiple EVPPI calculations are
##'   then performed, one for each component of \code{pars} defined in the above
##'   vector form.
##'
##'   \code{pars} must be specified if \code{inputs} is a matrix or data frame. This
##'   should then correspond to particular columns of \code{inputs}.    If
##'   \code{inputs} is a vector, this is assumed to define the single parameter
##'   of interest, then \code{pars} is not required.
##'   
##' @param method Character string indicating the calculation method.   The only methods currently implemented are based on nonparametric regression:
##'
##' \code{"gam"} for a generalized additive model implemented in the gam() function from the mgcv() package.
##' 
##' \code{"gp"} for a Gaussian process regression, as described by Strong et al. (2014) and implemented in the SAVI package (\url{http://savi.shef.ac.uk/SAVI/}).
##'
##' \code{"inla"} for an INLA/SPDE Gaussian process regression method, from Heath et al. (2016). 
##'
##' \code{"earth"} for a multivariate adaptive regression spline with the \pkg{earth} package (Milborrow, 2019). 
##'
##'
##' @param nsim Number of simulations from the model to use for calculating EVPPI.  The first \code{nsim} rows of the objects in \code{inputs} and \code{outputs} are used.
##'
##' @param verbose If \code{TRUE}, then print messages describing each step of the calculation.  Useful to see the progress of slow calculations.  Currently only supported by the \code{"inla"} EVPPI method. 
##'
##' @param tidy Only used if \code{outputs} is in cost-effectiveness form.  Then if
##'   \code{tidy=TRUE}, \code{\link{evppi}} returns a data frame containing a
##'   column \code{k} for the willingness-to-pay values and a column
##'   \code{evppi} for the EVPPIs.  Or if \code{tidy=FALSE} then the function
##'   returns a vector of the EVPPIs at each willingess-to-pay.
##'
##' @param ... Other arguments to control specific methods.
##'
##' For \code{method="gam"}:
##'
##' \code{gam_formula}: a character string giving the right hand side of the formula supplied to the \code{gam()} function. By default, this is a tensor product of all the parameters of interest, e.g. if \code{pars = c("pi","rho")}, then \code{gam_formula} defaults to \code{t(pi, rho, bs="cr")}.  The option \code{bs="cr"} indicates a cubic spline regression basis, which more computationally efficient than the default "thin plate" basis.  If there are four or more parameters of interest, then the additional argument \code{k=4} is supplied to \code{te()}, specifying a four-dimensional basis, which is currently the default in the SAVI package (\url{http://savi.shef.ac.uk/SAVI/}).
##'
##' For \code{method="gp"}:
##'
##' \code{gp_hyper_n}: number of samples to use to estimate the hyperparameters in the Gaussian process regression method.  By default, this is the minimum of the following three quantities: 30 times the number of parameters of interest, 250, and the number of simulations being used for calculating EVPPI.
##'
##' \code{maxSample} Maximum sample size to employ for \code{method="gp"}.  Only increase this from the default 5000 if your computer has sufficent memory to invert square matrices with this dimension.
##'
##' For \code{method="inla"}, as described in detail in Baio, Berardi and Heath (2017):
##'
##' \code{int.ord} (integer) maximum order of interaction terms to include in the regression predictor, e.g. if \code{int.ord=k} then all k-way interactions are used.  TODO handle effects and costs separately
##' 
#'  \code{cutoff} (default 0.3) controls the
#' density of the points inside the mesh in the spatial part of the mode. 
#' Acceptable values are typically in
#' the interval (0.1,0.5), with lower values implying more points (and thus
#' better approximation and greatercomputational time).
#'
#' \code{convex.inner} (default = -0.4) and \code{convex.outer} (default =
#' -0.7) control the boundaries for the mesh. 
#' These should be negative values and can be decreased (say to -0.7 and
#' -1, respectively) to increase the distance between the points and the outer
#' boundary, which also increases precision and computational time.
#'
#' \code{robust}. if \code{TRUE} then INLA will
#' use a t prior distribution for the coefficients of the linear predictor, rather than the default normal. 
#' 
#' \code{h.value} (default=0.00005) controls the accuracy of the INLA grid-search for the
#' estimation of the hyperparameters. 
#'  Lower values imply a more refined search
#' (and hence better accuracy), at the expense of computational speed.
#'
#' \code{plot_inla_mesh} (default \code{FALSE}) Produce a plot of the mesh. 
#'
#' \code{max.edge}  Largest allowed triangle edge length when constructing the mesh, passed to \code{\link[INLA]{inla.mesh.2d}}. 
#'
#'
#' @return If a single EVPPI calculation is being performed, when \code{pars} or
#'   \code{inputs} is a vector, and \code{outputs} is of "net benefit" form,
#'   then \code{\link{evppi}} simply returns the numeric value of the EVPPI.
#'   
#'   If \code{outputs} is of "cost-effectiveness analysis" form so that there is
#'   one EVPPI per willingness-to-pay value, then a data frame is returned with
#'   one column \code{k} identifying the willingness-to-pay, and another column
#'   \code{evppi} with the EVPPI.
#'   
#' If \code{pars} is a list, so that multiple EVPPI calculations are performed
#' with different parameters, then \code{\link{evppi}} also returns a data frame, with
#' another column \code{pars} identifying the parameter or parameters.
#' 
##' @references
##'
##' Strong, M., Oakley, J. E., & Brennan, A. (2014). Estimating multiparameter partial expected value of perfect information from a probabilistic sensitivity analysis sample: a nonparametric regression approach. Medical Decision Making, 34(3), 311-326.
##'
##' Heath, A., Manolopoulou, I., & Baio, G. (2016). Estimating the expected value of partial perfect information in health economic evaluations using integrated nested Laplace approximation. Statistics in medicine, 35(23), 4264-4280.
##'
##' Baio, G., Berardi, A., & Heath, A. (2017). Bayesian cost-effectiveness analysis with the R package BCEA. New York: Springer.
##'
##' Milborrow, S. (2019) earth: Multivariate Adaptive Regression Splines. R package version 5.1.2. Derived from mda:mars by Trevor Hastie and Rob Tibshirani. Uses Alan Miller's Fortran utilities with Thomas Lumley's leaps wrapper. https://CRAN.R-project.org/package=earth. 
##'
##' 
##' @export
evppi <- function(outputs,
                  inputs,
                  pars=NULL,
                  method=NULL,
                  nsim=NULL,
                  verbose=TRUE,
                  tidy=TRUE,
                  ...)
{
    inputs <- check_inputs(inputs, iname=deparse(substitute(inputs)))
    output_type <- check_outputs(outputs, inputs)
    if (is.list(pars))
        return(evppi_list(outputs, inputs, output_type, pars, method, nsim, verbose, ...))
    pars <- check_pars(pars, inputs)
    opts <- list(...)
    if (is.null(method))
        method <- default_evppi_method(pars)
    
    if (is.null(nsim)) nsim <- nrow(inputs)
    outputs <- subset_outputs(outputs, output_type, nsim)
    inputs <- inputs[1:nsim,,drop=FALSE]
    
    if (method %in% npreg_methods) {
        res <- evppi_npreg(outputs=outputs, inputs=inputs, output_type=output_type,
                    pars=pars, method=method, verbose=verbose, ...)
    } else stop("Other methods not implemented yet")
    if (output_type=="cea") 
        res <- if (tidy) data.frame(k=outputs$k, evppi=res) else res
    res
}

evppi_list <- function(outputs, inputs, output_type, pars, method, nsim, verbose, ...)
{
    npars <- length(pars)
    nouts <- if (output_type=="cea") length(outputs$k) else 1
    res <- data.frame(pars=character(npars*nouts)) 
    if (output_type=="cea") res$k <- rep(outputs$k, each=npars)
    res$evppi <- numeric(npars*nouts)
    indmat <- matrix(seq_len(npars*nouts), nrow=npars, ncol=nouts)
    for (i in seq_len(npars)){
        res$evppi[indmat[i,]] <- evppi(outputs, inputs, pars[[i]], method, nsim, verbose, tidy=FALSE, ...)
    }
    res$pars <- rep(sapply(pars, function(x)paste(x,collapse=",")), nouts)
    res
}

subset_outputs <- function(outputs, output_type, nsim){
    if (output_type == "nb") {
        outputs <- outputs[1:nsim,,drop=FALSE]
    } else if (output_type == "cea") {
        outputs$c <- outputs$c[1:nsim,,drop=FALSE]
        outputs$e <- outputs$e[1:nsim,,drop=FALSE]
    }
    outputs
}

npreg_methods <- c("gam", "gp", "inla", "earth")

evppi_npreg <- function(outputs, inputs, output_type, pars, method=NULL, verbose, ...){
    if (output_type == "nb")
        evppi_npreg_nb(nb=outputs, inputs=inputs, pars=pars, method=method,
                       verbose=verbose, ...)
    else if (output_type == "cea")
        evppi_npreg_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                        inputs=inputs, pars=pars, method=method,
                        verbose=verbose, ...)
}

evppi_npreg_nb <- function(nb, inputs, pars, method, verbose, ...){
    if (verbose) message("Fitting nonparametric regression") 
    fit <- fitted_npreg(nb, inputs=inputs, pars=pars, method=method,
                        verbose=verbose, ...)
    calc_evppi(fit)
}

evppi_npreg_cea <- function(costs, effects, wtp, inputs, pars, method, verbose, ...){
    nwtp <- length(wtp)
    res <- numeric(nwtp)
    if (verbose) message("Fitting nonparametric regression for costs") 
    cfit <- fitted_npreg(costs, inputs=inputs, pars=pars, method=method, verbose=verbose, ...)
    if (verbose) message("Fitting nonparametric regression for effects") 
    efit <- fitted_npreg(effects, inputs=inputs, pars=pars, method=method, verbose=verbose, ...)
    for (i in 1:nwtp){
        inbfit <- efit*wtp[i] - cfit
        res[i] <- calc_evppi(inbfit)
    }
    res
}

fitted_npreg <- function(nb, inputs, pars, method, verbose, ...){
    nopt <- ncol(nb)
    nsam <- nrow(nb)
    ## Transforming to incremental net benefit allows us to do one fewer regression
    ## Assume it doesn't matters which option is the baseline for this purpose
    inb <- nb[,1] - nb[, -1, drop=FALSE] #- nb[,1]
    fitted <- matrix(0, nrow=nsam, ncol=nopt)
    for (i in 1:(nopt-1)){
        if (verbose) message(sprintf("Decision option %s",i+1)) 
        fitted[,i+1] <- fitted_npreg_call(inb[,i], inputs, pars, method, verbose=verbose, ...) 
    }
    fitted
}

fitted_npreg_call <- function(y, inputs, pars, method, verbose, ...){
    if (method=="gam") {
        fitted <- fitted_gam(y, inputs, pars, ...)
    }
    if (method=="gp") {
        fitted <- fitted_gp(y, inputs, pars, ...)
    } 
    if (method=="inla") {
        fitted <- fitted_inla(y, inputs, pars, verbose, ...)
    } 
    if (method=="earth") {
        fitted <- fitted_earth(y, inputs, pars, ...)
    } 
    fitted
}

calc_evppi <- function(fit) {
    ## NAs are removed in BCEA here. Shouldn't we make users investigate them and remove by hand if they know what they are doing?  At least warn users if there are NAs in their samples
    mean(apply(fit, 1, max)) - max(colMeans(fit))
}

default_evppi_method <- function(pars){
    if (length(pars) <= 4) "gam" else "inla"
}

check_inputs <- function(inputs, iname=NULL){
    if (is.vector(inputs) && is.numeric(inputs)) {
        inputs <- data.frame(input = inputs)
        names(inputs) <- iname
    }
    if (!is.matrix(inputs) && !is.data.frame(inputs)){ 
        stop("`inputs` should be a numeric vector, matrix or data frame")
    }
    inputs 
}

check_outputs_matrix <- function(outputs, inputs, name){
    if (ncol(outputs) < 2)
        stop(sprintf("`%s` should have two or more columns", name)) # or else voi always zero
    if (nrow(outputs) != nrow(inputs))
        stop(sprintf("Number of rows of `%s` (%s) should equal the number of rows of `inputs` (%s)",
                     name, nrow(outputs), nrow(inputs)))
}

check_outputs <- function(outputs, inputs=NULL){
    if (is.matrix(outputs) || is.data.frame(outputs)){
        output_type <- "nb"
        if (!is.null(inputs)) # check not required for EVPI 
            check_outputs_matrix(outputs, inputs, "outputs")
    }
    else if (is.list(outputs)){
        output_type <- "cea"
        required_names <- c("c","e","k")
        for (i in required_names){
            if (!(i %in% names(outputs)))
                stop(sprintf("component named `(%s)` not found in `outputs` list", i))
        }
        if (!is.null(inputs)){
            check_outputs_matrix(outputs$c, inputs, "outputs$c")
            check_outputs_matrix(outputs$e, inputs, "outputs$e")
        }
        ## TODO Also check wtp
    }
    else stop("`outputs` should be a matrix, data frame or list, see help(evppi)")
    output_type
}


check_pars <- function(pars, inputs){
    if (is.null(pars)){
        if (ncol(inputs)==1)
            pars <- colnames(inputs)
        else stop("`pars` should be specified if there are two or more parameters in `inputs`")
    }
    if (!is.character(pars))
        stop("`pars` should be a character vector")
    badpars <- pars[!(pars %in% colnames(inputs))]
    if (length(badpars)>0){
        stop(sprintf("parameters of interest `%s` not found in columns of `inputs`",
                     paste(badpars,collapse=",")))
    }
    pars
}

