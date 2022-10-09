##' Calculate the expected value of partial perfect information from a decision model
##'
##' Calculate the expected value of partial perfect information from a decision model
##'
##' @param outputs This could take one of two forms
##'
##'   "net benefit" form: a matrix or data frame of samples from the uncertainty
##'   distribution of the expected net benefit.  The number of rows should equal
##'   the number of samples, and the number of columns should equal the number
##'   of decision options.
##'
##'   "cost-effectiveness analysis" form: a list with the following named
##'   components:
##'
##'   \code{"c"}: a matrix or data frame of samples from the distribution of
##'   costs.  There should be one column for each decision option.
##'
##'   \code{"e"}: a matrix or data frame of samples from the distribution of
##'   effects, likewise.
##'
##'   \code{"k"}: a vector of willingness-to-pay values.
##'
##' Objects of class \code{"bcea"}, as created by the \pkg{BCEA} package, are in
##' this "cost-effectiveness analysis" format, therefore they may be supplied as
##' the \code{outputs} argument.
##'
##' If \code{outputs} is a matrix or data frame, it is assumed to be of "net
##' benefit" form.  Otherwise if it is a list, it is assumed to be of "cost
##' effectiveness analysis" form.
##'
##' @param inputs Matrix or data frame of samples from the uncertainty
##'   distribution of the input parameters of the decision model.   The number
##'   of columns should equal the number of parameters, and the columns should
##'   be named.    This should have the same number of rows as there are samples
##'   in \code{outputs}, and each row of the samples in \code{outputs} should
##'   give the model output evaluated at the corresponding parameters.
##'
##' @param pars A character vector giving the parameters of interest, for which
##'   a single EVPPI calculation is required.  If the vector has multiple
##'   elements, then the joint expected value of perfect information on all these
##'   parameters together is calculated.
##'
##'   Alternatively, \code{pars} may be a list.  Multiple EVPPI calculations are
##'   then performed, one for each component of \code{pars} defined in the above
##'   vector form.
##'
##'   \code{pars} must be specified if \code{inputs} is a matrix or data frame.
##'   This should then correspond to particular columns of \code{inputs}.    If
##'   \code{inputs} is a vector, this is assumed to define the single parameter
##'   of interest, and then \code{pars} is not required.
##'
##' @param method Character string indicating the calculation method.  If one
##' string is supplied, this is used for all calculations.  A vector of different strings
##' can be supplied if a different method is desired for different list components
##' of \code{pars}. 
##'
##' The default methods are based on nonparametric regression:
##'
##' \code{"gam"} for a generalized additive model implemented in the \code{\link{gam}}
##' function from the \pkg{mgcv} package.  This is the default method for
##' calculating the EVPPI of 4 or fewer parameters.
##'
##' \code{"gp"} for a Gaussian process regression, as described by Strong et al.
##' (2014) and implemented in the \pkg{SAVI} package
##' (\url{http://savi.shef.ac.uk/SAVI/}).  This is the default method for calculating the EVPPI
##' of more than 4 parameters.
##'
##' \code{"inla"} for an INLA/SPDE Gaussian process regression method, from
##' Heath et al. (2016).
##'
##' \code{"bart"} for Bayesian additive regression trees, using the \pkg{dbarts} package.
##' Particularly suited for joint EVPPI of many parameters. 
##' 
##' \code{"earth"} for a multivariate adaptive regression spline with the
##' \pkg{earth} package (Milborrow, 2019).
##' 
##' \code{"so"} for the method of Strong and Oakley (2013).  Only supported
##' for single parameter EVPPI. 
##' 
##' \code{"sal"} for the method of Sadatsafavi et al. (2013).  Only supported 
##' for single parameter EVPPI.
##' 
##' @param se If this is \code{TRUE}, calculate a standard error for the EVPPI
##'  if possible.  Currently only supported for \code{method="gam"}, and
##' \code{method="bart"}.  (In the latter method it is more correctly called
##' a posterior standard deviation).
##'
##' @param B Number of parameter replicates for calculating the standard error.
##' Only applicable to \code{method="gam"}.  For \code{method="bart"} the
##' analogous quantity is the number of MCMC samples, which is controlled by
##' the \code{ndpost} argument to \code{\link[dbarts]{bart}}, which can be
##' passed as an argument to \code{\link{evppi}}. 
##'   
##' @param nsim Number of simulations from the model to use for calculating
##'   EVPPI.  The first \code{nsim} rows of the objects in \code{inputs} and
##'   \code{outputs} are used.
##'
##' @param verbose If \code{TRUE}, then messages are printed describing each step of
##'   the calculation, if the method supplies these.  Useful to see the progress
##'   of slow calculations.
##'
##' @param check If \code{TRUE}, then extra information about the estimation
##' is saved inside the object that this function returns.  This currently
##' only applies to the regression-based methods \code{"gam"} and \code{"earth"}
##' where the fitted regression model objects are saved.  This allows use
##' of the \code{\link{check_regression}} function, which produces some
##' diagnostic checks of the regression models.
##'
##' @param ... Other arguments to control specific methods.
##'
##'   For \code{method="gam"}:
##'   
##' * \code{gam_formula}: a character string giving the right hand side of the
##' formula supplied to the \code{gam()} function. By default, this is a tensor
##' product of all the parameters of interest, e.g. if \code{pars =
##' c("pi","rho")}, then \code{gam_formula} defaults to \code{t(pi, rho,
##' bs="cr")}.  The option \code{bs="cr"} indicates a cubic spline regression
##' basis, which is more computationally efficient than the default "thin plate"
##' basis.  If there are four or more parameters of interest, then the
##' additional argument \code{k=4} is supplied to \code{te()}, specifying a
##' four-dimensional basis, which is currently the default in the SAVI package
##' (\url{http://savi.shef.ac.uk/SAVI/}).
##' 
##' If there are spaces in the variable names in \code{inputs}, then these should
##' be converted to underscores before forming an explicit \code{gam_formula}.
##' 
##'
##' For \code{method="gp"}:
##'
##' * \code{gp_hyper_n}: number of samples to use to estimate the hyperparameters
##' in the Gaussian process regression method.  By default, this is the minimum
##' of the following three quantities: 30 times the number of parameters of
##' interest, 250, and the number of simulations being used for calculating
##' EVPPI.
##'
##' * \code{maxSample}: Maximum sample size to employ for \code{method="gp"}.  Only
##' increase this from the default 5000 if your computer has sufficent memory to
##' invert square matrices with this dimension.
##'
##' For \code{method="inla"}, as described in detail in Baio, Berardi and Heath
##' (2017):
##'
##' * \code{int.ord} (integer) maximum order of interaction terms to include in
##' the regression predictor, e.g. if \code{int.ord=k} then all k-way
##' interactions are used.  Currently this applies to both effects and costs.
##' 
#' * \code{cutoff} (default 0.3) controls the
#' density of the points inside the mesh in the spatial part of the mode. 
#' Acceptable values are typically in
#' the interval (0.1,0.5), with lower values implying more points (and thus
#' better approximation and greatercomputational time).
#'
#' * \code{convex.inner} (default = -0.4) and \code{convex.outer} (default = -0.7)
#' control the boundaries for the mesh. These should be negative values and can
#' be decreased (say to -0.7 and -1, respectively) to increase the distance
#' between the points and the outer boundary, which also increases precision and
#' computational time.
#'
#' * \code{robust}. if \code{TRUE} then INLA will use a t prior distribution for
#' the coefficients of the linear predictor, rather than the default normal.
#'
#' * \code{h.value} (default=0.00005) controls the accuracy of the INLA
#' grid-search for the estimation of the hyperparameters. Lower values imply a
#' more refined search (and hence better accuracy), at the expense of
#' computational speed.
#'
#' * \code{plot_inla_mesh} (default \code{FALSE}) Produce a plot of the mesh.
#'
#' * \code{max.edge}  Largest allowed triangle edge length when constructing the
#' mesh, passed to \code{\link[INLA]{inla.mesh.2d}}.
#' 
#' * \code{pfc_struc} Variance structure to pass to \code{pfc} in the \pkg{ldr}
#' package for principal fitted components. The default \code{"AIC"} selects the
#' one that fits best given two basis terms.  Change this to, e.g. \code{"iso"},
#' \code{"aniso"} or \code{"unstr"} if an "Error in eigen..." is obtained.
#'
##' For \code{method="so"}:
##'
##' * \code{n.blocks} Number of blocks to split the sample into. Required.
##'
##' For \code{method="sal"}:
##'
##' * \code{n.seps} Number of separators (default 1). 
#'
#' @return A data frame with a column \code{pars} indicating the parameter(s)
#'   and a column \code{evppi} giving the corresponding EVPPI. 
#'
#'   If \code{outputs} is of "cost-effectiveness analysis" form so that there is
#'   one EVPPI per willingness-to-pay value, then a column \code{k} identifies the 
#'   willingness-to-pay.
#'   
#'   If standard errors are requested, then the standard errors are returned in 
#'   the column \code{se}.
#'   
##' @references
##'
##' Strong, M., Oakley, J. E., & Brennan, A. (2014). Estimating multiparameter
##' partial expected value of perfect information from a probabilistic
##' sensitivity analysis sample: a nonparametric regression approach. Medical
##' Decision Making, 34(3), 311-326.
##'
##' Heath, A., Manolopoulou, I., & Baio, G. (2016). Estimating the expected
##' value of partial perfect information in health economic evaluations using
##' integrated nested Laplace approximation. Statistics in medicine, 35(23),
##' 4264-4280.
##'
##' Baio, G., Berardi, A., & Heath, A. (2017). Bayesian cost-effectiveness
##' analysis with the R package BCEA. New York: Springer.
##'
##' Milborrow, S. (2019) earth: Multivariate Adaptive Regression Splines. R
##' package version 5.1.2. Derived from mda:mars by Trevor Hastie and Rob
##' Tibshirani. Uses Alan Miller's Fortran utilities with Thomas Lumley's leaps
##' wrapper. https://CRAN.R-project.org/package=earth.
##'
##' Strong, M., & Oakley, J. E. (2013). An efficient method for computing
##' single-parameter partial expected value of perfect information. Medical
##' Decision Making, 33(6), 755-766. Chicago
##'
##' Sadatsafavi, M., Bansback, N., Zafari, Z., Najafzadeh, M., & Marra, C.
##' (2013). Need for speed: an efficient algorithm for calculation of
##' single-parameter expected value of partial perfect information. Value in
##' Health, 16(2), 438-448.
##'
##' 
##' @export
evppi <- function(outputs,
                  inputs,
                  pars=NULL,
                  method=NULL,
                  se=FALSE,
                  B=1000,
                  nsim=NULL,
                  verbose=FALSE,
                  check=FALSE,
                  ...)
{
    inputs <- check_inputs(inputs, iname=deparse(substitute(inputs)))
    outputs <- check_outputs(outputs, inputs)

    if (!is.list(pars))
        pars <- list(pars)
    for (i in seq_along(pars)){
        pars[[i]] <- check_pars(pars[[i]], inputs)
        if (is.null(names(pars)) || identical(names(pars)[i], "") || is.na(names(pars)[i]))
            names(pars)[i] <- paste(pars[[i]], collapse=",")
    }
    npars <- length(pars)
    if (is.null(nsim)) nsim <- nrow(inputs)
    outputs <- subset_outputs(outputs, nsim)
    inputs <- inputs[1:nsim,,drop=FALSE]

    if (is.null(method))
        methods <- sapply(pars, default_evppi_method)
    if (length(method) > 0) methods <- rep(method, length.out=npars)
    eres <- vector(npars, mode="list") 
    for (i in seq_len(npars)){
        if (methods[i] %in% npreg_methods) {
            evppi_fn <- evppi_npreg
        } else if (methods[i]=="so") {
            evppi_fn <- evppi_so
        } else if (methods[i]=="sal") {
            evppi_fn <- evppi_sal
        } else stop("Other methods not implemented yet")
        eres[[i]] <- evppi_fn(outputs=outputs, inputs=inputs, pars=pars[[i]], 
                              method=methods[i], se=se, B=B, verbose=verbose,
                              ...)
    }
    res <- do.call("rbind", eres)
    nwtp <- if (inherits(outputs, "nb")) 1 else length(outputs$k)
    res <- cbind(pars=rep(names(pars), each = nwtp), res)
    if (check){
        attr(res, "models") <- lapply(eres, function(x)attr(x, "models"))
        names(attr(res, "models")) <- res$pars
    }
    attr(res, "methods") <- methods
    attr(res, "outputs") <- class(outputs)[1]
    class(res) <- c("evppi", attr(res,"class"))
    res
}

## could do fancier S3 stuff with implementing subset operator, but too much
## faff

subset_outputs <- function(outputs, ...){
    UseMethod("subset_outputs", outputs)
}

subset_outputs.nb <- function(outputs, nsim){
    outputs <- outputs[1:nsim,,drop=FALSE]
    class(outputs) <- c("nb", attr(outputs, "class"))
    outputs 
}

subset_outputs.cea <- function(outputs, nsim){
    outputs$c <- outputs$c[1:nsim,,drop=FALSE]
    outputs$e <- outputs$e[1:nsim,,drop=FALSE]
    class(outputs) <- c("cea", attr(outputs, "class"))
    outputs
}

default_evppi_method <- function(pars){
    if (length(pars) <= 4) "gam" else "gp"
}

check_inputs <- function(inputs, iname=NULL){
    if (is.vector(inputs) && is.numeric(inputs)) {
        inputs <- data.frame(input = inputs)
        names(inputs) <- gsub(" ", "", iname)
    }
    if (!is.matrix(inputs) && !is.data.frame(inputs)){ 
        stop("`inputs` should be a numeric vector, matrix or data frame")
    }
    as.data.frame(inputs) 
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
        class(outputs) <- c("nb", attr(outputs, "class"))
        if (!is.null(inputs)) # check not required for EVPI 
            check_outputs_matrix(outputs, inputs, "outputs")
    }
    else if (is.list(outputs)){
        class(outputs) <- c("cea", attr(outputs, "class"))
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
    outputs
}


## TODO which of these are needed for EVSI 

check_pars <- function(pars, inputs, evppi=TRUE){
    if (is.null(pars) && evppi){
        if (ncol(inputs)==1)
            pars <- colnames(inputs)
        else stop("`pars` should be specified if there are two or more parameters in `inputs`")
    }
    if (!is.null(pars) && !is.character(pars))
        stop("`pars` should be a character vector")
    badpars <- pars[!(pars %in% colnames(inputs))]
    if (length(badpars)>0){
        stop(sprintf("parameters of interest `%s` not found in columns of `inputs`",
                     paste(badpars,collapse=",")))
    }
    pars
}

clean_pars <- function(pars) {
    parsc <- gsub(" ", "_", pars)
    r_specials <- c("letters","month.abb","month.name","pi")
    for (i in seq_along(r_specials)){
        inds <- parsc == r_specials[i]
        if (any(inds)){
            for (j in which(inds)){
                stop(sprintf("Parameter name `%s` is also the name of a R internal constant. This should be changed to another name to allow the `gam` method for VoI calculation to be used", parsc[j]))
            }
        }
    }
    parsc
}

## TODO TEST AND INTEGRATE THIS
remove_constant_linear_cols <- function(inputs, verbose=TRUE){
  inputs <- as.matrix(inputs)
  p <- ncol(inputs)
  sets <- seq_len(p)
  constantParams <- (apply(inputs, 2, var) == 0)
  if (sum(constantParams) == p)
    stop("EVPI is zero as all parameters are constant")
  if (sum(constantParams) > 0)
    sets <- sets[-which(constantParams)] # remove constants
  inputs <- inputs[, sets, drop=FALSE] # now with constants removed

  rankifremoved <- sapply(1:NCOL(inputs), function (x) qr(inputs[, -x])$rank)
  while(length(unique(rankifremoved)) > 1) {
    linearCombs <- which(rankifremoved == max(rankifremoved))
    if (verbose) 
        print(paste("Linear dependence: removing column", colnames(inputs)[max(linearCombs)]))
    inputs <- inputs[, -max(linearCombs), drop=FALSE]
    sets <- sets[-max(linearCombs)]
    rankifremoved <- sapply(1:NCOL(inputs), function(x) qr(inputs[, -x])$rank)
  }  
  if(qr(inputs)$rank == rankifremoved[1]) {
    inputs <- inputs[, -1, drop=FALSE] # special case only lincomb left
    sets <- sets[-1]
    if (verbose) 
        print(paste("Linear dependence: removing column", colnames(inputs)[1]))
  }
  inputs[, sets, drop=FALSE]
}

