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
##' Objects of class \code{"bcea"}, as created by the \pkg{BCEA} package, may be supplied as the \code{outputs} argument, as they are lists with the required components. 
##'
##' If \code{outputs} is a matrix or data frame it is assumed to be of "net benefit" form.  Otherwise if it is a list, it is assumed to be of "cost effectiveness analysis" form.
##'
##' The exact names and formats of all these arguments is up for discussion! 
##'
##' @param inputs Matrix or data frame of samples from the uncertainty distribution of the input parameters of the decision model.   The number of columns should equal the number of parameters, and the columns should be named.    This should have the same number of rows as there are samples in \code{outputs}, and each row of the samples in \code{outputs} should give the model output evaluated at the corresponding parameters.
##'
##' @param poi A character vector giving the parameters of interest, for which the EVPPI is required.   This should correspond to particular columns of \code{inputs}.  [ I thought about allowing numeric indices for this, but perhaps better to encourage good practice by mandating that the parameters are named ] 
##'
##' @param method Character string indicating the calculation method.
##'
##' "gam" for a generalized additive model implemented in the gam() function from the mgcv() package.
##' "gp" for a Gaussian process regression TODO MORE 
##'
##' "inla" for an INLA/SPDE Gaussian process regression  TODO MORE 
##'
##' "earth" for a multivariate adaptive regression spline with the \pkg{earth} package  TODO MORE 
##'
##' TODO earth, single par methods
##'
##' @param nsim Number of simulations from the model to use for calculating EVPPI.  The first \code{nsim} rows of the objects in \code{inputs} and \code{outputs} are used. 
##'
##' @param ... Other arguments required by specific methods
##'
##' \code{gam_formula}: a character string giving the right hand side of the formula supplied to the gam() function, when \code{method="gam"}. By default, this is a tensor product of all the parameters of interest, e.g. if \code{poi = c("pi","rho")}, then \code{gam_formula} defaults to \code{t(pi, rho, bs="cr")}.  The option \code{bs="cr"} indicates a cubic spline regression basis, which more computationally efficient than the default "thin plate" basis.  If there are four parameters of interest, then the additional argument \code{k=4} is supplied to \code{te()}, specifying a four-dimensional basis.   [ This is the default in SAVI ] 
##' 
##' @export
evppi <- function(outputs,
                  inputs,
                  poi,
                  method=NULL,
                  nsim=NULL,
                  ...)
{
    check_inputs(inputs)
    output_type <- check_outputs(outputs, inputs)
    check_poi(poi, inputs)
    opts <- list(...)
    if (is.null(method))
        method <- default_evppi_method(poi)
    if (is.null(nsim)) nsim <- nrow(inputs)
    if (method %in% npreg_methods) { 
        if (output_type == "nb")
            evppi_npreg_nb(nb=outputs, inputs=inputs, poi=poi, method=method, nsim=nsim, ...)
        else if (output_type == "cea")
            evppi_npreg_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                            inputs=inputs, poi=poi, method=method, nsim=nsim, ...)
    } else stop("Other methods not implemented yet")
}

npreg_methods <- c("gam", "gp", "inla", "earth")

evppi_npreg_nb <- function(nb, inputs, poi, method, nsim, ...){
    ## TODO move the subsetting outside when do methods other than np regression
    nb <- nb[1:nsim,,drop=FALSE]
    inputs <- inputs[1:nsim,,drop=FALSE]
    fit <- fitted_npreg(nb, inputs=inputs, poi=poi, method=method, ...)
    calc_evppi(fit)
}

evppi_npreg_cea <- function(costs, effects, wtp, inputs, poi, method, nsim, ...){
    costs <- costs[1:nsim,,drop=FALSE]
    effects <- effects[1:nsim,,drop=FALSE]
    inputs <- inputs[1:nsim,,drop=FALSE]
    nwtp <- length(wtp)
    res <- numeric(nwtp)
    for (i in 1:nwtp){
        cfit <- fitted_npreg(costs, inputs=inputs, poi=poi, method=method, ...)
        efit <- fitted_npreg(effects, inputs=inputs, poi=poi, method=method, ...)
        inbfit <- efit*wtp[i] - cfit
        res[i] <- calc_evppi(inbfit)
    }
    res
}

fitted_npreg <- function(nb, inputs, poi, method, ...){
    nopt <- ncol(nb)
    nsam <- nrow(nb)
    ## Transforming to incremental net benefit allows us to do one fewer regression
    ## Assume it doesn't matters which option is the baseline for this purpose
    inb <- nb[,1] - nb[, -1, drop=FALSE] #- nb[,1]
    fitted <- matrix(0, nrow=nsam, ncol=nopt)
    for (i in 1:(nopt-1)){
        fitted[,i+1] <- fitted_npreg_call(inb[,i], inputs, poi, method, ...) 
    }
    fitted
}

fitted_npreg_call <- function(y, inputs, poi, method, ...){
## could make this neater with do.call? 
    if (method=="gam") {
        fitted <- fitted_gam(y, inputs, poi, ...)
    } else 
    if (method=="gp") {
        fitted <- fitted_gp(y, inputs, poi, ...)
    } 
    if (method=="inla") {
        fitted <- fitted_inla(y, inputs, poi, ...)
    } 
    if (method=="earth") {
        fitted <- fitted_earth(y, inputs, poi, ...)
    } 
    fitted
}

calc_evppi <- function(fit) {
    ## NAs are removed in BCEA here. Shouldn't we make users investigate them and remove by hand if they know what they are doing?  At least warn users if there are NAs in their samples
    mean(apply(fit, 1, max)) - max(colMeans(fit))
}

default_evppi_method <- function(poi){
    if (length(poi) <= 4) "gam" else "INLA"
}

check_inputs <- function(inputs){
    if (!is.matrix(inputs) && !is.data.frame(inputs))
        stop("`inputs` should be a matrix or data frame")
}

check_outputs_matrix <- function(outputs, inputs, name){
    if (ncol(outputs) < 2)
        stop(sprintf("`%s` should have two or more columns", name)) # or else voi always zero
    if (nrow(outputs) != nrow(inputs))
        stop(sprintf("Number of rows of `%s` (%s) should equal the number of rows of `inputs` (%s)",
                     name, nrow(outputs), nrow(inputs)))
}

check_outputs <- function(outputs, inputs){
    if (is.matrix(outputs) || is.data.frame(outputs)){
        output_type <- "nb"
        check_outputs_matrix(outputs, inputs, "outputs")
    }
    else if (is.list(outputs)){
        output_type <- "cea"
        required_names <- c("c","e","k")
        for (i in required_names){
            if (!(i %in% names(outputs)))
                stop(sprintf("component named `(%s)` not found in `outputs` list", i))
        }
        check_outputs_matrix(outputs$c, inputs, "outputs$c")
        check_outputs_matrix(outputs$e, inputs, "outputs$e")
        ## TODO Also check wtp
    }
    else stop("`outputs` should be a matrix, data frame or list, see help(evppi)")
    output_type
}


check_poi <- function(poi, inputs){
    if (!is.character(poi))
        stop("`poi` should be a character vector")
    badpoi <- poi[!(poi %in% colnames(inputs))]
    if (length(badpoi)>0){
        stop(sprintf("parameters of interest `%s` not found in columns of `inputs`",
                     paste(badpoi,collapse=",")))
    }
}

