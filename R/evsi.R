##' Calculate the expected value of sample information from a decision-analytic model
##'
##' Calculate the expected value of sample information from a decision-analytic model
##'
##' @inheritParams evppi
##' 
##' @param datagen_fn Function to sample predicted data from a proposed future study.  This should have the following specification:
##'
##' 1. the function's first argument should be a data frame of parameter simulations, with one row per simulation and one column per parameter.  The parameters in this data frame must all be found in \code{inputs}.
##'
##' 2. the function should return a data frame.
##'
##' 3. the returned data frame should have number of rows equal to the number of parameter simulations in \code{inputs}.
##'
##' 4. if \code{inputs} is considered as a sample from the posterior, then \code{datagen_fn} returns a corresponding sample from the posterior predictive distribution, which includes two sources of uncertainty: (a) uncertainty about the parameters and (b) sampling variation in observed data given fixed parameter values.
##'
##' 5. the function can optionally have more than one argument. If so, these additional arguments should be given default values in the definition of \code{datagen_fn}.  These arguments might be used to define sample sizes for a proposed study.
##'
##' DEVELOPERS - EXAMPLES OF THIS CURRENTLY in \code{tests/tests_slow}
##'
##' @param n Sample size of future study - optional argument to datagen_fn - facilitates calculating EVSI for multiple sample sizes.  TODO how will this work for randomised trials with multiple arms? 
##'
##' @param method Character string indicating the calculation method.
##'
##' All the nonparametric regression methods supported for \code{\link{evppi}}, that is \code{"gam","gp","earth","inla"}, can also be used for EVSI calculation by regressing on a summary statistics of the predicted data (Strong et al 201?).   Defaults to \code{"gam"}.
##'
##' \code{"is"} for importance sampling (Menzies 2016)
##'
##' \code{"mm"} for moment matching (Heath et al 2018)
##'
##' Note that the  \code{"is"} and \code{"mm"} (and Jalal) methods are used in conjunction with nonparametric regression, thus the \code{gam_formula} argument can be supplied to \code{evsi} to specify this regression - see \code{\link{evppi}}. 
##' 
##' TODO Heath and Jalal methods
##'
##' @param likelihood Likelihood function, required (and only required) for the importance sampling method.  This should take arguments:
##'
##' the first: a data frame with columns defined by the number of outcomes in the data, and with names matching the names of the data frame returned by \code{datagen_fn}. 
##'
##' the second: a data frame of parameter values, whose names should all correspond to variables in \code{inputs}.
##'
##' The function should return a vector whose length matches the number of rows of the data frame given as the second argument. 
##'
##' POINT TO AN EXAMPLE WHICH WILL MAKE ALL THIS CLEARER.  EXAMPLE OF THIS CURRENTLY in \code{tests/tests_slow} but need a vignette eventually
##'
##' Note the definition of the likelihood should agree with the definition of \code{datagen_fn} to define a consistent sampling distribution for the data.   CLARIFY.  [ eventually we'll want some built-in common examples where people don't have to specify either datagen_fn or likelihood ]
##'
##' @param analysis_model Function which fits a Bayesian model to the generated data. TODO work out format, output, JAGS dependencies, etc.  Required for \code{method="mm"} (and Jalal method if n0 not given??) 
##'
##' @param model Function which evaluates the decision-analytic model, given parameter values.  TODO sort out format for nb, c, e output.   TODO should all function arguments have \code{_fn} suffix, or not? be consistent here
##'
##' @param Q Number of quantiles to use in \code{method="mm"}. 
##'
##' @param poi Parameters of interest, that is, those which are informed by the data in the future study.  Required (and only required) for the methods which involve an intermediate EVPPI calculation, that is the \code{"is"} and \code{"mm"} TODO OTHER methods.
##'
##' This should bee a character vector naming particular columns of \code{inputs}.  It should consist of the variables used in the definition of \code{datagen_fn} (and \code{likelihood} if used TODO ALSO in \code{analysis_model} and \code{model}?) and only these variables.
##'
##' @param npreg_method Method to use to calculate the EVPPI, for those methods that require it.    STATE SUPPORTED VALUES
##'
##' @param nsim Number of simulations from the model to use for calculating EVPPI.  The first \code{nsim} rows of the objects in \code{inputs} and \code{outputs} are used. 
##'
##' @param ... Other arguments required by specific methods
##'
##' @references
##'
##' Strong, M., Oakley, J. E., Brennan, A., & Breeze, P. (2015). Estimating the expected value of sample information using the probabilistic sensitivity analysis sample: a fast, nonparametric regression-based method. Medical Decision Making, 35(5), 570-583.
##' 
##' Menzies, N. A. (2016). An efficient estimator for the expected value of sample information. Medical Decision Making, 36(3), 308-320.
##'
##' Heath, A., Manolopoulou, I., & Baio, G. (2018). Efficient Monte Carlo estimation of the expected value of sample information using moment matching. Medical Decision Making, 38(2), 163-173.
##'
##' Jalal, H., & Alarid-Escudero, F. (2018). A Gaussian approximation approach for value of information analysis. Medical Decision Making, 38(2), 174-188.
##' 
##'
##' @export
evsi <- function(outputs,
                 inputs,
                 datagen_fn,
                 n=100,
                 method=NULL, # TODO speficy gam here or npreg? 
                 likelihood=NULL,
                 analysis_model=NULL,
                 model=NULL,
                 Q=30,
                 poi=NULL,
                 npreg_method="gam",
                 nsim=NULL,
                 verbose=TRUE, 
                 ...)
{
    check_inputs(inputs)
    output_type <- check_outputs(outputs, inputs)
    if (is.null(method))
        method <- default_evsi_method()

    ## Take subset of full PSA sample 
    if (is.null(nsim)) nsim <- nrow(inputs)
    if (output_type == "nb") {
        outputs <- outputs[1:nsim,,drop=FALSE]
    } else if (output_type == "cea") {
        outputs$c <- outputs$c[1:nsim,,drop=FALSE]
        outputs$e <- outputs$e[1:nsim,,drop=FALSE]
    }
    inputs <- inputs[1:nsim,,drop=FALSE]

    ## Could use any nonparametric regression method to regress on a summary statistic, in identical way to EVPPI estimation.
    
    if (method %in% npreg_methods) { 
        evsi_npreg(outputs=outputs, inputs=inputs, output_type=output_type,
                   datagen_fn=datagen_fn, n=n, likelihood=likelihood, method=method, verbose=verbose, ...)
    } else if (method=="is") {
        evsi_is(outputs=outputs, inputs=inputs, output_type=output_type,
                poi=poi, datagen_fn=datagen_fn, n=n, likelihood=likelihood,
                npreg_method=npreg_method, verbose=verbose, ...)
    } else if (method=="mm") {
        evsi_mm(outputs=outputs, inputs=inputs, output_type=output_type,
                poi=poi, datagen_fn=datagen_fn, n=n, Q=Q, 
                analysis_model=analysis_model,
                model=model, 
                npreg_method=npreg_method,
                verbose=verbose, ...)
    }
    else stop("Other methods not implemented yet")
}

evsi_npreg <- function(outputs, inputs, output_type, datagen_fn, n, method=NULL, verbose, ...){
    Tdata <- generate_data(inputs, datagen_fn, n)
    if (output_type == "nb")
        evppi_npreg_nb(nb=outputs, inputs=Tdata, poi=names(Tdata), method=method, verbose=verbose, ...)
    else if (output_type == "cea")
        evppi_npreg_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                        inputs=Tdata, poi=names(Tdata), method=method, verbose=verbose, ...)
}

generate_data <- function(inputs, datagen_fn, n=150){
    check_datagen_fn(datagen_fn, inputs)
    datagen_fn(inputs=inputs, n=n)
}

default_evsi_method <- function(){
    "gam" # TODO think about this 
}

check_datagen_fn <- function(datagen_fn, inputs){
    if (!is.function(datagen_fn)) stop("`datagen_fn` should be a function")
    ## If there's more than one argument, check that those args have
    ## defaults (e.g. sample sizes). Give error for now if not, but
    ## consider relaxing if this becomes a problem
    extra_args <- formals(datagen_fn)[-1]
    if (length(extra_args)>0){
        no_defaults <- sapply(extra_args, is.symbol)
        if (any(no_defaults)){
            stop(sprintf("Arguments \"%s\" of `datagen_fn` do not have default values", paste(names(no_defaults),collapse=",")))
        }
    } 
    ret <- datagen_fn(inputs)
    if (!is.data.frame(ret)) stop("`datagen_fn` should return a data frame")
    parnames <- names(ret)[names(ret) %in% names(inputs)]
    if (length(parnames)>0) {
        stop(sprintf("`datagen_fn` returns variables with the same names as parameters (%s).  It should return simulated data", paste(parnames, collapse=",")))
    }
    if (nrow(ret) != nrow(inputs)){
        stop("`datagen_fn` returns a data frame with %s rows. There should be %s rows, the same number of rows as `inputs`")
    }
}
