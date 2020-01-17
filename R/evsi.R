##' Calculate the expected value of sample information from a decision-analytic model
##'
##' Calculate the expected value of sample information from a decision-analytic model
##'
##' @inheritParams evppi
##' 
##' @param rfn Function to sample predicted data from a proposed future study.
##'
##' This should have first argument defined by a data frame of parameter simulations, with one row per simulation and one column per parameter.  The names of the parameters must all correspond to parameters in \code{inputs}.
##'
##' The function should return a data frame with the number of rows equal to the number of parameter simulations.  Thus if \code{inputs} is considered as a sample from the posterior, then \code{rfn} returns a corresponding sample from the posterior predictive distribution, which includes two sources of uncertainty: (a) uncertainty about the parameters and (b) sampling variation in observed data given fixed parameter values.
##'
##'
##' @param n Sample size of future study - optional argument to rfn - facilitates calculating EVSI for multiple sample sizes.  TODO how will this work for randomised trials with multiple arms? 
##'
##' @param method Character string indicating the calculation method.
##'
##' All the nonparametric regression methods supported for \code{\link{evppi}}, that is \code{"gam","gp","earth","inla"}, can also be used for EVSI calculation by regressing on a summary statistics of the predicted data (Strong et al 201?).   Defaults to \code{"gam"}.
##'
##' \code{"is"} for importance sampling (Menzies 200?)
##'
##' \code{"mm"} for moment matching (Heath 200?)
##'
##' TODO Heath and Jalal methods
##'
##' @param likelihood Likelihood function, required (and only required) for the importance sampling method.  This should take arguments:
##'
##' the first: a data frame with columns defined by the number of outcomes in the data, and with names matching the names of the data frame returned by \code{rfn}. 
##'
##' the second: a data frame of parameter values, whose names should all correspond to variables in \code{inputs}.
##'
##' The function should return a vector whose length matches the number of rows of the data frame given as the second argument. 
##'
##' POINT TO AN EXAMPLE WHICH WILL MAKE ALL THIS CLEARER
##'
##' Note the definition of the likelihood should agree with the definition of \code{rfn} to define a consistent sampling distribution for the data.   CLARIFY.  [ eventually we'll want some built-in common examples where people don't have to specify either rfn or likelihood ]
##'
##' @param analysis_model Function which fits a Bayesian model to the generated data. TODO work out format, output, JAGS dependencies, etc.  Required for \code{method="mm"} (and Jalal method if n0 not given??) 
##'
##' @param model Function which evaluates the decision-analytic model, given parameter values.  TODO sort out format for nb, c, e output
##'
##' @param Q Number of quantiles to use in \code{method="mm"}. 
##'
##' @param poi Parameters of interest, that is, those which are informed by the data in the future study.  Required (and only required) for the methods which involve an intermediate EVPPI calculation, that is the \code{"is"} and \code{"mm"} TODO OTHER methods.
##'
##' This should bee a character vector naming particular columns of \code{inputs}.  It should consist of the variables used in the definition of \code{rfn} (and \code{likelihood} if used TODO ALSO in \code{analysis_model} and \code{model}?) and only these variables.
##'
##' @param npreg_method Method to use to calculate the EVPPI, for those methods that require it.    STATE SUPPORTED VALUES
##'
##' @param nsim Number of simulations from the model to use for calculating EVPPI.  The first \code{nsim} rows of the objects in \code{inputs} and \code{outputs} are used. 
##'
##' @param ... Other arguments required by specific methods 
##'
##' @export
evsi <- function(outputs,
                 inputs,
                 rfn,
                 n=100,
                 method=NULL, # TODO speficy gam here or npreg? 
                 likelihood=NULL,
                 analysis_model=NULL,
                 model=NULL,
                 Q=30,
                 poi=NULL,
                 npreg_method="gam",
                 nsim=NULL,
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
                   rfn=rfn, n=n, likelihood=likelihood, method=method, ...)
    } else if (method=="is") {
        evsi_is(outputs=outputs, inputs=inputs, output_type=output_type,
                poi=poi, rfn=rfn, n=n, likelihood=likelihood, npreg_method=npreg_method, ...)
    } else if (method=="mm") {
        evsi_mm(outputs=outputs, inputs=inputs, output_type=output_type,
                poi=poi, rfn=rfn, n=n, Q=Q, 
                analysis_model=analysis_model,
                model=model, 
                npreg_method=npreg_method, ...)
    }
    else stop("Other methods not implemented yet")
}

evsi_npreg <- function(outputs, inputs, output_type, rfn, n, method=NULL, ...){
    Tdata <- generate_data(inputs, rfn, n)
    if (output_type == "nb")
        evppi_npreg_nb(nb=outputs, inputs=Tdata, poi=names(Tdata), method=method, ...)
    else if (output_type == "cea")
        evppi_npreg_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                        inputs=Tdata, poi=names(Tdata), method=method, ...)
}

generate_data <- function(inputs, rfn, n=150){
    ## TODO check rfn is a function with correct args, and check its output
    ## do we need this in a separate function anyway?  better name for rfn? 
    rfn(inputs=inputs, n=n)
}

default_evsi_method <- function(){
    "gam" # TODO think about this 
}
