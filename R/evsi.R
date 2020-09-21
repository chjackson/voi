##' Calculate the expected value of sample information from a decision-analytic model
##'
##' Calculate the expected value of sample information from a decision-analytic model
##'
##' @inheritParams evppi
##'
##' @param study Name of one of the built-in study types supported by this package for EVSI calculation.  If this is supplied, then the columns of \code{inputs} that correspond to the parameters governing the study data should be identified in \code{poi}.  
##'
##' Currently supported studies are
##'
##' \code{"trial_binary"} Two-arm trial with a binary outcome.   Requires two parameters: the probability of the outcome in arm 1 and 2 respectively.  The sample size is the same in each arm, specifed in the \code{n} argument to \code{evsi()}, and the binomial outcomes are named \code{X1} and \code{X2} respectively. 
##' 
##' Either \code{study} or \code{datagen_fn} should be supplied to \code{evsi()}.  
##'
##' 
##' @param datagen_fn If the proposed study is not one of the built-in types supported, it can be specified in this argument as an R function to sample predicted data from the study.  This function should have the following specification:
##'
##' 1. the function's first argument should be a data frame of parameter simulations, with one row per simulation and one column per parameter.  The parameters in this data frame must all be found in \code{inputs}.
##'
##' 2. the function should return a data frame.
##'
##' 3. the returned data frame should have number of rows equal to the number of parameter simulations in \code{inputs}.
##'
##' 4. if \code{inputs} is considered as a sample from the posterior, then \code{datagen_fn(inputs)} returns a corresponding sample from the posterior predictive distribution, which includes two sources of uncertainty: (a) uncertainty about the parameters and (b) sampling variation in observed data given fixed parameter values.
##'
##' 5. the function can optionally have more than one argument. If so, these additional arguments should be given default values in the definition of \code{datagen_fn}.  These arguments might be used to define sample sizes for a proposed study.
##'
##' Examples of this are currently in the \code{tests/tests_slow} package directory.
##'
##' @param pars Character vector identifying which columns of \code{inputs} are the parameters required to generate data from the proposed study.  Required if the proposed study is specified through the \code{study} argument, but not if it is specified through the \code{datagen_fn} argument.
##'
##' For example, if \code{study = "trial_binary"} is specified, then \code{pars} should be a vector of two elements naming the probability of the outcome in arm 1 and arm 2 respectively.
##' 
##' The \code{pars} argument is also required for the methods which involve an intermediate EVPPI calculation, that is the \code{"is"} and \code{"mm"}.  It should consist of the variables used in the definition of \code{datagen_fn} (and \code{likelihood} if used TODO ALSO in \code{analysis_model} and \code{model}?) and only these variables.
##'
##' @param n Sample size of future study - optional argument to datagen_fn - facilitates calculating EVSI for multiple sample sizes.  TODO if we want to design trials with multiple unbalanced arms, we'll need more than one argument. 
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
##' The Heath et al. and Jalal et al. methods are under development.
##'
##' @param likelihood Likelihood function, required (and only required) for the importance sampling method.  This should have two arguments as follows:
##'
##' 1. a data frame of predicted data. Columns are defined by the number of outcomes in the data, and names matching the data frame returned by \code{datagen_fn}. 
##'
##' 2. a data frame of parameter values, whose names should all correspond to variables in \code{inputs}.
##'
##' The function should return a vector whose length matches the number of rows of the parameters data frame given as the second argument.   Each element of the vector gives the likelihood of the corresponding set of parameters, given the data in the first argument.
##'
##' Examples of this are currently in \code{tests/tests_slow} and \code{tests/testthat} in the package directory. 
##'
##' Note the definition of the likelihood should agree with the definition of \code{datagen_fn} to define a consistent sampling distribution for the data.
##'
##' @param analysis_model Function which fits a Bayesian model to the generated data.   Under development (need to decide format, output, JAGS dependencies, etc.). Required for \code{method="mm"} (and Jalal method if n0 not given)
##'
##' @param analysis_options List of arguments required by \code{analysis_model}.  Under development - for \code{method="mm"} and Jalal method. 
##'
##' @param decision_model Function which evaluates the decision-analytic model, given parameter values.  Under development - for \code{method="mm"} and Jalal method.  Need to decide the required format for nb, c, e output.
##'
##' @param Q Number of quantiles to use in \code{method="mm"} (under development). 
##'
##' @param npreg_method Method to use to calculate the EVPPI, for those methods that require it.  This is passed to \code{\link{evppi}} as the \code{method} argument. 
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
                 study=NULL,
                 datagen_fn=NULL,
                 pars=NULL,
                 n=100,
                 method=NULL, # TODO speficy gam here or npreg? 
                 likelihood=NULL,
                 analysis_model=NULL,
                 analysis_options=NULL,
                 decision_model=NULL,
                 Q=30,
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

    datagen_fn <- form_datagen_fn(study, datagen_fn, inputs)
    ## Could use any nonparametric regression method to regress on a summary statistic, in identical way to EVPPI estimation.
    
    if (method %in% npreg_methods) { 
        evsi_npreg(outputs=outputs, inputs=inputs, output_type=output_type,
                   datagen_fn=datagen_fn, pars=pars, n=n, 
                   method=method, verbose=verbose, ...)
    } else if (method=="is") {
        likelihood <- form_likelihood(study, likelihood, inputs, datagen_fn, pars)
        evsi_is(outputs=outputs, inputs=inputs, output_type=output_type,
                pars=pars, datagen_fn=datagen_fn, n=n, likelihood=likelihood,
                npreg_method=npreg_method, verbose=verbose, ...)
    } else if (method=="mm") {
        evsi_mm(outputs=outputs, inputs=inputs, output_type=output_type,
                pars=pars, datagen_fn=datagen_fn, n=n, Q=Q, 
                analysis_model=analysis_model,
                analysis_options=analysis_options,
                decision_model=decision_model, 
                npreg_method=npreg_method,
                verbose=verbose, ...)
    }
    else stop("Other methods not implemented yet")
}

evsi_npreg <- function(outputs, inputs, output_type, datagen_fn, pars, n, method=NULL, se=FALSE, B=500, verbose, ...){
    Tdata <- generate_data(inputs, datagen_fn, n, pars)
    if (output_type == "nb")
        evppi_npreg_nb(nb=outputs, inputs=Tdata, pars=names(Tdata), method=method, se=se, B=B, verbose=verbose, ...)
    else if (output_type == "cea")
        evppi_npreg_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                        inputs=Tdata, pars=names(Tdata), method=method, se=se, B=B, verbose=verbose, ...)
}

generate_data <- function(inputs, datagen_fn, n=150, pars){
    check_datagen_fn(datagen_fn, inputs, pars)
    datagen_fn(inputs=inputs, n=n, pars=pars)
}

default_evsi_method <- function(){
    "gam" # TODO think about this 
}

form_datagen_fn <- function(study, datagen_fn, inputs){
    if (!is.null(study)){
        if (!is.character(study) || (!(study %in% studies_builtin)))
            stop("``study` should be a character string matching one of the supported study designs")
        else datagen_fn <- get(sprintf("datagen_%s", study))
    } else {
        if (is.null(datagen_fn)) stop("`datagen_fn` should be supplied if `study` is not supplied")
        if (!is.function(datagen_fn)) stop("`datagen_fn` should be a function")
        formals(datagen_fn) <- c(formals(datagen_fn), list(pars=NULL))
        check_datagen_fn(datagen_fn, inputs)
    }
    datagen_fn
}

check_datagen_fn <- function(datagen_fn, inputs, pars=NULL){
    ## If there's more than one argument, check that those args have
    ## defaults (e.g. sample sizes). Give error for now if not, but
    ## consider relaxing if this becomes a problem
    extra_args <- formals(datagen_fn)[-1]
    extra_args <- extra_args[names(extra_args) != "pars"]
    if (length(extra_args)>0){
        no_defaults <- sapply(extra_args, is.symbol)
        if (any(no_defaults)){
            stop(sprintf("Arguments \"%s\" of `datagen_fn` do not have default values", paste(names(no_defaults),collapse=",")))
        }
    }
    ret <- datagen_fn(inputs, pars=pars)
    if (!is.data.frame(ret)) stop("`datagen_fn` should return a data frame")
    parnames <- names(ret)[names(ret) %in% names(inputs)]
    if (length(parnames)>0) {
        stop(sprintf("`datagen_fn` returns variables with the same names as parameters (%s).  It should return simulated data", paste(parnames, collapse=",")))
    }
    if (nrow(ret) != nrow(inputs)){
        stop(sprintf("`datagen_fn` returns a data frame with %s rows. There should be %s rows, the same number of rows as `inputs`", nrow(ret), nrow(inputs)))
    }
}
