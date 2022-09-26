##' Calculate the expected value of sample information from a decision-analytic
##' model
##'
##' Calculate the expected value of sample information from a decision-analytic
##' model
##'
##' @inheritParams evppi
##'
##' @param study Name of one of the built-in study types supported by this
##'   package for EVSI calculation.  If this is supplied, then the columns of
##'   \code{inputs} that correspond to the parameters governing the study data
##'   should be identified in \code{pars}.
##'
##'   Current built-in studies are
##'
##'   \code{"binary"} A study with a binary outcome observed on one sample of
##'   individuals.   Requires one parameter: the probability of the outcome.
##'   The sample size is specifed in the \code{n} argument to \code{evsi()}, and
##'   the binomially-distributed outcome is named \code{X1}.
##'
##' \code{"trial_binary"} Two-arm trial with a binary outcome.   Requires two
##' parameters: the probability of the outcome in arm 1 and 2 respectively.  The
##' sample size is the same in each arm, specifed in the \code{n} argument to
##' \code{evsi()}, and the binomial outcomes are named \code{X1} and \code{X2}
##' respectively.
##'
##' \code{"normal_known"} A study of a normally-distributed outcome, with a
##' known standard deviation, on one sample of individuals.  Likewise the sample
##' size is specified in the \code{n} argument to \code{evsi()}.  The standard
##' deviation defaults to 1, and can be changed by specifying \code{sd} as a
##' component of the \code{aux_pars} argument, e.g. 
##' \code{evsi(..., aux_pars=list(sd=2))}.
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
##' @param pars Character vector identifying which columns of \code{inputs} are the parameters required to generate data from the proposed study.  Required if the proposed study is specified through the \code{study} argument, but not if it is specified through the \code{datagen_fn} argument.
##'
##' For example, if \code{study = "trial_binary"} is specified, then \code{pars} should be a vector of two elements naming the probability of the outcome in arm 1 and arm 2 respectively.
##' 
##' The \code{pars} argument is also required for the methods which involve an intermediate EVPPI calculation, that is, \code{method="is"} and \code{method="mm"}.  It should consist of the variables used in the definition of \code{datagen_fn} (and \code{likelihood}, \code{analysis_fn} and \code{model_fn}, if these are used) and only these variables.
##'
##' @param n Sample size of future study - optional argument to datagen_fn - facilitates calculating EVSI for multiple sample sizes.  If more than one quantity is required to describe the sample size (e.g. trials with unbalanced arms) currently you will have to write a different `datagen_fn` for each different sample size that you want the EVSI for.
##'
##' @param aux_pars A list of additional fixed arguments to supply to the function to generate the data, 
##' whether that is a built-in or user-defined function, e.g. \code{evsi(..., aux_pars = list(sd=2))} to change the 
##' fixed standard deviation in the \code{"normal_known"} model.  
##' 
##' @param method Character string indicating the calculation method.
##'
##' All the nonparametric regression methods supported for \code{\link{evppi}}, that is \code{"gam","gp","earth","inla"}, can also be used for EVSI calculation by regressing on a summary statistic of the predicted data (Strong et al 2015).   Defaults to \code{"gam"}.
##'
##' \code{"is"} for importance sampling (Menzies 2016)
##'
##' \code{"mm"} for moment matching (Heath et al 2018) (experimental and only partially implemented) 
##'
##' Note that the  \code{"is"} and \code{"mm"} methods are used in conjunction with nonparametric regression, thus the \code{gam_formula} argument can be supplied to \code{evsi} to specify this regression - see \code{\link{evppi}}. 
##'
##' @param likelihood Likelihood function, required (and only required) for the importance sampling method when a study design other than one of the built-in ones is used.  This should have two arguments as follows:
##'
##' 1. a data frame of predicted data. Columns are defined by the number of outcomes in the data, and names matching the data frame returned by \code{datagen_fn}. 
##'
##' 2. a data frame of parameter values, whose names should all correspond to variables in \code{inputs}.
##'
##' The function should return a vector whose length matches the number of rows of the parameters data frame given as the second argument.   Each element of the vector gives the likelihood of the corresponding set of parameters, given the data in the first argument.  An example is given in the vignette.
##'
##' Note the definition of the likelihood should agree with the definition of \code{datagen_fn} to define a consistent sampling distribution for the data.
##'
##' @param analysis_fn Function which fits a Bayesian model to the generated data.   Required for \code{method="mm"} if a study design other than one of the built-in ones is used.  This should be a function that takes the following arguments:
##'
##' `data`: A data frame with names matching the output of `datagen_fn`
##'
##' `args`: A list with constants required in the Bayesian analysis, e.g. prior parameters, or options for the analysis, e.g. number of MCMC simulations.
##'
##' `pars` Names of the parameters whose posterior is being sampled. 
##'
##' The function should return a data frame with names matching `pars`, containing a sample from the posterior distribution of the parameters given data supplied through `data`, and prior supplied through `args`. 
##'
##' @param analysis_args List of arguments required for the Bayesian analysis of the predicted data, e.g. definitions of the prior and options to control sampling.  Only used in \code{method="mm"}.  This is required if the study design is one of the built-in ones specified in \code{study}.  If a custom design is specifed through \code{analysis_fn}, then any constants needed in `analysis_fn` can either be supplied in `analysis_args`, or hard-coded in `analysis_fn` itself.
##'
##' For the built-in designs, the lists should have the following named components.  In addition for each study, to a component named `n` (the study sample size) should be supplied (TODO can't we read this from the evsi arg??).  An optional component `niter` in each case defines the posterior sample size (default 1000).
##'
##' `study="binary"`: `a` and `b`: Beta shape parameters 
##'
##' `study="trial_binary"`: `a1` and `b1`: Beta shape parameters for the prior for the first arm,  `a2` and `b2`: Beta shape parameters for the prior for the second arm. 
##'
##' `study="normal_known"`: `prior_mean`, `prior_sd` (prior mean and standard deviation) and `sampling_sd` (SD of an individual-level normal observation, so that the sampling SD of the mean outcome over the study is `sampling_sd/sqrt(n)`. 
##'
##' @param model_fn Function which evaluates the decision-analytic model, given parameter values.  Required for \code{method="mm"}.  See \code{\link{evppi_mc}} for full specification. 
##'
##' @param par_fn Function to simulate values from the uncertainty distributions of parameters needed by the decision-analytic model.  Should take one argument and return a data frame with one row for each simulated value, and one column for each parameter.  See \code{\link{evppi_mc}} for full specification. 
##'
##' @param Q Number of quantiles to use in \code{method="mm"}.
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
##' @export
evsi <- function(outputs,
                 inputs,
                 study=NULL,
                 datagen_fn=NULL,
                 pars=NULL,
                 n=100,
                 aux_pars=NULL, 
                 method=NULL,
                 likelihood=NULL,
                 analysis_fn=NULL,
                 analysis_args=NULL,
                 model_fn=NULL,
                 par_fn=NULL,
                 Q=30,
                 npreg_method="gam",
                 nsim=NULL,
                 verbose=FALSE, 
                 check=FALSE,
                 ...)
{
    check_inputs(inputs)
    check_ss(n)
    outputs <- check_outputs(outputs, inputs)
    check_pars(pars, inputs, evppi=FALSE)
    if (is.null(method))
        method <- default_evsi_method()

    if (is.null(nsim)) nsim <- nrow(inputs)
    outputs <- subset_outputs(outputs, nsim)
    inputs <- inputs[1:nsim,,drop=FALSE]

    datagen_fn <- form_datagen_fn(study, datagen_fn, inputs, aux_pars)
    ## Could use any nonparametric regression method to regress on a summary statistic, in identical way to EVPPI estimation.
    
    if (method %in% npreg_methods) { 
        res <- evsi_npreg(outputs=outputs, inputs=inputs, 
                   datagen_fn=datagen_fn, pars=pars, n=n, 
                   method=method, verbose=verbose, check=check, 
                   aux_pars = aux_pars, ...)
    } else if (method=="is") {
        likelihood <- form_likelihood(study, likelihood, inputs, datagen_fn, pars)
        res <- evsi_is(outputs=outputs, inputs=inputs, 
                pars=pars, datagen_fn=datagen_fn, n=n,
                aux_pars=aux_pars, likelihood=likelihood,
                npreg_method=npreg_method, verbose=verbose, ...)
    } else if (method=="mm") {
        res <- evsi_mm(outputs=outputs, inputs=inputs, 
                       pars=pars, datagen_fn=datagen_fn,
                       study=study, 
                       analysis_fn=analysis_fn,
                       analysis_args=analysis_args,
                       model_fn=model_fn, 
                       par_fn=par_fn, 
                       n=n, Q=Q,
                       npreg_method=npreg_method,
                       verbose=verbose, ...)
    }
    else stop("Other methods not implemented yet")
    attr(res, "method") <- method
    class(res) <- c("evsi", attr(res,"class"))
    res
}

evsi_npreg <- function(outputs, inputs, datagen_fn, pars, n, method=NULL, se=FALSE, B=500, verbose, check, aux_pars=NULL, ...){
    nn <- length(n)
    res <- vector(nn, mode="list")
    for (i in seq_along(n)){
        Tdata <- generate_data(inputs, datagen_fn, n[i], pars, aux_pars)
        res[[i]] <- evppi_npreg(outputs=outputs, inputs=Tdata, 
                           pars=names(Tdata), method=method, se=se, B=B, verbose=verbose, ...)
        names(res[[i]])[names(res[[i]])=="evppi"] <- "evsi"
        rownames(res[[i]]) <- NULL
    }
    resall <- do.call("rbind", res)
    resall <- cbind(n=n, resall)
    if (check){
      attr(resall, "models") <- lapply(res, function(x)attr(x, "models"))
      names(attr(resall,"models")) <- as.character(resall$n)
    }
    resall
}

generate_data <- function(inputs, datagen_fn, n=150, pars, aux_pars=NULL){
    check_datagen_fn(datagen_fn, inputs, pars)
    args <- list(inputs=inputs, n=n, pars=pars)
    args <- c(args, aux_pars)
    do.call(datagen_fn, args)
}

default_evsi_method <- function(){
    "gam"
}

check_study <- function(study) {
    if (!is.character(study) || (!(study %in% studies_builtin)))
      stop("``study` should be a character string matching one of the supported study designs")
}

form_datagen_fn <- function(study, datagen_fn, inputs, aux_pars=NULL){
  if (!is.null(study)) {
    check_study(study)
    datagen_fn <- get(sprintf("datagen_%s", study))
  } else {
    if (is.null(datagen_fn)) stop("`datagen_fn` should be supplied if `study` is not supplied")
    if (!is.function(datagen_fn)) stop("`datagen_fn` should be a function")
    formals(datagen_fn) <- c(formals(datagen_fn), list(pars=NULL))
    check_datagen_fn(datagen_fn, inputs, aux_pars)
  }
  datagen_fn
}

check_datagen_fn <- function(datagen_fn, inputs, pars=NULL, aux_pars=NULL){
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

form_analysis_fn <- function(study, analysis_fn){
  if (!is.null(study)){
    check_study(study)
    analysis_fn <- get(sprintf("analysis_%s", study))
  } else {
    if (is.null(analysis_fn)) stop("`analysis_fn` should be supplied if `study` is not supplied")
    if (!is.function(analysis_fn)) stop("`analysis_fn` should be a function")
#    formals(analysis_fn) <- c(formals(analysis_fn), list(pars=NULL))
#     check_analysis_fn(analysis_fn, inputs, aux_pars) # TODO 
  }
  analysis_fn
}

check_ss <- function(n){
    if (!is.numeric(n))
        stop("sample size `n` should be a numeric vector")
    if (any(n < 0))
        stop("sample sizes `n` should all be positive integers")
}

form_analysis_args <- function(analysis_args, study, n){
  if (!is.null(study)){
  if (study %in% studies_builtin){
    if (!is.list(analysis_args))
      stop("analysis_args should be supplied as a named list if using one of the built-in study designs")
    analysis_args$n <- n
  }
  }
  analysis_args
}
