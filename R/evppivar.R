##' Calculate the expected value of partial perfect information for an estimation problem
##'
##' Calculate the expected value of partial perfect information for an estimation problem.  This computes the expected reduction in variance in some quantity of interest with perfect information about a parameter or parameters of interest.   
##'
##' @param outputs a vector of values for the quantity of interest, sampled from the uncertainty distribution of this quantity that is induced by the uncertainty about the parameters.
##'
##' Typically this will come from a Monte Carlo sample, where we first sample from the uncertainty distributions of the parameters, and then compute the quantity of interest as a function of the parameters.  It might also be produced by a Markov Chain Monte Carlo sample from the joint distribution of parameters and outputs. 
##'
##'
##' @inheritParams evppi
##' 
##' @references
##' Jackson, C., Presanis, A., Conti, S., & De Angelis, D. (2019). Value of information:
##' Sensitivity analysis and research design in Bayesian evidence synthesis.
##' Journal of the American Statistical Association, 114(528), 1436-1449.
##'
##' 
##' @export
evppivar <- function(outputs,
                  inputs,
                  pars=NULL,
                  method=NULL,
                  nsim=NULL,
                  verbose=TRUE,
                  ...)
{
    inputs <- check_inputs(inputs, iname=deparse(substitute(inputs)))
    check_outputs_vector(outputs, inputs)
    pars <- check_pars(pars, inputs)
    opts <- list(...)
    if (is.null(method))
        method <- default_evppi_method(pars)

    if (is.null(nsim)) nsim <- nrow(inputs)
    outputs <- outputs[1:nsim]
    inputs <- inputs[1:nsim,,drop=FALSE]
    
    if (method %in% npreg_methods) {
        evppivar_npreg(outputs=outputs, inputs=inputs, 
                    pars=pars, method=method, verbose=verbose, ...)
    } else stop("Other methods not implemented yet")
}

check_outputs_vector <- function(outputs, inputs){
    if (!is.numeric(outputs)) stop("`outputs` should be a numeric vector")
    if (length(outputs) != nrow(inputs))
        stop(sprintf("Length of `outputs` (%s) should equal the number of rows of `inputs` (%s)",
                     length(outputs), nrow(inputs)))
}

evppivar_npreg <- function(outputs, inputs, pars, method, verbose, ...){
    fitted <- fitted_npreg_fn(method)(y=outputs, inputs=inputs, pars=pars, 
                                      verbose=verbose, ...)
    var(fitted)
} 
