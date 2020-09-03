##' Calculate the expected value of sample information for an estimation problem
##'
##' Calculate the expected value of sample information for an estimation problem.  This computes the expected reduction in variance in some quantity of interest from a study of a certain design that informs the parameters of interest.
##'
##' @param outputs a vector of values for the quantity of interest, sampled from the uncertainty distribution of this quantity that is induced by the uncertainty about the parameters.
##'
##' @param method See \code{\link{evsi}}, only nonparametric regression methods are
##' currently supported in \code{\link{evsivar}}.
##'
##' @inheritParams evsi
##'
##' @references Jackson, 
##' Jackson, C., Presanis, A., Conti, S., & De Angelis, D. (2019). Value of information:
##' Sensitivity analysis and research design in Bayesian evidence synthesis.
##' Journal of the American Statistical Association, 114(528), 1436-1449.
##'
##' @export
evsivar <- function(outputs,
                 inputs,
                 study=NULL,
                 datagen_fn=NULL,
                 poi=NULL,
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
    inputs <- check_inputs(inputs, iname=deparse(substitute(inputs)))
    check_outputs_vector(outputs, inputs)
    if (is.null(method))
        method <- default_evsi_method()

    ## Take subset of full PSA sample 
    if (is.null(nsim)) nsim <- nrow(inputs)
    outputs <- outputs[1:nsim]
    inputs <- inputs[1:nsim,,drop=FALSE]

    datagen_fn <- form_datagen_fn(study, datagen_fn, inputs)
    if (method %in% npreg_methods) { 
        evsivar_npreg(outputs=outputs, inputs=inputs, 
                      datagen_fn=datagen_fn, poi=poi, n=n, 
                      method=method, verbose=verbose, ...)
    } 
    else stop("Other methods not implemented yet")
}

evsivar_npreg <- function(outputs, inputs, output_type, datagen_fn, poi, n, method=NULL, verbose, ...){
    Tdata <- generate_data(inputs, datagen_fn, n, poi)
    evppivar_npreg(outputs=outputs, inputs=Tdata, poi=names(Tdata), 
                   method=method, verbose=verbose, ...)
}

generate_data <- function(inputs, datagen_fn, n=150, poi){
    check_datagen_fn(datagen_fn, inputs, poi)
    datagen_fn(inputs=inputs, n=n, poi=poi)
}
