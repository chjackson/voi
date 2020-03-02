evsi_mm <- function(outputs, inputs, output_type, poi, datagen_fn, analysis_model, model, n=100, likelihood, Q, npreg_method="gam", verbose, ...){
    if (output_type=="nb"){
        evsi_mm_nb(nb=outputs, inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                   Q=Q, 
                   analysis_model=analysis_model,
                   model=model, 
                   npreg_method=npreg_method,
                   verbose=verbose, ...)
    }
    else if (output_type=="cea"){
        evsi_mm_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                    inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                    Q=Q, 
                    analysis_model=analysis_model,
                    model=model, 
                    npreg_method=npreg_method,
                    verbose=verbose, ...)
    }
}

evsi_mm_nb <- function(nb, inputs, poi, datagen_fn, n, Q=Q, 
                    analysis_model,
                    model, 
                    npreg_method,
                    verbose, ...){

###   TODO from EVSI package

    ## Get quantiles of input parameters from PSA sample 

    ## Generate future data given Q quantiles

    ## Fit Bayesian model to get posterior

    ## Run HE model again with posterior sample

    ## Loop, calculate variance, -> EVSI.

}


evsi_mm_cea <- function(costs, effects, wtp, inputs, poi, datagen_fn, n, Q=Q, 
                    analysis_model,
                    model, 
                    npreg_method,
                    verbose, ...){
}
