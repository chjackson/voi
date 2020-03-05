
## Importance sampling method for calculating EVSI (Menzies)

evsi_is <- function(outputs, inputs, output_type, poi, datagen_fn, n=100, likelihood, npreg_method="gam", verbose, ...){
    if (output_type=="nb"){
        evsi_is_nb(nb=outputs, inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                   likelihood=likelihood, npreg_method=npreg_method, verbose=verbose, ...)
    }
    else if (output_type=="cea"){
        evsi_is_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                    inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                    likelihood=likelihood, npreg_method=npreg_method, verbose=verbose, ...)
    }
}

evsi_is_nb <- function(nb, inputs, poi, datagen_fn, n, likelihood, npreg_method, verbose, ...){
    nbfit <- prepost_evsi_is(nb, inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                             likelihood=likelihood, npreg_method=npreg_method,
                             verbose=verbose, ...)
    calc_evppi(nbfit)
}

evsi_is_cea <- function(costs, effects, wtp, inputs, poi, datagen_fn, n, likelihood, npreg_method, verbose, ...){ 
    nwtp <- length(wtp)
    res <- numeric(nwtp)
    cfit <- prepost_evsi_is(costs, inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                            likelihood=likelihood, npreg_method=npreg_method, verbose=verbose, ...)
    efit <- prepost_evsi_is(effects, inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                            likelihood=likelihood, npreg_method=npreg_method, verbose=verbose, ...)
    for (i in 1:nwtp){
        nbfit <- efit*wtp[i] - cfit
        res[i] <- calc_evppi(nbfit)
    }
    res
}

prepost_evsi_is <- function(nb, inputs, poi, datagen_fn, n=100, likelihood, npreg_method="gam", verbose, ...){
    simdat <- generate_data(inputs, datagen_fn=datagen_fn, n=n, poi=poi)
    ## nb or ce? 
    if (is.null(npreg_method))
        npreg_method <- default_evppi_method(poi)
    if (verbose) cat("Calculating EVPPI...\n")
    y <- fitted_npreg(nb, inputs=inputs, poi=poi, method=npreg_method, verbose=verbose, ...)

    if (verbose) cat("Calculating EVSI...\n")
    nsam <- nrow(inputs)
    nout <- ncol(y) # TODO make sure 1D handled 
    prepost <- matrix(nrow=nsam, ncol=nout)
    for (i in 1:nsam){
        ll <- likelihood(simdat[i,], inputs, poi=poi) # vector of length nsim 
        w <- ll/sum(ll)
        for (j in 1:nout) { 
            ## to vectorise would need nsam x nsam storage 
            ## could do in C? though could we still work with user's lik fn?
            prepost[i,j] <- w %*% y[,j]
## ultra-verbose option
#            if (verbose)
#                message(sprintf("sample %s/%s, w1=%s, w2=%s", i, nsam, w[1],w[2]))
        }
    }
    prepost
}


form_likelihood <- function(study, likelihood, inputs, datagen_fn, poi){
    if (!is.null(study))
        likelihood <- get(sprintf("likelihood_%s", study))
    else {
        if (is.null(likelihood)) stop("`likelihood` should be supplied for method=\"is\"")
        if (!is.function(likelihood)) stop("`likelihood` should be a function")
        formals(likelihood) <- c(formals(likelihood), list(poi=NULL))
        check_likelihood(likelihood, inputs, datagen_fn, poi)
    }
    likelihood
}

check_likelihood <- function(likelihood, inputs, datagen_fn, poi){
    ## check that when likelihood is called with first two arguments data frames with
    ## 1. names matching output of datagen_fn
    ## 2. names matching inputs 
    ## returns output: vector length equal to nrow(inputs)
    data_sim <- datagen_fn(inputs, poi=poi)
    ret <- likelihood(data_sim, inputs)
    if (!is.vector(ret) | !is.numeric(ret))
        stop("likelihood function should return a numeric vector")
    if (length(ret) != nrow(inputs))
        stop(sprintf("likelihood function returns a vector of length %s, should be length %s, the number of rows in `inputs`", length(ret), nrow(inputs)))
}
