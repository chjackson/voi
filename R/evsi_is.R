
## Importance sampling method for calculating EVSI (Menzies)

evsi_is <- function(outputs, inputs, pars, pars_datagen,
                    datagen_fn, n=100, aux_pars=aux_pars,
                    likelihood, npreg_method="gam", verbose, ...){
    nn <- length(n)
    res <- vector(nn, mode="list")
    for (i in seq_along(n)){
        res[[i]] <- data.frame(
            n = n[i], 
          evsi = evsi_is_singlen(outputs, inputs, pars, pars_datagen,
                                 datagen_fn=datagen_fn, n=100, aux_pars=aux_pars,
                                 likelihood=likelihood, npreg_method=npreg_method, verbose=verbose,
                                 ...)
        )
        if (inherits(outputs, "cea")) res[[i]] <- cbind(k = outputs$k, res[[i]])
    }
    do.call("rbind", res)
}

evsi_is_singlen <- function(outputs, inputs, pars, pars_datagen,
                            datagen_fn, n=100, aux_pars=aux_pars, likelihood, npreg_method="gam", verbose, ...){
    UseMethod("evsi_is", outputs)    
}
    
evsi_is.nb <- function(outputs, inputs, pars, pars_datagen,
                       datagen_fn, n, aux_pars, likelihood, npreg_method, verbose, ...){
  nbfit <- prepost_evsi_is(outputs, inputs=inputs, pars=pars, pars_datagen=pars_datagen,
                           datagen_fn=datagen_fn, n=n,
                           aux_pars=aux_pars, 
                           likelihood=likelihood, npreg_method=npreg_method,
                           verbose=verbose, ...)
    calc_evppi(nbfit)
}

evsi_is.cea <- function(outputs, inputs, pars, pars_datagen,
                        datagen_fn, n, aux_pars, likelihood, npreg_method, verbose, ...){ 
    wtp <- outputs$k 
    nwtp <- length(wtp)
    res <- numeric(nwtp)
    cfit <- prepost_evsi_is(outputs$c, inputs=inputs, pars=pars, pars_datagen=pars_datagen,
                            datagen_fn=datagen_fn, n=n,
                             aux_pars=aux_pars, 
                            likelihood=likelihood, npreg_method=npreg_method, verbose=verbose, ...)
    efit <- prepost_evsi_is(outputs$e, inputs=inputs, pars=pars, pars_datagen=pars_datagen,
                            datagen_fn=datagen_fn, n=n,
                            aux_pars=aux_pars, 
                            likelihood=likelihood, npreg_method=npreg_method, verbose=verbose, ...)
    for (i in 1:nwtp){
        nbfit <- efit*wtp[i] - cfit
        res[i] <- calc_evppi(nbfit)
    }
    res
}

prepost_evsi_is <- function(nb, inputs, pars, pars_datagen,
                            datagen_fn, n=100, 
                             aux_pars=aux_pars, 
                            likelihood, npreg_method="gam", verbose, ...){
    simdat <- generate_data(inputs, datagen_fn=datagen_fn, n=n, pars=pars_datagen, aux_pars=aux_pars)
    ## nb or ce? 
    if (is.null(npreg_method))
        npreg_method <- default_evppi_method(pars)
    if (verbose) cat("Calculating EVPPI...\n")
    y <- fitted_npreg(nb, inputs=inputs, pars=pars, method=npreg_method, se=FALSE, verbose=verbose, ...)

    if (verbose) cat("Calculating EVSI...\n")
    nsam <- nrow(inputs)
    nout <- ncol(y) # TODO make sure 1D handled 
    prepost <- matrix(nrow=nsam, ncol=nout)
    for (i in 1:nsam){
        ## TODO handle aux_pars 
        ll <- likelihood(simdat[i,,drop=FALSE], inputs, pars=pars) # vector of length nsim 
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


form_likelihood <- function(study, likelihood, inputs, datagen_fn, pars){
    if (!is.null(study))
        likelihood <- get(sprintf("likelihood_%s", study))
    else {
        if (is.null(likelihood)) stop("`likelihood` should be supplied for method=\"is\"")
        if (!is.function(likelihood)) stop("`likelihood` should be a function")
        formals(likelihood) <- c(formals(likelihood), list(pars=NULL))
        check_likelihood(likelihood, inputs, datagen_fn, pars)
    }
    likelihood
}

check_likelihood <- function(likelihood, inputs, datagen_fn, pars){
    ## check that when likelihood is called with first two arguments data frames with
    ## 1. names matching output of datagen_fn
    ## 2. names matching inputs 
    ## returns output: vector length equal to nrow(inputs)
    data_sim <- datagen_fn(inputs, pars=pars)
    ret <- likelihood(data_sim, inputs=inputs)
    if (!is.vector(ret) | !is.numeric(ret))
        stop("likelihood function should return a numeric vector")
    if (length(ret) != nrow(inputs))
        stop(sprintf("likelihood function returns a vector of length %s, should be length %s, the number of rows in `inputs`", length(ret), nrow(inputs)))
}
