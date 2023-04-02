
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
                                 datagen_fn=datagen_fn, n=n, aux_pars=aux_pars,
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

prepost_evsi_is <- function(out, inputs, pars, pars_datagen,
                            datagen_fn, n=100, 
                             aux_pars=aux_pars, 
                            likelihood, npreg_method="gam", verbose, ...){
    if (verbose) message("Generating data...")
    simdat <- generate_data(inputs, datagen_fn=datagen_fn, n=n, pars=pars_datagen, aux_pars=aux_pars)
    if (is.null(npreg_method))
        npreg_method <- default_evppi_method(pars)
    if (verbose) message("Calculating EVPPI...")
    y <- fitted_npreg(out, inputs=inputs, pars=pars, method=npreg_method, se=FALSE, verbose=verbose, ...)

    if (verbose) message("Calculating EVSI...")
    nsam <- nrow(inputs)
    nout <- ncol(y) # this doesn't handle 1D (for evsivar?)
    prepost <- matrix(nrow=nsam, ncol=nout)
    pb <- progress::progress_bar$new(total = nsam) 
    for (i in 1:nsam){
      ll <- eval_likelihood(likelihood, Y=simdat[i,,drop=FALSE], inputs=inputs,
                            n=n, pars=pars, aux_pars=aux_pars)
        w <- ll/sum(ll)
        for (j in 1:nout) { 
            ## This is slow, though to vectorise would need nsam x nsam storage 
            ## Could do in C? though could we still work with user's lik fn?
            prepost[i,j] <- w %*% y[,j]
        }
      pb$tick()
    }
    prepost
}

eval_likelihood <- function(likelihood, Y, inputs, n=100, pars, aux_pars=NULL){
    args <- list(Y=Y, inputs=inputs, n=n, pars=pars)
    args <- c(args, aux_pars)
    do.call(likelihood, args)
}

form_likelihood <- function(study, likelihood, inputs, datagen_fn, pars){
    if (!is.null(study)){
      if (!is.null(likelihood))
        warning("Ignoring `likelihood`, since a built-in study design was requested")
      likelihood <- get(sprintf("likelihood_%s", study))
    }
    else {
        if (is.null(likelihood)) stop("`likelihood` should be supplied for method=\"is\"")
        if (!is.function(likelihood)) stop("`likelihood` should be a function")
        if (length(formals(likelihood)) < 2) stop("`likelihood` should have at least two arguments")
        if (!identical(names(formals(likelihood))[1:2], c("Y","inputs")))
          stop("The first two arguments of `likelihood` should be named `Y` and `inputs`")
        if (!("n" %in% names(formals(likelihood))))
          formals(likelihood) <- c(formals(likelihood), list(n=100))
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
    ret <- likelihood(data_sim[1,,drop=FALSE], inputs=inputs)
    if (!is.vector(ret) | !is.numeric(ret))
        stop("likelihood function should return a numeric vector")
    if (length(ret) != nrow(inputs))
        stop(sprintf("likelihood function returns a vector of length %s, should be length %s, the number of rows in `inputs`", length(ret), nrow(inputs)))
}
