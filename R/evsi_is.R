
## Importance sampling method for calculating EVSI (Menzies)

evsi_is <- function(outputs, inputs, output_type, poi, rfn, n=100, likelihood, npreg_method="gam", ...){
    if (output_type=="nb"){
        evsi_is_nb(nb=outputs, inputs=inputs, poi=poi, rfn=rfn, n=n,
                   likelihood=likelihood, npreg_method=npreg_method, ...)
    }
    else if (output_type=="cea"){
        evsi_is_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                    inputs=inputs, poi=poi, rfn=rfn, n=n,
                    likelihood=likelihood, npreg_method=npreg_method, ...)
    }
}

evsi_is_nb <- function(nb, inputs, poi, rfn, n, likelihood, npreg_method, ...){
    nbfit <- prepost_evsi_is(nb, inputs=inputs, poi=poi, rfn=rfn, n=n,
                             likelihood=likelihood, npreg_method=npreg_method, ...)
    calc_evppi(nbfit)
}

evsi_is_cea <- function(costs, effects, wtp, inputs, poi, rfn, n, likelihood, npreg_method, ...){ 
    nwtp <- length(wtp)
    res <- numeric(nwtp)
    cfit <- prepost_evsi_is(costs, inputs=inputs, poi=poi, rfn=rfn, n=n,
                            likelihood=likelihood, npreg_method=npreg_method, ...)
    efit <- prepost_evsi_is(effects, inputs=inputs, poi=poi, rfn=rfn, n=n,
                            likelihood=likelihood, npreg_method=npreg_method, ...)
    for (i in 1:nwtp){
        nbfit <- efit*wtp[i] - cfit
        res[i] <- calc_evppi(nbfit)
    }
    res
}

prepost_evsi_is <- function(nb, inputs, poi, rfn, n=100, likelihood, npreg_method="gam", ...){
    simdat <- generate_data(inputs, rfn=rfn, n=n)
    ## nb or ce? 
    if (is.null(npreg_method))
        npreg_method <- default_evppi_method(poi)
    y <- fitted_npreg(nb, inputs=inputs, poi=poi, method=npreg_method, ...)
    cat("done evppi\n")
    
    nsam <- nrow(inputs)
    nout <- ncol(y) # TODO make sure 1D handled 
    prepost <- matrix(nrow=nsam, ncol=nout)
    for (i in 1:nsam){
        ll <- likelihood(simdat[i,], inputs) # vector of length nsim 
        w <- ll/sum(ll)
        for (j in 1:nout) { 
            ## to vectorise would need nsam x nsam storage 
            ## could do in C? though could we still work with user's lik fn?
            prepost[i,j] <- w %*% y[,j]
            cat(sprintf("w1=%s, w2=%s, [%s]\n",w[1],w[2],i))
        }
    }
    prepost
}
