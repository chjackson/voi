
npreg_methods <- c("gam", "gp", "inla", "earth", "bart")

evppi_npreg <- function(outputs, ...){
    UseMethod("evppi_npreg", outputs)
}

##' @noRd
evppi_npreg.nb <- function(outputs, inputs, pars, method, se, B, verbose, ...){
    if (verbose) message("Fitting nonparametric regression") 
    fit <- fitted_npreg(outputs, inputs=inputs, pars=pars, method=method, se=se, B=B, 
                        verbose=verbose, ...)
    res <- data.frame(evppi = calc_evppi(fit))
    if (se){
        if (verbose) message("Calculating replicate EVPPIs from simulation")
        evppi_rep <- numeric(B)
        for (i in 1:B){
            ## This bit could be faster.  Is there a vectorised way?  Naive apply doesn't help
            evppi_rep[i] <- calc_evppi(attr(fit, "rep")[i,,])   
        }
        res$se <-  sd(evppi_rep)
    }
    attr(res, "models") <- attr(fit, "models")
    res
}

##' @noRd
evppi_npreg.cea <- function(outputs, inputs, pars, method, se, B, verbose, ...){
    wtp <- outputs$k
    if (verbose) message("Fitting nonparametric regression for costs") 
    cfit <- fitted_npreg(outputs$c, inputs=inputs, pars=pars, method=method, se=se, B=B, verbose=verbose, ...)
    if (verbose) message("Fitting nonparametric regression for effects") 
    efit <- fitted_npreg(outputs$e, inputs=inputs, pars=pars, method=method, se=se, B=B, verbose=verbose, ...)
    calc_evppi_ce(cfit, efit, wtp, se=se, B=B, verbose=verbose)
}

calc_evppi_ce <- function(cfit, efit, wtp, se=FALSE, B=0, verbose=FALSE){
    nwtp <- length(wtp)
    res <- resse <- numeric(nwtp)
    for (i in 1:nwtp){
        evppi_rep <- numeric(B)
        inbfit <- efit*wtp[i] - cfit
        if (se){
            inbrep <- attr(efit, "rep")*wtp[i] - attr(cfit, "rep")
            if (verbose) message("Calculating replicate EVPPIs from simulation")
            for (j in 1:B){
                evppi_rep[j] <- calc_evppi(inbrep[j,,])   
            }
            resse[i] <- sd(evppi_rep)
        }
        res[i] <- calc_evppi(inbfit)
    }
    res <- data.frame(k=wtp, evppi=res)
    attr(res, "models") <- list(c=attr(cfit, "models"), e=attr(efit, "models"))
    if (se) res$se <- resse
    res
}

check_ref <- function(ref, nb){
  if (is.character(ref)){
    ref_num <- match(ref, colnames(nb))
    if (is.na(ref_num))
      stop(sprintf("reference decision option ref=\"%s\" does not appear in the supplied column names of `outputs`: %s",
                   ref,
                   paste(paste0("\"", colnames(nb), "\""),collapse=",")))
  }
  else {
    if (!(ref %in% 1:ncol(nb)))
      stop(sprintf("reference decision option `ref` should either be a string matching one of the column names of `outputs`, or an integer <= %s indicating the corresponding column number of `outputs`", ncol(nb)))
    ref_num <- ref
  }
  ref_num
}

fitted_npreg <- function(nb, inputs, pars, method, se=FALSE, B=NULL, verbose, ...){
    nopt <- ncol(nb)
    nsim <- nrow(nb)
    ## Transforming to incremental net benefit allows us to do one fewer regression
    ref <- list(...)$ref # reference decision option
    if (is.null(ref)) ref <- 1
    else ref <- check_ref(ref, nb)
    inb <- nb[, -ref, drop=FALSE] - nb[,ref]
    fitted <- matrix(0, nrow=nsim, ncol=nopt)
    if (se) {
        if (method=="bart"){
            B <- list(...)$ndpost
            if (is.null(B)) B <- formals(dbarts::bart)$ndpost
        }
        fitted_rep <- array(0, dim=c(B, nsim, nopt))
    }
    models <- vector(nopt-1, mode="list")
    for (i in 1:(nopt-1)){
        if (verbose) message(sprintf("Decision option %s",i+1)) 
        fit <- fitted_npreg_fn(method)(y=inb[,i], inputs=inputs, pars=pars, verbose=verbose, se=se, ...) 
        fitted[,i+1] <- as.vector(fit)
        if (se){
            fitted_rep[,,i+1] <- fitted_npreg_rep_call(method, attr(fit,"model"), B, verbose)
        }
        models[[i]] <- attr(fit, "model")
    }
    if (se) attr(fitted, "rep") <- fitted_rep
    attr(fitted, "models") <- models
    fitted
}

fitted_npreg_fn <- function(method){
    switch(method, 
           gam = fitted_gam,
           gp = fitted_gp,
           inla = fitted_inla,
           earth = fitted_earth,
           bart = fitted_bart)
}

fitted_npreg_rep_call <- function(method, model, B, verbose=FALSE){
    if (verbose) message("Simulating parameters to calculate standard errors")
    if (method=="gam") {
        frep <- fitted_rep_gam(model, B)
    } else if (method=="earth") {
        frep <- fitted_rep_earth(model, B)
    } else if (method=="gp") {
        frep <- fitted_rep_gp(model, B)
    } else if (method=="bart") {
        frep <- fitted_rep_bart(model)
    }
    else stop(sprintf("Standard errors not available for method = \"%s\"",method))
    frep
}

calc_evppi <- function(fit) {
    ## NAs are removed in BCEA here. Shouldn't we make users investigate them and remove by hand if they know what they are doing?  At least warn users if there are NAs in their samples
    mean(apply(fit, 1, max)) - max(colMeans(fit))
}

