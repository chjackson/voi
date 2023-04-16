
fitted_bart <- function(y, inputs, pars, verbose=FALSE, ...){
    opts <- list(...)
    model <- dbarts::bart(x.train=inputs[,pars,drop=FALSE], y.train=y, verbose=verbose, ...)
    model$y <- y
    model$rhat_mean <- check_bart_conv(model)
    res <- as.numeric(fitted(model))
    attr(res, "model") <- model
    res
}


## Convergence check for mean of fitted values, ie for 
## mubar = (average of mu|X over X in data)
## where mu is the expectation that BART estimates
## Assume this is sufficient to ensure EVPPI estimate has "converged" 

check_bart_conv <- function(model){
    sam <- dbarts::extract(model) # 1000 MCMC samples for BART fit  x  nsam fitted values to evaluate convergence of
    sam.df <- data.frame(mean = rowMeans(sam))
    summ <- summary(posterior::as_draws(sam.df))
    summ$rhat
}

fitted_rep_bart <- function(model) {
  as.matrix(dbarts::extract(model))
}
