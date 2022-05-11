## TODO standard errors 

fitted_bart <- function(y, inputs, pars, verbose=FALSE, ...){
    opts <- list(...)
    model <- dbarts::bart(x.train=inputs[,pars,drop=FALSE], y.train=y, ...)
    model$y <- y
    model$rhat_mean <- check_bart_conv(model)
    res <- as.numeric(fitted(model))
    attr(res, "model") <- model
    res
}


## Convergence check for mean of fitted values
## Ideally should check convergence of EVPPI instead.
## For that, would need to extract MCMC sample of fitted values for each net benefit 
## then combine in list.   Expect good enough to do mean

check_bart_conv <- function(model){
    sam <- dbarts::extract(model) # 1000 MCMC samples for BART fit  x  nsam fitted values to evaluate convergence of
    sam.df <- data.frame(mean = rowMeans(sam))
    summ <- summary(posterior::as_draws(sam.df))
    summ$rhat
}
