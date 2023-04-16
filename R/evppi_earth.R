
fitted_earth <- function(y, inputs, pars, verbose=FALSE, ...){
  opts <- list(...)
  earth_formula <- opts$earth_formula
  if (is.null(earth_formula)) {
    args <- list(y=y, x=inputs[,pars,drop=FALSE])
  } else {
    earth_formula <- formula(sprintf("y ~ %s", earth_formula))
    args <- list(formula=earth_formula, data = inputs[,pars,drop=FALSE])
  }
  if (!is.null(opts$se) && opts$se){
    if (is.null(opts$nfold)) opts$nfold <- 10 
    if (is.null(opts$ncross)) opts$ncross <- 30 
    if (is.null(opts$varmod.method)) opts$varmod.method <- "lm"
  }
  opts$se <- NULL
  model <- do.call(earth::earth, c(args, opts))
  model$y <- y
  res <- as.numeric(model$fitted)
  attr(res, "model") <- model
  res
}

check_plot_earth <- function(mod){
  graphics::par(mfrow=c(2,2))
  plot(mod)
}

check_stats_earth <- function(mod){
  list(gcv = mod$gcv)
}

fitted_rep_earth <- function(model, B) { 
  nobs <- length(model$fitted)
  fitted_rep <- matrix(nrow=B, ncol=nobs)
  se <- sqrt(as.numeric(model$varmod$model.var))
  for (i in 1:B){
    fitted_rep[i,] <- model$fitted + rnorm(nobs, mean=0, sd=se)
  }
  fitted_rep
}
