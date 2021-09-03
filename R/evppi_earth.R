
## TODO standard errors 

fitted_earth <- function(y, inputs, pars, verbose=FALSE, ...){
    opts <- list(...)
    earth_formula <- opts$earth_formula
    if (is.null(earth_formula)) { 
        model <- earth::earth(y=y, x=inputs[,pars,drop=FALSE], ...)
    } else {
        earth_formula <- formula(sprintf("y ~ %s", earth_formula))
        model <- earth::earth(formula=earth_formula, data = inputs[,pars,drop=FALSE], ...)
    }
    as.numeric(model$fitted)
}
