
## TODO standard errors 

fitted_earth <- function(y, inputs, poi, ...){
    opts <- list(...)
    earth_formula <- opts$earth_formula
    if (is.null(earth_formula)) { 
        model <- earth::earth(y=y, x=inputs[,poi,drop=FALSE], ...)
    } else {
        earth_formula <- formula(sprintf("y ~ %s", earth_formula))
        model <- earth::earth(formula=earth_formula, data = inputs[,poi,drop=FALSE], ...)
    }
    model$fitted
}
