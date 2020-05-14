## TODO standard errors 

fitted_gam <- function(y, inputs, poi, ...){
    opts <- list(...)
    gam_formula <- opts$gam_formula
    if (is.null(gam_formula))
        gam_formula <- default_gam_formula(poi)
    gam_formula <- formula(sprintf("y ~ %s", gam_formula))
    model <- mgcv::gam(gam_formula, data = inputs)
    model$fitted
}

default_gam_formula <- function(poi){
    karg <- if (length(poi) >=4) ", k=4" else ""
    sprintf("te(%s, bs='cr'%s)", paste(poi, collapse=", "), karg)
}
