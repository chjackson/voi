## TODO standard errors 

fitted_gam <- function(y, inputs, pars, ...){
    opts <- list(...)
    gam_formula <- opts$gam_formula
    if (is.null(gam_formula))
        gam_formula <- default_gam_formula(pars)
    gam_formula <- formula(sprintf("y ~ %s", gam_formula))
    model <- mgcv::gam(gam_formula, data = inputs)
    model$fitted
}

default_gam_formula <- function(pars){
    karg <- if (length(pars) >=4) ", k=4" else ""
    sprintf("te(%s, bs='cr'%s)", paste(pars, collapse=", "), karg)
}
