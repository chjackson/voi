
fitted_gam <- function(y, inputs, pars, ...){
    opts <- list(...)
    gam_formula <- opts$gam_formula
    if (is.null(gam_formula))
        gam_formula <- default_gam_formula(pars)
    gam_formula <- formula(sprintf("y ~ %s", gam_formula))
    model <- mgcv::gam(gam_formula, data = inputs)
    res <- model$fitted
    attr(res, "model") <- model
    res
}

default_gam_formula <- function(pars){
    karg <- if (length(pars) >=4) ", k=4" else ""
    sprintf("te(%s, bs='cr'%s)", paste(pars, collapse=", "), karg)
}

fitted_rep_gam <- function(model, B) { 
    beta_rep <- mvtnorm::rmvnorm(B, coef(model), vcov(model))
    fitted_rep <- beta_rep %*% t(predict(model,type="lpmatrix"))
    fitted_rep
}
