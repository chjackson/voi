##' Calculate the expected value of partial perfect information from a decision-analytic model
##'
##' Calculate the expected value of partial perfect information from a decision-analytic model
##'
##' @param nb Matrix or data frame of samples from the uncertainty distribution of the expected net benefit.  The number of rows should equal the number of samples, and the number of columns should equal the number of decision options. 
##'
##' The exact names of all these arguments is up for discussion! 
##' 
##' (if we want to transform the NB to the INB to make certain calculation methods more efficient, we can do that internally. this would need a "baseline" treatment.  would the user need to supply this?  should default to first column at least)
##'
##' Should also be easy to compute for different willingness to pay values.  So should think about how to separate out costs and effects
##'
##' @param pars Matrix or data frame of samples from the uncertainty distribution of the parameters of the decision-analytic model.   The number of columns should equal the number of parameters, and the columns should be named.    This should have the same number of rows as \code{nb}, and each row of \code{nb} should give the net benefit function evaluated at the corresponding parameters.
##'
##' @param poi A character vector giving the parameters of interest, for which the EVPPI is required.   This should correspond to particular columns of \code{pars}.  [ I thought about allowing numeric indices for this, but perhaps better to encourage good practice by mandating that the parameters are named ] 
##'
##' @param method Character string indicating the calculation method.
##'
##' Should support gam, gp, inla at least 
##'
##' @param ... Other arguments required by specific methods 
##' 
##' @export
evppi <- function(nb,
                  pars,
                  poi,
                  method,
                  ...)
{
    check_nbpars(nb, pars)
    check_poi(poi, pars)
    nsam <- nrow(nb)
    npars <- nrow(pars)
}


check_nbpars <- function(nb, pars){
    if (!is.matrix(nb) && !is.data.frame(nb))
        stop("`nb` should be a matrix or data frame")
    if (!is.matrix(pars) && !is.data.frame(pars))
        stop("`pars` should be a matrix or data frame")
    if (nrow(nb) != nrow(pars))
        stop(sprintf("Number of rows of `nb` (%s) should equal the number of rows of `pars` (%s)", nrow(nb), nrow(pars)))
}

check_poi <- function(poi, pars){
    if (!is.character(poi))
        stop("`poi` should be a character vector")
    badpoi <- poi[!(poi %in% colnames(pars))]
    if (length(badpoi)>0){
        stop(sprintf("parameters of interest `%s` not found in columns of `pars`",
                     paste(badpoi,collapse=",")))
    }
}

