##' Calculate the expected value of sample information under a decision-analytic model
##'
##' Calculate the expected value of sample information under a decision-analytic model
##'
##' @inheritParams evppi
##' 
##' @param rfn Function to sample data from a proposed future study
##'
##' @param method Character string indicating the calculation method.
##'
##' Should support Strong, Heath, Jalal, Menzies at least (but named after method not person, as in convoi 1 paper)
##'
##' @param ... Other arguments required by specific methods 
##'
##' @export
evsi <- function(nb,
                  pars,
                  rfn,
                  method,
                  ...)
{
    check_nbpars(nb, pars)
}
