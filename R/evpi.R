##' Calculate the expected value of perfect information from a decision model
##' 
##' Calculate the expected value of perfect information from a decision model using standard Monte Carlo simulation
##'
##' @inheritParams evppi
##'
##' @return The expected value of perfect information, either a single value, or a vector of values for different willingness-to-pay.
##' 
##' @export
evpi <- function(outputs,
                 nsim=NULL)
{
    outputs <- check_outputs(outputs)
    if (is.null(nsim)) nsim <- if (inherits(outputs, "nb")) nrow(outputs) else nrow(outputs$e)
    outputs <- subset_outputs(outputs, nsim)
    if (inherits(outputs, "nb")){
        res <- mean(apply(outputs, 1, max)) - max(colMeans(outputs))
    } else if (inherits(outputs, "cea")){
        nwtp <- length(outputs$k)
        res <- numeric(length(nwtp))
        for (i in 1:nwtp){
            nb <- outputs$e * outputs$k[i]  -  outputs$c
            res[i] <- mean(apply(nb, 1, max)) - max(colMeans(nb))
        }
    }
    res
}
