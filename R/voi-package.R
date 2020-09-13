##' Methods to calculate the Expected Value of Information 
##'
##' \code{\link{evppi}} calculates the expected value of partial perfect information.  Currently this implements a selection of nonparametric regression methods.
##'
##' \code{\link{evsi}} calculates the expected value of sample information.   Currently this implements the same set of nonparametric regression methods as in \code{\link{evppi}}, and a method based on importance sampling.
##'
##' @references
##'
##' Heath, A., Manolopoulou, I., & Baio, G. (2017). A review of methods for analysis of the expected value of information. Medical Decision Making, 37(7), 747-758.
##' 
##' Heath, A., Kunst, N., Jackson, C., Strong, M., Alarid-Escudero, F., Goldhaber-Fiebert, J. D., Baio, G. Menzies, N.A, Jalal, H. (2020). Calculating the Expected Value of Sample Information in Practice: Considerations from 3 Case Studies. Medical Decision Making, 40(3), 314-326.
##'
##' Kunst, N., Wilson, E. C., Glynn, D., Alarid-Escudero, F., Baio, G., Brennan, A., Fairley, M., Glynn, D., Goldhaber-Fiebert, J. D., Jackson, C., Jalal, H., Menzies, N. A., Strong, M., Thom, H., Heath, A. (2020). Computing the Expected Value of Sample Information Efficiently: Practical Guidance and Recommendations for Four Model-Based Methods. Value in Health, 3(6), 734-742.
##'
##' @name voi-package
##' @aliases voi-package
##' @docType package
##'
##' @importFrom grDevices dev.off
##' @importFrom graphics points
##' @importFrom stats formula as.formula dist dnorm formula optim sd var rbeta rbinom rlnorm rnorm quantile dbinom
##' @importFrom utils select.list
##' 
NULL
