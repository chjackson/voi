##' Methods to calculate the Expected Value of Information 
##'
##' @description
##'
##' \code{\link{evppi}} calculates the expected value of partial perfect information from a decision-analytic model.  The default, recommended computation methods are based on nonparametric regression.  \code{\link{evpi}} is also provided for the expected value of perfect information.
##'
##' \code{\link{evsi}} calculates the expected value of sample information.   Currently this implements the same set of nonparametric regression methods as in \code{\link{evppi}}, and methods based on moment matching and importance sampling.  \code{\link{enbs}} can then be used to calculate and optimise the expected net benefit of sampling for a simple study with a fixed upfront cost and per-participant costs.
##'
##' \code{\link{evppi}} and \code{\link{evsi}} both require a sample of inputs and outputs from a Monte Carlo probabilistic analysis of a decision-analytic model.
##'
##' Analogous functions \code{\link{evppivar}} and \code{\link{evsivar}} calculate the EVPPI and EVSI for models used for estimation rather than decision-making.   The value of information is measured by expected reductions in variance of an uncertain model output of interest.
##'
##' A pure "brute-force" Monte Carlo method for EVPPI calculation is provided in \code{\link{evppi_mc}}, though this is usually computationally impractical.
##'
##' The \href{https://chjackson.github.io/voi/articles/voi.html}{package overview / Get Started vignette} gives worked examples of the use of all of these functions.
##'
##' @references
##'
##' Heath, A., Kunst, N., & Jackson, C. (eds.). (2024). Value of Information for Healthcare Decision-Making. CRC Press.
##'
##' Heath, A., Manolopoulou, I., & Baio, G. (2017). A review of methods for analysis of the expected value of information. Medical Decision Making, 37(7), 747-758.
##' 
##' Heath, A., Kunst, N., Jackson, C., Strong, M., Alarid-Escudero, F., Goldhaber-Fiebert, J. D., Baio, G. Menzies, N.A, Jalal, H. (2020). Calculating the Expected Value of Sample Information in Practice: Considerations from 3 Case Studies. Medical Decision Making, 40(3), 314-326.
##'
##' Kunst, N., Wilson, E. C., Glynn, D., Alarid-Escudero, F., Baio, G., Brennan, A., Fairley, M., Glynn, D., Goldhaber-Fiebert, J. D., Jackson, C., Jalal, H., Menzies, N. A., Strong, M., Thom, H., Heath, A. (2020). Computing the Expected Value of Sample Information Efficiently: Practical Guidance and Recommendations for Four Model-Based Methods. Value in Health, 3(6), 734-742.
##'
##' @name voi-package
##' @aliases voi-package
##'
##' @importFrom grDevices dev.off
##' @importFrom graphics points hist par
##' @importFrom ggplot2 ggplot aes geom_point ylab xlab geom_linerange geom_histogram
##' @importFrom gridExtra grid.arrange
##' @importFrom stats formula as.formula dist dnorm formula optim sd var rbeta rbinom rlnorm rnorm quantile dbinom coef vcov predict AIC cov lm fitted reorder IQR residuals update pnorm
##' @importFrom utils select.list combn
##' @importFrom progress progress_bar
##' @importFrom Matrix nearPD
##' 
"_PACKAGE"
