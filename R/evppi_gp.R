## Gaussian process method for estimating EVPPI
## Code from SAVI 
## [Aren't there also packaged GP regression methods ?  kernlab?  gaupro, GPfit? are they more efficient?]

fitted_gp <- function(y, inputs, pars, verbose=FALSE, ...){
    args <- list(...)
    model <- gp(y=y, X=inputs[,pars], m=args$gp_hyper_n, maxSample=args$maxSample, verbose=verbose)
    res <- model$fitted
    attr(res, "model") <- model
    res 
}

fitted_rep_gp <- function(model, B, ...){
  mvtnorm::rmvnorm(B, model$fitted, model$V)
}

## Code below adapted from SAVI, copyright (c) 2014, 2015 the SAVI authors
## https://github.com/Sheffield-Accelerated-VoI/SAVI-package/blob/master/inst/SAVI/scripts_GPfunctions.R
## * Changes to make it more clean/modular. 
## * Check for constant/collinearity temporarily moved out 
## * Added gp_hyper_n option to select the number of iterations used for estimating the hyperparameters.  BCEA [code from older SAVI version?] uses all by default, while SAVI uses a small number for efficiency.  SAVI value taken as the default here. 

dinvgamma <- function(x, alpha, beta) {
  (beta ^ alpha) / gamma(alpha) * x ^ (-alpha - 1) * exp(-beta / x)
}

cor.Gaussian <- function(X, phi, m) { #  needs sorting...... probably with dist()
  txbuild1 <- function(h) exp(-rowSums(t((t(X) - h) / phi) ^ 2))
  apply(as.matrix(as.matrix(X)[1:m, ]), 1, txbuild1)
}

makeA.Gaussian <- function(X, phi) { # function to make A matrix with the Gaussian correlation function
  n <- NROW(X)
  if(length(phi) > 1) {
    b <- diag(phi ^ (-2))
  } else {
    b <- phi ^ (-2) }
  R <- X %*% as.matrix(b) %*% t(X)
  
  S <- matrix(diag(R), nrow = n, ncol = n)
  exp(2 * R - S - t(S))
}

# function to calculate posterior density


post.density <- function(hyperparams, NB, input.m) {
  
  input.m <- as.matrix(input.m, drop = FALSE)
  
  N <- nrow(input.m)
  p <- NCOL(input.m)
  H <- cbind(1, input.m)
  q <- ncol(H)
  
  a.sigma <- 0.001; b.sigma <- 0.001  ##  hyperparameters for IG prior for sigma^2
  a.nu <- 0.001; b.nu <- 1            ##  hyperparameters for IG prior for nu
  delta <- exp(hyperparams)[1:p]
  nu <- exp(hyperparams)[p + 1]
  
  A <- makeA.Gaussian(input.m, delta)
  Astar <- A + nu * diag(N)
  T <- chol(Astar)
  y <- backsolve(t(T), NB, upper.tri = FALSE)
  x <- backsolve(t(T), H, upper.tri = FALSE)
  tHAstarinvH <- t(x) %*% (x) + 1e-7* diag(q)
  betahat <- solve(tHAstarinvH) %*% t(x) %*% y
  residSS <- y %*% y -t(y) %*% x %*% betahat - t(betahat) %*% t(x) %*% y +
    t(betahat) %*% tHAstarinvH %*% betahat
  prior <- prod(dnorm(log(delta), 0, sqrt(1e5))) * dinvgamma(nu, a.nu, b.nu)
  l <- -sum(log(diag(T))) - 1 / 2 * log(det(tHAstarinvH)) -
    (N - q + 2 * a.sigma) / 2 * log(residSS / 2 + b.sigma) + log(prior)
  return(l)
}


estimate.hyperparameters <- function(NB, inputs, verbose=FALSE) {
  p <- NCOL(inputs)
  hyperparameters <- NA
    initial.values <- rep(0, p + 1)
    repeat {
        if (verbose) 
            print(paste("calling optim function for net benefit"))
      log.hyperparameters <- optim(initial.values, fn=post.density, 
                                   NB=NB, input.m=inputs,
                                   method="Nelder-Mead",
                                   control=list(fnscale=-1, maxit=10000, trace=0))$par
      if (sum(abs(initial.values - log.hyperparameters)) < 0.05) {
        hyperparameters <- exp(log.hyperparameters)
        break
      }	
      initial.values <- log.hyperparameters
    }
  return(hyperparameters)
}


##' Fit a Gaussian process regression
##'
##' @param y Vector of outcome data
##' 
##' @param X Matrix of inputs
##'
##' @param hyper Hyperparameter values. If set to \code{"est"} (the default), these are estimated. 
##'
##' @param m number of samples to use to estimate the hyperparameters. By default, this is the minimum
##' of the following three quantities: 30 times the number of predictors in \code{X}, 
##' interest, 250, and the length of \code{y}. 
##'
##' @param maxSample Maximum sample size to employ. If datasets larger than this are supplied, they are
##' truncated to this size.
##'
##' @param verbose Progress messages: currently not fully implemented.
##'
##' @return A list with the following components
##'
##' \code{fitted} The fitted values
##'
##' \code{V} The covariance matrix of the fitted values
##'
##' \code{residuals} The residuals
##'
##' \code{hyper} The hyperparameters (delta and nu) used in the fit. Lower values of delta give
##' less smooth regression fits, and nu is an independent measurement error variance (or "nugget").
##'
##' Prediction outside the input data points is not implemented. 
##'
##' @noRd
gp <- function(y, X, hyper="est", m=NULL, maxSample=5000,  verbose=FALSE) {
  maxSample <- min(maxSample, length(y)) # to avoid trying to invert huge matrix
  X <- as.matrix(X)[1:maxSample, , drop=FALSE]
  y <- y[1:maxSample]
  colmin <- apply(X, 2, min)
  colmax <- apply(X, 2, max)
  colrange <- colmax - colmin
  X <- sweep(X, 2, colmin, "-")
  X <- sweep(X, 2, colrange, "/")
  N <- nrow(X)
  p <- ncol(X)
  H <- cbind(1, X)
  q <- ncol(H)

    if(is.null(m)){
        m <- min(30 * p, 250)
        m <- min(length(y), m)
    }
    setForHyperparamEst <- 1:m # sample(1:N, m, replace=FALSE)
  if (identical(hyper,"est"))
    hyper <- estimate.hyperparameters(y[setForHyperparamEst], 
                                              X[setForHyperparamEst, ], verbose=verbose)
  else {
    if (!is.numeric(hyper) || (length(hyper !=2)))
      stop("`hyper` should be a numeric vector of length two, or \"est\"")
  }

    delta.hat <- hyper[1:p]
    nu.hat <- hyper[p+1]
    A <- makeA.Gaussian(X, delta.hat)
    Astar <- A + nu.hat * diag(N)
    Astarinv <- chol2inv(chol(Astar))
    rm(Astar); gc()
    AstarinvY <- Astarinv %*% y
    tHAstarinv <- t(H) %*% Astarinv
    tHAHinv <- solve(tHAstarinv %*% H + 1e-7* diag(q))
    betahat <- tHAHinv %*% (tHAstarinv %*% y)
    Hbetahat <- H %*% betahat
    resid <- y - Hbetahat
    g.hat <- Hbetahat+A %*% (Astarinv %*% resid)
    AAstarinvH <- A %*% t(tHAstarinv)
    sigmasqhat <- as.numeric(t(resid) %*% Astarinv %*% resid)/(N - q - 2)
    V <- sigmasqhat*(nu.hat * diag(N) - nu.hat ^ 2 * Astarinv +
                            (H - AAstarinvH) %*% (tHAHinv %*% t(H - AAstarinvH)))
    rm(A, Astarinv, AstarinvY, tHAstarinv, tHAHinv, betahat, Hbetahat, resid, sigmasqhat);gc()

  fitted <- unlist(g.hat)
  list(y=y, fitted=fitted, V=V, residuals=y-fitted, hyper=hyper)
}


gp.check <- function(mod){
    ## qqplot? 
    ## histogram of residuals
    graphics::hist(mod$residuals, main="Histogram of residuals", xlab="Residuals")
    ## residuals vs fitted values
    plot(mod$fitted, mod$residuals, xlab="Fitted values", ylab="Residuals")
    ## response vs fitted values 
    plot(mod$fitted, mod$y, xlab="Fitted values", ylab="Response")
}

check_plot_gp <- function(mod){
    graphics::par(mfrow=c(2,2))
    gp.check(mod)
}

check_stats_gp <- function(mod){
    invisible()
}
