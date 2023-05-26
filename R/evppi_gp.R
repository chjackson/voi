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

makeA.Gaussian.old <- function(X, phi) { # function to make A matrix with the Gaussian correlation function
  n <- NROW(X)
  if(length(phi) > 1) {
    b <- diag(phi ^ (-2))
  } else {
    b <- phi ^ (-2) }
  R <- X %*% as.matrix(b) %*% t(X)
  
  S <- matrix(diag(R), nrow = n, ncol = n)
  exp(2 * R - S - t(S))
}

##' Form a squared exponential correlation matrix between two sets of predictors for a Gaussian process regression
##'
##' @param X First set of predictors.  Matrix with number of columns given by the number of
##' predictors in the model, and number of rows given by the number of alternative values for these. 
##'
##' @param Xstar Second set of predictors, in the same format. 
##'
##' @param phi correlation parameter. Vector with length given by the number of predictors.
##'
##' @return A matrix with \code{nrow(X)} and \code{nrow(Xstar)} columns, with \eqn{r},\eqn{s} entry
##' given by the correlation between the
##' \eqn{r}th element of \code{X} and the \eqn{s}th element of \code{Xstar}.
##'
##' @noRd
makeA.Gaussian <- function(X, Xstar, phi){
  if (!is.matrix(X)) stop("`X` should be a matrix")
  if (!is.matrix(Xstar)) stop("`Xstar` should be a matrix")
  n <- nrow(X)
  m <- nrow(Xstar)
  if (!(n>0)) stop("X should have 1 or more rows")
  if (!(m>0)) stop("Xstar should have 1 or more rows")
  X_rep <- X[rep(1:n, m), , drop=FALSE]
  Xstar_rep <- Xstar[rep(1:m, each=n), , drop=FALSE]
  scale_rep <- matrix(phi, nrow=n*m, ncol=length(phi), byrow=TRUE)
  dists <- rowSums(((X_rep - Xstar_rep) / scale_rep)^2)
  matrix(exp(-dists), nrow=n, ncol=m)
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
  
  A <- makeA.Gaussian(input.m, input.m, delta)
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
##' Fit a Gaussian process regression.  This is a simple
##' implementation, in pure R, that is designed to be sufficient for
##' doing VoI calculations, but not polished and tuned for any other
##' purpose.
##'
##' @param y Vector of outcome data
##' 
##' @param X Matrix of inputs
##'
##' @param Xpred Matrix of inputs at which predictions are wanted
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
##' @param verbose Progress messages (not thoroughly implemented).
##'
##' @return A list with the following components
##'
##' \code{fitted} The fitted values at the input data points 
##'
##' \code{pred} The fitted values at \code{Xpred}

##' \code{V} The covariance matrix of the fitted values 
##'
##' \code{residuals} The residuals
##'
##' \code{hyper} The hyperparameters (delta and nu) used in the fit. Lower values of delta give
##' less smooth regression fits, and nu is an independent measurement error variance (or "nugget").
##'
##' The variance of the predicted points is not implemented.
##'
##' @noRd
gp <- function(y, X, Xpred=NULL, hyper="est", m=NULL, maxSample=5000,  verbose=FALSE) {
  maxSample <- min(maxSample, length(y)) # to avoid trying to invert huge matrix
  X <- as.matrix(X)[1:maxSample, , drop=FALSE]
  y <- y[1:maxSample]

  standardiseX <- function(X){
    colmin <- apply(X, 2, min)
    colmax <- apply(X, 2, max)
    colrange <- colmax - colmin
    X <- sweep(X, 2, colmin, "-")
    X <- sweep(X, 2, colrange, "/")
    X
  }
  X <- standardiseX(X)

  p <- ncol(X)
  if (is.null(Xpred))
    Xpred <- X
  else {
    Xpred <- as.matrix(Xpred)
    Xpred <- standardiseX(Xpred)
  }  
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
    if (!is.numeric(hyper) || (length(hyper) !=2))
      stop("`hyper` should be a numeric vector of length two, or \"est\"")
  }

    delta.hat <- hyper[1:p]
    nu.hat <- hyper[p+1]
    A <- makeA.Gaussian(X, X, delta.hat)
    Apred <- makeA.Gaussian(Xpred, X, delta.hat)
    N <- nrow(X)
    Astar <- A + nu.hat * diag(N)
    Astarinv <- chol2inv(chol(Astar))
    rm(Astar); gc()
    AstarinvY <- Astarinv %*% y
    tHAstarinv <- t(H) %*% Astarinv
    tHAHinv <- solve(tHAstarinv %*% H + 1e-7* diag(q))
    betahat <- tHAHinv %*% (tHAstarinv %*% y)
    Hbetahat <- H %*% betahat
    resid <- y - Hbetahat

    g.hat <- Hbetahat  +  A %*% (Astarinv %*% resid)

    Hpred <- cbind(1, Xpred)
    Hbetahatpred <- Hpred %*% betahat
    pred <- Hbetahatpred  +  Apred %*% (Astarinv %*% resid)

    AAstarinvH <- A %*% t(tHAstarinv)

    ## get the variance V of the fitted values 
    sigmasqhat <- as.numeric(t(resid) %*% Astarinv %*% resid)/(N - q - 2)
    V <- sigmasqhat*(nu.hat * diag(N) -
                     nu.hat ^ 2 * Astarinv +
                     (H - AAstarinvH) %*% (tHAHinv %*% t(H - AAstarinvH)))

    rm(A, Astarinv, AstarinvY, tHAstarinv, tHAHinv, Hbetahat, resid);gc()

  fitted <- unlist(g.hat)
  pred <- unlist(pred)
  list(y=y, fitted=fitted, V=V, pred=pred, residuals = y - fitted, hyper=hyper, beta=betahat, sigmasq=sigmasqhat)
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
  oldpar <- graphics::par(no.readonly=TRUE)
  on.exit(par(oldpar))
  graphics::par(mfrow=c(2,2))
  gp.check(mod)
}

check_stats_gp <- function(mod){
    invisible()
}
