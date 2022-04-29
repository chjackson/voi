## Gaussian process method for estimating EVPPI
## Code from SAVI 
## [Aren't there also packaged GP regression methods ?  gaupro, GPfit? are they more efficient?]

fitted_gp <- function(y, inputs, pars, verbose=FALSE, ...){
    args <- list(...)
    res <- gpFunc(NB=y, sets=pars, s=1000, input.parameters=inputs, m=args$gp_hyper_n, maxSample=args$maxSample, session=NULL, verbose=verbose)$fitted
    attr(res, "model") <- data.frame(y=y, fitted=res, residuals=y-res)
    res 
}

## TODO a diagnostic check method that is helpful for voi and works for any number of pars 
## residuals 

gp.check <- function(mod){
    ## qqplot? 
    ## histogram of residuals
    hist(mod$residuals, main="Histogram of residuals", xlab="Residuals")
    ## residuals vs fitted values
    plot(mod$fitted, mod$residuals, xlab="Fitted values", ylab="Residuals")
    ## response vs fitted values 
    plot(mod$fitted, mod$y, xlab="Fitted values", ylab="Response")
}

## Code below taken from SAVI, copyright (c) 2014, 2015 the SAVI authors
## https://github.com/Sheffield-Accelerated-VoI/SAVI-package/blob/master/inst/SAVI/scripts_GPfunctions.R

## Changes made by CJ for voi package
## 
## * Shiny progress messages commented out 
## * Edited to only accept a vector of INBs, since we iterate over multiple INBs outside this file.
## * gpFunc returns the fitted values and variances, since EVPPI is computed from these outside this file
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


estimate.hyperparameters <- function(NB, inputs, session, verbose=FALSE) {
  
  p <- NCOL(inputs)
  
  hyperparameters <- NA

#  progress1 <- shiny::Progress$new(session, min=1, max=D)
#  on.exit(progress1$close())
#  progress1$set(message = 'Estimating GP hyperparameters',
#                detail = 'Please wait...')
#  progress1$set(value = 1)
  

#    progress1$set(value = d)
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



gpFunc <- function(NB, sets, s=1000, input.parameters, m=NULL, maxSample=5000,  session, verbose=FALSE) {
  
#  input.parameters <- cache$params
  paramSet <- cbind(input.parameters[, sets])
  constantParams <- (apply(paramSet, 2, var) == 0)

  #remove constants
  if (sum(constantParams) == length(sets)) return(list(EVPI=0, SE=0)) # if all regressors are constant
  if (sum(constantParams) > 0) sets <- sets[-which(constantParams)] # remove constants
  
  # check for linear dependence and remove 
  paramSet <- cbind(cbind(input.parameters)[, sets]) # now with constants removed
  rankifremoved <- sapply(1:NCOL(paramSet), function (x) qr(paramSet[, -x])$rank)
  while(length(unique(rankifremoved)) > 1) {
    linearCombs <- which(rankifremoved == max(rankifremoved))
    # print(linearCombs)
    if (verbose) 
        print(paste("Linear dependence: removing column", colnames(paramSet)[max(linearCombs)]))
    paramSet <- cbind(paramSet[, -max(linearCombs)])
    sets <- sets[-max(linearCombs)]
    rankifremoved <- sapply(1:NCOL(paramSet), function(x) qr(paramSet[, -x])$rank)
  }  
  if(qr(paramSet)$rank == rankifremoved[1]) {
    paramSet <- cbind(paramSet[, -1]) # special case only lincomb left
    sets <- sets[-1]
    if (verbose) 
        print(paste("Linear dependence: removing column", colnames(paramSet)[1]))
  }
  
  inputs.of.interest <- sets
  p <- length(inputs.of.interest)
  
  maxSample <- min(maxSample, length(NB)) # to avoid trying to invert huge matrix
  
  input.matrix <- as.matrix(input.parameters[1:maxSample, inputs.of.interest, drop=FALSE])
    NB <- NB[1:maxSample]
  colmin <- apply(input.matrix, 2, min)
  colmax <- apply(input.matrix, 2, max)
  colrange <- colmax - colmin
  input.matrix <- sweep(input.matrix, 2, colmin, "-")
  input.matrix <- sweep(input.matrix, 2, colrange, "/")
  N <- nrow(input.matrix)
  p <- ncol(input.matrix)
  H <- cbind(1, input.matrix)
  q <- ncol(H)

    if(is.null(m)){
        m <- min(30 * p, 250)
        m <- min(length(NB), m)
    }
    setForHyperparamEst <- 1:m # sample(1:N, m, replace=FALSE)
    hyperparameters <- estimate.hyperparameters(NB[setForHyperparamEst], 
                                              input.matrix[setForHyperparamEst, ], session, verbose=verbose)
    
#  progress1 <- shiny::Progress$new(session, min=1, max=D)
  #on.exit(progress1$close())
#  progress1$set(message = 'Calculating conditional expected net benefits',
#                detail = 'Please wait...')
#  progress1$set(value = 1)
  
#    progress1$set(value = d)
#    print(paste("estimating g.hat for incremental NB for option", d, "versus 1"))

    delta.hat <- hyperparameters[1:p]
    nu.hat <- hyperparameters[p+1]
    A <- makeA.Gaussian(input.matrix, delta.hat)
    Astar <- A + nu.hat * diag(N)
    Astarinv <- chol2inv(chol(Astar))
    rm(Astar); gc()
    AstarinvY <- Astarinv %*% NB
    tHAstarinv <- t(H) %*% Astarinv
    tHAHinv <- solve(tHAstarinv %*% H + 1e-7* diag(q))
    betahat <- tHAHinv %*% (tHAstarinv %*% NB)
    Hbetahat <- H %*% betahat
    resid <- NB - Hbetahat
    g.hat <- Hbetahat+A %*% (Astarinv %*% resid)
    AAstarinvH <- A %*% t(tHAstarinv)
    sigmasqhat <- as.numeric(t(resid) %*% Astarinv %*% resid)/(N - q - 2)
    V <- sigmasqhat*(nu.hat * diag(N) - nu.hat ^ 2 * Astarinv +
                            (H - AAstarinvH) %*% (tHAHinv %*% t(H - AAstarinvH)))
    rm(A, Astarinv, AstarinvY, tHAstarinv, tHAHinv, betahat, Hbetahat, resid, sigmasqhat);gc()

#  progress1$close()
#  perfect.info <- mean(do.call(pmax, g.hat)) 
#  baseline <- max(unlist(lapply(g.hat, mean)))
  
 # partial.evpi <- perfect.info - baseline
  
#  print("Computing standard error via Monte Carlo")
#  tilde.g <- matrix(0, nrow=s, ncol=N)     
  
#  progress2 <- shiny::Progress$new(session, min=1, max=D)
  #on.exit(progress2$close())
#  progress2$set(message = 'Calculating Standard Error',
#                detail = 'Please wait...')
#  progress2$set(value = 1)

#    progress2$set(value = d)
#    tilde.g <- MASS::mvrnorm(s, g.hat[[d]][1:(min(N, 1000))], V[[d]][1:(min(N, 1000)), 1:(min(N, 1000))])
#  progress2$close()
  
#  sampled.perfect.info <- rowMeans(do.call(pmax, tilde.g))
#  sampled.baseline <- do.call(pmax, lapply(tilde.g, rowMeans)) 
#  rm(tilde.g);gc()
#  sampled.partial.evpi <- sampled.perfect.info - sampled.baseline
#  SE <- sd(sampled.partial.evpi)
  # g.hat.short <- lapply(g.hat,function(x) x[1:(min(N, 1000))])
  # mean.partial.evpi <- mean(do.call(pmax, g.hat.short)) - max(unlist(lapply(g.hat.short,mean)))
  # upward.bias <- mean(sampled.partial.evpi) - mean.partial.evpi 
                                        #   return(list(EVPI=partial.evpi, SE=SE))
    list(fitted=unlist(g.hat), V=V)
}

