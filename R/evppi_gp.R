### TODO acknowledge source of code.  Mark Strong? SAVI?

### TODO different code is in SAVI.  check difference.
## TODO standard errors 

## isn't there also packaged GP regression methods ?  gaupro, GPfit? are they more efficient?


fitted_gp <- function(y, inputs, poi, ...){
    fit.gp(poi, inputs, y, n.sim=length(y))$fitted
}

###GP Fitting
post.density <- function(hyperparams, parameter, x, input.matrix) {
    dinvgamma <- function(x, alpha, beta) {
        (beta^alpha)/gamma(alpha) * x^(-alpha - 1) *
          exp(-beta/x)
    }
    N <- length(x)
    p <- length(parameter)
    H <- cbind(1, input.matrix)
    q <- ncol(H)
    a.sigma <- 0.001
    b.sigma <- 0.001
    a.nu <- 0.001
    b.nu <- 1
    delta <- exp(hyperparams)[1:p]
    nu <- exp(hyperparams)[p + 1]
    A <- exp(-(as.matrix(dist(t(t(input.matrix)/delta),
                              upper = TRUE, diag = TRUE))^2))
    Astar <- A + nu * diag(N)
    T <- chol(Astar)
    y <- backsolve(t(T),(x), upper.tri = FALSE)
    x. <- backsolve(t(T), H, upper.tri = FALSE)
    tHAstarinvH <- t(x.) %*% (x.)
    betahat <- solve(tHAstarinvH) %*% t(x.) %*% y
    residSS <- y %*% y - t(y) %*% x. %*% betahat - t(betahat) %*%
      t(x.) %*% y + t(betahat) %*% tHAstarinvH %*% betahat
    prior <- prod(dnorm(log(delta), 0, sqrt(1e+05))) *
      dinvgamma(nu, a.nu, b.nu)
    l <- -sum(log(diag(T))) - 1/2 * log(det(tHAstarinvH)) -
      (N - q + 2 * a.sigma)/2 * log(residSS/2 + b.sigma) +
      log(prior)
    names(l) <- NULL
    return(l)
}


estimate.hyperparameters <- function(x, input.matrix, parameter,n.sim) {
    p <- length(parameter)
    initial.values <- rep(0, p + 1)
    repeat {
        log.hyperparameters <- optim(initial.values,
                                     fn = post.density,parameter=parameter, x = x[1:n.sim],
                                     input.matrix = input.matrix[1:n.sim, ],
                                     method = "Nelder-Mead", control = list(fnscale = -1,
                                                                            maxit = 10000, trace = 0))$par
        if (sum(abs(initial.values - log.hyperparameters)) <
            0.01) {
            hyperparameters <- exp(log.hyperparameters)
            break
        }
        initial.values <- log.hyperparameters
    }
    return(hyperparameters)
}



fit.gp <- function(parameter, inputs, x, n.sim) {
    tic <- proc.time()
    p <- length(parameter)
    input.matrix <- as.matrix(inputs[, parameter, drop = FALSE])
    colmin <- apply(input.matrix, 2, min)
    colmax <- apply(input.matrix, 2, max)
    colrange <- colmax - colmin
    input.matrix <- sweep(input.matrix, 2, colmin, "-")
    input.matrix <- sweep(input.matrix, 2, colrange,
                          "/")
    N <- nrow(input.matrix)
    H <- cbind(1, input.matrix)
    q <- ncol(H)
    hyperparameters <- estimate.hyperparameters(x = x,input.matrix = input.matrix, parameter = parameter, n.sim = n.sim)
    delta.hat <- hyperparameters[1:p]
    nu.hat <- hyperparameters[p + 1]
    A <- exp(-(as.matrix(dist(t(t(input.matrix)/delta.hat),
                              upper = TRUE, diag = TRUE))^2))
    Astar <- A + nu.hat * diag(N)
    Astarinv <- chol2inv(chol(Astar))
    rm(Astar)
    gc()
    AstarinvY <- Astarinv %*% x
    tHAstarinv <- t(H) %*% Astarinv
    tHAHinv <- solve(tHAstarinv %*% H)
    betahat <- tHAHinv %*% (tHAstarinv %*% x)
    Hbetahat <- H %*% betahat
    resid <- x - Hbetahat
    fitted<- Hbetahat + A %*% (Astarinv %*%
                               resid)
    AAstarinvH <- A %*% t(tHAstarinv)
    sigmasqhat <- as.numeric(t(resid) %*% Astarinv %*%
                             resid)/(N - q - 2)
    rm(A, Astarinv, AstarinvY, tHAstarinv, tHAHinv,
       Hbetahat, resid, sigmasqhat)
    gc()
    toc <- proc.time() - tic
    time <- toc[3]
    names(time) = "Time to fit GP regression (seconds)"
    list(fitted = fitted,time = time, fit=NULL,formula = NULL)
}
