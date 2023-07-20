## BUILT-IN STUDY DESIGNS
## Not exported for user use - they will break because of the extra "pars" argument

## Single-arm study of a binary outcome

datagen_binary <- function(inputs, n=100, pars){
    nsim <- nrow(inputs)
    validate_prob(inputs[,pars[1]], pars[1])
    data.frame(
        X1 = rbinom(nsim, size=n, prob=inputs[,pars[1]])
    )
}

likelihood_binary <- function(Y, inputs, n=100, pars){
    validate_prob(inputs[,pars[1]], pars[1])
    loglik <-
        dbinom(Y[,"X1"], size=n, inputs[,pars[1]], log=TRUE) 
    exp(loglik)
}

analysis_binary <- function(data, args, pars, ...){
    if (is.null(args$niter))
      args$niter <- 1000
    check_analysis_args(args, c("a", "b"))
    res <- data.frame(rbeta(args$niter,
                            shape1 = data$X1 + args$a,
                            shape2 = args$n - data$X1 + args$b))
    names(res) <- pars
    res
}

analysis_trial_binary <- function(data, args, pars,...){
    if (is.null(args$niter))
      args$niter <- 1000
    check_analysis_args(args, c("a1", "a2", "b1", "b2"))
    p1 <- rbeta(args$niter,
                shape1 = data$X1 + args$a1,
                shape2 = args$n - data$X1 + args$b1)
    p2 <- rbeta(args$niter,
                shape1 = data$X2 + args$a2,
                shape2 = args$n - data$X2 + args$b2)
    res <- data.frame(p1, p2)
    names(res) <- pars
    res
}

## Could we do one that estimates a log odds ratio from a study with two binary outcomes???
## But if we did a study of two treatments we'd get the absolute outcomes too.

## Two-arm trial of a binary outcome

datagen_trial_binary <- function(inputs, n=100, pars){
    nsim <- nrow(inputs)
    validate_prob(inputs[,pars[1]], pars[1])
    validate_prob(inputs[,pars[2]], pars[2])
    data.frame(
        X1 = rbinom(nsim, size=n, prob=inputs[,pars[1]]),
        X2 = rbinom(nsim, size=n, prob=inputs[,pars[2]])
    )
}

likelihood_trial_binary <- function(Y, inputs, n=100, pars){
    validate_prob(inputs[,pars[1]], pars[1])
    validate_prob(inputs[,pars[2]], pars[2])
    loglik <-
        dbinom(Y[,"X1"], size=n, inputs[,pars[1]], log=TRUE) + 
        dbinom(Y[,"X2"], size=n, inputs[,pars[2]], log=TRUE) 
    exp(loglik)
}


## Single-arm study of a normal outcome with known variance (supplied as an argument) 
## Return an estimate of the mean from a study of size n 

datagen_normal_known <- function(inputs, n=100, pars, sd=1){
    data.frame(
      X1 = rnorm(nrow(inputs),
                 mean = inputs[,pars[1]],
                 sd = sd/sqrt(n))
    )
}

likelihood_normal_known <- function(Y, inputs, n=100, pars, sd=1){
    mu <- inputs[,pars[1]]
    dnorm(Y[,"X1"], mu, sd/sqrt(n))
}

analysis_normal_known <- function(data, args, pars,...){
  if (is.null(args$niter))
    args$niter <- 1000
  check_analysis_args(args, c("prior_mean", "prior_sd", "sampling_sd"))
  
  xbar_var <- args$sampling_sd^2 / args$n
  w <- args$prior_sd^2 / (args$prior_sd^2 + xbar_var)
  post_mean <- w*data$X1 + (1-w)*args$prior_mean
  post_var <- 1 / (1 / args$prior_sd^2 + 1 / xbar_var)
  res <- data.frame(rnorm(args$niter, post_mean, sqrt(post_var)))
  names(res) <- pars
  res
}


studies_builtin <- c("binary","trial_binary","normal_known")

check_analysis_args <- function(args, required){
  for (i in required){
    if (is.null(args[[i]]))
      stop(sprintf("`%s` not supplied in `analysis_args`", i))
  }
}

validate_prob <- function(p, name){
  if (any((p<0)|(p>1)))
    stop(sprintf("%s should represent a probability, but it has values that are greater than 1 or less than 0", name))
}
