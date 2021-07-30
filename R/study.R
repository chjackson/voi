## BUILT-IN STUDY DESIGNS

## Single-arm study of a binary outcome

datagen_binary <- function(inputs, n=100, pars){
    nsim <- nrow(inputs)
    data.frame(
        X1 = rbinom(nsim, size=n, prob=inputs[,pars[1]])
    )
}

likelihood_binary <- function(Y, inputs, n=100, pars){
    loglik <-
        dbinom(Y[,"X1"], size=n, inputs[,pars[1]], log=TRUE) 
    exp(loglik)
}


## Two-arm trial of a binary outcome

datagen_trial_binary <- function(inputs, n=100, pars){
    nsim <- nrow(inputs)
    data.frame(
        X1 = rbinom(nsim, size=n, prob=inputs[,pars[1]]),
        X2 = rbinom(nsim, size=n, prob=inputs[,pars[2]])
    )
}

likelihood_trial_binary <- function(Y, inputs, n=100, pars){
    loglik <-
        dbinom(Y[,"X1"], size=n, inputs[,pars[1]], log=TRUE) + 
        dbinom(Y[,"X2"], size=n, inputs[,pars[2]], log=TRUE) 
    exp(loglik)
}


## Single-arm study of a normal outcome with known variance (supplied as an argument) 
## Return an estimate of the mean from a study of size n 

datagen_normal_known <- function(inputs, n=100, pars, sd=1){
    nsim <- nrow(inputs)
    mu <- inputs[,pars[1]]
    data.frame(
        X1 = rnorm(nsim, mu, sd/sqrt(n))
    )
}

likelihood_normal_known <- function(Y, inputs, n=100, pars, sd=1){
    mu <- inputs[,pars[1]]
    dnorm(Y[,"X1"], mu, sd/sqrt(n))
}


studies_builtin <- c("binary","trial_binary","normal_known")
