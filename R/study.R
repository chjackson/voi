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

studies_builtin <- c("binary","trial_binary")
