## Built-in study designs 

datagen_binary <- function(inputs, n=100, poi){
    nsim <- nrow(inputs)
    data.frame(
        X1 = rbinom(nsim, size=n, prob=inputs[,poi[1]])
    )
}

likelihood_binary <- function(Y, inputs, n=100, poi){
    loglik <-
        dbinom(Y[,"X1"], size=n, inputs[,poi[1]], log=TRUE) 
    exp(loglik)
}

datagen_trial_binary <- function(inputs, n=100, poi){
    nsim <- nrow(inputs)
    data.frame(
        X1 = rbinom(nsim, size=n, prob=inputs[,poi[1]]),
        X2 = rbinom(nsim, size=n, prob=inputs[,poi[2]])
    )
}

likelihood_trial_binary <- function(Y, inputs, n=100, poi){
    loglik <-
        dbinom(Y[,"X1"], size=n, inputs[,poi[1]], log=TRUE) + 
        dbinom(Y[,"X2"], size=n, inputs[,poi[2]], log=TRUE) 
    exp(loglik)
}

studies_builtin <- c("binary","trial_binary")
