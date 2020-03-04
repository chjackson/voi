## Built-in study designs 

datagen_trial_binary <- function(inputs, n=100, poi){
    nsim <- nrow(inputs)
    data.frame(
        X1 = rbinom(nsim, size=n, prob=inputs[,poi[1]]),
        X2 = rbinom(nsim, size=n, prob=inputs[,poi[2]])
    )
}

studies_builtin <- list("trial_binary" = datagen_trial_binary)
