inputs <- chemo_pars

test_that("Errors in data generating function are handled",{
    expect_error(evsi(chemo_nb, chemo_pars, datagen_fn="foo"), "should be a function")
    
    datagen_fn <- function(inputs, n){
        nsim <- nrow(inputs)
        data.frame(x1 = rbinom(nsim, size=n, prob=inputs[,"pi1"]))
    }
    expect_error(evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn), "do not have default")

    datagen_fn <- function(inputs, n=100){1}
    expect_error(evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn), "should return a data frame")

    datagen_fn <- function(inputs, n=100){
        nsim <- nrow(inputs)
        data.frame(pi1 = rbinom(nsim, size=n, prob=inputs[,"pi1"]),
                   x2 = rbinom(nsim, size=n, prob=inputs[,"pi2"])
                   )
    }
    expect_error(evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn), "same names as parameters")
    
    datagen_fn <- function(inputs, n=100){
        nsim <- nrow(inputs)
        data.frame(x1 = rbinom(nsim-1, size=n, prob=inputs[,"pi1"]))
    }
    expect_error(evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn), "same number of rows as `inputs`")

})


