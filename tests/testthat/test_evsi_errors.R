inputs <- chemo_pars

test_that("Errors in data generating function are handled",{
    expect_error(evsi(chemo_nb, chemo_pars,  datagen_fn="foo"), 
                 "should be a function")
    
    datagen_fn <- function(inputs, n){
        nsim <- nrow(inputs)
        data.frame(x1 = rbinom(nsim, size=n, prob=inputs[,"p_side_effects_t1"]))
    }
    expect_error(evsi(chemo_nb, chemo_pars,  datagen_fn=datagen_fn), 
                 "do not have default")

    datagen_fn <- function(inputs, n=100){1}
    expect_error(evsi(chemo_nb, chemo_pars,  datagen_fn=datagen_fn), 
                 "should return a data frame")

    datagen_fn <- function(inputs, n=100){
        nsim <- nrow(inputs)
        data.frame(p_side_effects_t1 = rbinom(nsim, size=n, prob=inputs[,"p_side_effects_t1"]),
                   x2 = rbinom(nsim, size=n, prob=inputs[,"p_side_effects_t2"])
                   )
    }
    expect_error(evsi(chemo_nb, chemo_pars,  datagen_fn=datagen_fn), 
                 "same names as parameters")
    
    datagen_fn <- function(inputs, n=100){
        nsim <- nrow(inputs)
        data.frame(x1 = rbinom(nsim-1, size=n, prob=inputs[,"p_side_effects_t1"]))
    }
    expect_error(evsi(chemo_nb, chemo_pars,  datagen_fn=datagen_fn), 
                 "same number of rows as `inputs`")

})

test_that("Errors in sample size are handled", {
    expect_error(evsi(chemo_nb, chemo_pars, study="binary", pars="p_side_effects_t1", n="foo"),
                 "should be a numeric vector")
    expect_error(evsi(chemo_nb, chemo_pars, study="binary", pars="p_side_effects_t1", n=c(10, 100, -1)),
                 "should all be positive integers")
})

test_that("Errors in analysis function are handled",{
})
