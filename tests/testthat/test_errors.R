context("Error handling")

test_that("Basic inputs are of the right format",{
    poi <- "effect"

    nb <- "wrong"
    inputs <- data.frame(baseline=1:11, effect=1:11)
    expect_error(evppi(nb, inputs, poi), "should be a matrix, data frame")

    nb <- data.frame(treatment=1:10, control=1:10)
    inputs <- "wrong"
    expect_error(evppi(nb, inputs, poi), "should be a matrix or data frame")

    nb <- data.frame(treatment=1:10, control=1:10)
    inputs <- data.frame(baseline=1:11, effect=1:11)
    expect_error(evppi(nb, inputs, poi=3), "should equal the number of rows")

    nb <- data.frame(treatment=1:10, control=1:10)
    inputs <- data.frame(baseline=1:10, effect=1:10)
    expect_error(evppi(nb, inputs, poi=3), "should be a character vector")

    expect_error(evppi(nb, inputs, poi="badpar"), "not found in columns")

})


test_that("Errors in data generating function", { 
    expect_error(evsi(chemo_nb, chemo_pars, poi=poi, nsim=1000, datagen_fn="foo", verbose=FALSE),
             "`datagen_fn` should be a function")
    expect_error(evsi(chemo_nb, chemo_pars, poi=poi, nsim=1000, verbose=FALSE),
                 "`datagen_fn` should be supplied if `study` is not supplied")

    ## TODO no defaults in sample size args
    ## TODO doesn't return a data frame 
    ## TODO returns variables with same names as parameters 
    ## TODO returns wrong number of rows
})

example_datagen_fn <- function(inputs, n=150){
    nsim <- nrow(inputs)
    with(inputs, { 
        X.SE1 <- rbinom(nsim, size=n, prob=pi1)
        X.SE2 <- rbinom(nsim, size=n, prob=pi2)
        data.frame(X.SE1, X.SE2)
    })
}

test_that("Errors for importance sampling method",{
expect_error(evsi(chemo_nb, chemo_pars, method="is", poi=poi, nsim=1000, 
                  datagen_fn=example_datagen_fn, likelihood="foo", verbose=FALSE),
             "`likelihood` should be a function")
expect_error(evsi(chemo_nb, chemo_pars, method="is", poi=poi, nsim=1000, 
                  datagen_fn=example_datagen_fn, verbose=FALSE),
             "`likelihood` should be supplied")
})


