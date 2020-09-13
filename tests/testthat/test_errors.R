context("Error handling")

test_that("Basic inputs are of the right format",{
    pars <- "effect"

    nb <- "wrong"
    inputs <- data.frame(baseline=1:11, effect=1:11)
    expect_error(evppi(nb, inputs, pars), "matrix, data frame or list")

    nb <- data.frame(treatment=1:10, control=1:10)
    inputs <- "wrong"
    expect_error(evppi(nb, inputs, pars), "vector, matrix or data frame")

    nb <- data.frame(treatment=1:10, control=1:10)
    inputs <- data.frame(baseline=1:11, effect=1:11)
    expect_error(evppi(nb, inputs, pars=3), "should equal the number of rows")

    nb <- data.frame(treatment=1:10, control=1:10)
    inputs <- data.frame(baseline=1:10, effect=1:10)
    expect_error(evppi(nb, inputs, pars=3), "should be a character vector")

    expect_error(evppi(nb, inputs, pars="badpar"), "not found in columns")

})


test_that("Errors in data generating function", { 
    expect_error(evsi(chemo_nb, chemo_pars, pars=pars, nsim=1000, datagen_fn="foo", verbose=FALSE),
             "`datagen_fn` should be a function")
    expect_error(evsi(chemo_nb, chemo_pars, pars=pars, nsim=1000, verbose=FALSE),
                 "`datagen_fn` should be supplied if `study` is not supplied")

    datagen_wrong <- function(inputs, n){}
    expect_error(evsi(chemo_nb, chemo_pars, pars=pars, datagen_fn=datagen_wrong, nsim=1000, verbose=FALSE),
                 "do not have default values")
        
    datagen_wrong <- function(inputs, n=100){
        array(dim=dim(inputs))
    }
    expect_error(evsi(chemo_nb, chemo_pars, pars=pars, datagen_fn=datagen_wrong, nsim=1000, verbose=FALSE),
                 "should return a data frame")

    datagen_wrong <- function(inputs, n=100){
        ret <- data.frame(rnorm(nrow(inputs)))
        names(ret) <- names(inputs)[1]
        ret
    }
    expect_error(evsi(chemo_nb, chemo_pars, pars=pars, datagen_fn=datagen_wrong, nsim=1000, verbose=FALSE),
                 "returns variables with the same names as parameters")

    datagen_wrong <- function(inputs, n=100){
        data.frame(X = rnorm(nrow(inputs) + 1))
    }
    expect_error(evsi(chemo_nb, chemo_pars, pars=pars, datagen_fn=datagen_wrong, nsim=1000, verbose=FALSE),
                 "same number of rows as `inputs`")

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
expect_error(evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
                  datagen_fn=example_datagen_fn, likelihood="foo", verbose=FALSE),
             "`likelihood` should be a function")
expect_error(evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
                  datagen_fn=example_datagen_fn, verbose=FALSE),
             "`likelihood` should be supplied")
})

test_that("Errors in likelihood for importance sampling method", { 
    lik_wrong <- function(Y, inputs, n=100, pars){ "foo" }
    expect_error(evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
                      datagen_fn=example_datagen_fn, likelihood=lik_wrong, verbose=FALSE),
                 "likelihood function should return a numeric vector")
    
    lik_wrong <- function(Y, inputs, n=100, pars){ matrix(1:4, nrow=2) }
    expect_error(evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
                      datagen_fn=example_datagen_fn, likelihood=lik_wrong, verbose=FALSE),
                 "likelihood function should return a numeric vector")

    lik_wrong <- function(Y, inputs, n=100, pars){
        rep(1, nrow(inputs) + 1)
    }
    expect_error(evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
                      datagen_fn=example_datagen_fn, likelihood=lik_wrong, verbose=FALSE),
                 "likelihood function returns a vector of length")

})
