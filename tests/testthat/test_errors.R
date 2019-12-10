context("Error handling")

test_that("Basic inputs are of the right format",{
    poi <- "effect"

    nb <- "wrong"
    pars <- data.frame(baseline=1:11, effect=1:11)
    expect_error(evppi(nb, pars, poi), "should be a matrix or data frame")

    nb <- data.frame(treatment=1:10, control=1:10)
    pars <- "wrong"
    expect_error(evppi(nb, pars, poi), "should be a matrix or data frame")

    nb <- data.frame(treatment=1:10, control=1:10)
    pars <- data.frame(baseline=1:11, effect=1:11)
    expect_error(evppi(nb, pars, poi=3), "should equal the number of rows")

    nb <- data.frame(treatment=1:10, control=1:10)
    pars <- data.frame(baseline=1:10, effect=1:10)
    expect_error(evppi(nb, pars, poi=3), "should be a character vector")

    expect_error(evppi(nb, pars, poi="badpar"), "not found in columns")

})



