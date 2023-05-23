## Toy model : INB(p1,p2) = p1+p2, p1~N(0,1), p2~N(0,2)

model_fn <- function(p1, p2) {c(0, p1 + p2)}
s1 <- 1
s2 <- 2
par_fn <- function(n){data.frame(p1 = rnorm(n, 0, s1), p2 = rnorm(n, 0, s2))}

test_that("Monte Carlo EVPPI, net benefit format",{
  set.seed(1)
  expect_equal(evppi_mc(model_fn, par_fn, pars="p1", nouter=100, ninner=50)$evppi, 0.3546835, tol=1e-01)
  expect_equal(evppi_mc(model_fn, par_fn, pars="p2", nouter=100, ninner=50)$evppi, 0.6015171, tol=1e-01)
})

test_that("Errors in Monte Carlo EVPPI",{
  
  expect_error(evppi_mc(model_fn, par_fn, pars="p2", nouter=100, ninner="foo"), "should be a number")
  expect_error(evppi_mc(model_fn, par_fn, pars="p2", nouter="woo", ninner=100), "should be a number")
  expect_error(evppi_mc(model_fn, par_fn, pars="p2", nouter=-1, ninner=100), "greater than 1")
  expect_error(evppi_mc(model_fn, par_fn, pars="p3", nouter=10, ninner=100), "not found in arguments")

  bad_modelfn <- function(p1, p2) { matrix(1:9, nrow=3)}  
  expect_error(evppi_mc(bad_modelfn, par_fn, pars="p2", nouter=10, ninner=100), "should have two rows")
  bad_modelfn <- function(p1, p2) { 1 }  
  expect_error(evppi_mc(bad_modelfn, par_fn, pars="p2", nouter=10, ninner=100), "more than one decision")
  bad_modelfn <- function(p1, p2) { function()1 }  
  expect_error(evppi_mc(bad_modelfn, par_fn, pars="p2", nouter=10, ninner=100), 
               "should return a vector, matrix or data frame")
    
  bad_parfn <- function(n){1}
  expect_error(evppi_mc(model_fn, bad_parfn, pars="p2", nouter=10, ninner=100), "should return a named vector")
  bad_parfn <- function(n){rbind(rnorm(n,0,s1), rnorm(n,0,s2))}
  expect_error(evppi_mc(model_fn, bad_parfn, pars="p2", nouter=10, ninner=100), "columns should be named")
  bad_parfn <- function(n){function()1}
  expect_error(evppi_mc(model_fn, bad_parfn, pars="p2", nouter=10, ninner=100), "should return")
  bad_parfn <- function(n){cbind(p1=rnorm(n,0,s1), p3=rnorm(n,0,s2))}
  expect_error(evppi_mc(model_fn, bad_parfn, pars="p2", nouter=10, ninner=100), "parameters of `model_fn` were not found")
  bad_parfn <- function(n){stop()}
  expect_error(evppi_mc(model_fn, bad_parfn, pars="p2", nouter=10, ninner=100), "returned the following error")
  expect_error(evppi_mc(model_fn, par_fn_foo, pars="p1", nouter=100, ninner=50)$evppi, "not found")
  bad_parfn <- function(n){c(p1=1,p2=1)}
  expect_error(evppi_mc(model_fn, bad_parfn, pars="p2", nouter=10, ninner=100), 
               "for n>1 should return a matrix or data frame")
  bad_parfn <- function(n){data.frame(p1=1,p2=1)}
  expect_error(evppi_mc(model_fn, bad_parfn, pars="p2", nouter=10, ninner=100),
               "should have n rows")
  
})

## Toy CEA model with effects and no costs 
model_fn <- function(p1, p2) {rbind(e = c(0, p1 + p2), 
                                    c = c(0, 0))}
s1 <- 1
s2 <- 2
par_fn <- function(n){data.frame(p1 = rnorm(n, 0, s1), p2 = rnorm(n, 0, s2))}

test_that("Monte Carlo EVPPI, CEA format",{
  set.seed(1)
  expect_equal(evppi_mc(model_fn, par_fn, pars="p1", nouter=50, ninner=50, k=1)$evppi, 0.245, tol=1e-01)
  expect_equal(evppi_mc(model_fn, par_fn, pars="p1", nouter=50, ninner=50, k=c(1,2))$evppi, c(0.271, 0.096), tol=1e-01)
})

test_that("Errors in Monte Carlo EVPPI for CEA format",{
  expect_error(evppi_mc(model_fn, par_fn, pars="p1", nouter=100, ninner=50), "at least one willingness-to-pay")
})


test_that("Monte Carlo EVPPI: Additional arguments to model_fn",{
  model_fn <- function(p1, p2, basenb) {c(basenb, p1 + p2)}
  expect_error(evppi_mc(model_fn, par_fn, pars="p1", nouter=100, ninner=50), 
               "parameters of `model_fn` were not found")
  ## can either supply in mfargs
  expect_error(evppi_mc(model_fn, par_fn, pars="p1", nouter=10, ninner=10, mfargs=list(basenb=0)), NA)
  ## or fall back to a default value in model_fn
  model_fn <- function(p1, p2, basenb=0) {c(basenb, p1 + p2)}
  expect_error(evppi_mc(model_fn, par_fn, pars="p1", nouter=10, ninner=10), NA)
  ## note that mfargs will override any default value
  expect_error(evppi_mc(model_fn, par_fn, pars="p1", nouter=10, ninner=10, mfargs=list(basenb=0)), NA)
})


test_that("Monte Carlo EVPPI: Correlated parameters",{
  model_fn <- function(p1, p2) {c(0, p1 + p2)}
  par_fn_corr <- function(n, p1=NULL){
    p1_new <- if (is.null(p1)) rnorm(n, 1, 1) else p1
    data.frame(p1 = p1_new,
               p2 = rnorm(n, p1_new, 2))
  }
  set.seed(1)
  expect_equal(evppi_mc(model_fn, par_fn_corr, pars="p1", nouter=100, ninner=50)$evppi, 0.237, tol=1e-01)
})
