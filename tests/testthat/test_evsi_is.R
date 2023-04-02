pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"
rho <- "logor_side_effects"
pars <- c(pi1,pi2)
gamf <- paste0("s(",pars,")", collapse=" + ") 

chemo_datagen_fn <- function(inputs, n=150){
  nsim <- nrow(inputs)
  with(inputs, { 
    X.SE1 <- rbinom(nsim, size=n, prob=p_side_effects_t1)
    X.SE2 <- rbinom(nsim, size=n, prob=p_side_effects_t2)
    data.frame(X.SE1, X.SE2)
  })
}

lik_chemo <- function(Y, inputs, n=150){
  loglik <-
    dbinom(Y[,"X.SE1"], size=n, inputs[,pi1], log=TRUE) +
    dbinom(Y[,"X.SE2"], size=n, inputs[,pi2], log=TRUE) 
  exp(loglik)
}

test_that("EVSI importance sampling method: user vs built in likelihood", {
  
  ## Note these are not identical because of RNG seeding.
  ## With a user-defined datagen_fn, the function is evaluated to validate it,
  ## which resets the seed.
  set.seed(1)
  expect_equal(
    evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
         datagen_fn=chemo_datagen_fn, likelihood=lik_chemo, n=100000, 
         gam_formula=gamf, verbose=FALSE)$evsi,
    322.7, tol=1e-01)
  
  set.seed(1)
  expect_equal(
    evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
         study = "trial_binary", n=100000, 
         gam_formula=gamf, verbose=FALSE)$evsi,
    323.1, tol=1e-01)
})


test_that("EVSI importance sampling method: likelihood without n argument", {
  chemo_datagen_fn <- function(inputs){
    nsim <- nrow(inputs)
    with(inputs, { 
      X.SE1 <- rbinom(nsim, size=100000, prob=p_side_effects_t1)
      X.SE2 <- rbinom(nsim, size=100000, prob=p_side_effects_t2)
      data.frame(X.SE1, X.SE2)
    })
  }
  
  datagen_fn_non <- function(inputs){
    nsim <- nrow(inputs)
    with(inputs, { 
      X.SE1 <- rbinom(nsim, size=100000, prob=p_side_effects_t1)
      X.SE2 <- rbinom(nsim, size=100000, prob=p_side_effects_t2)
      data.frame(X.SE1, X.SE2)
    })
  }   ## No check is made that the n supplied here is the same! 
  lik_chemo_non <- function(Y, inputs){
    loglik <-
      dbinom(Y[,"X.SE1"], size=100000, inputs[,pi1], log=TRUE) +
      dbinom(Y[,"X.SE2"], size=100000, inputs[,pi2], log=TRUE) 
    exp(loglik)
  }
  set.seed(1)
  expect_equal(
    evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
         datagen_fn=datagen_fn_non, likelihood=lik_chemo_non, 
         gam_formula=gamf, verbose=FALSE)$evsi,
    323.0, tol=1e-01)
})

test_that("EVSI importance sampling method: likelihood without n argument", {
  expect_warning(
    evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
         datagen_fn=chemo_datagen_fn, likelihood=lik_chemo, 
         study="trial_binary",
         gam_formula=gamf, verbose=FALSE),
    "Ignoring `likelihood`")
  expect_warning(
    evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
         datagen_fn=chemo_datagen_fn, 
         study="trial_binary",
         gam_formula=gamf, verbose=FALSE),
    "Ignoring `datagen_fn`")
})


test_that("EVSI importance sampling method: errors", {
  lik_chemo_wrong <- function(data, inputs, n=150){
    loglik <-
      dbinom(data[,"X.SE1"], size=n, inputs[,pi1], log=TRUE) +
      dbinom(data[,"X.SE2"], size=n, inputs[,pi2], log=TRUE) 
    exp(loglik)
  }
  expect_error(
    evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000,
         datagen_fn=chemo_datagen_fn, likelihood=lik_chemo_wrong),
    "should be named `Y` and `inputs`")
  lik_chemo_wrong <- function(Y){
    loglik <- 1
  }
  expect_error(
    evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000,
         datagen_fn=chemo_datagen_fn, likelihood=lik_chemo_wrong),
    "at least two arguments")
})
