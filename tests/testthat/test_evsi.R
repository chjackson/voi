pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"
rho <- "logor_side_effects"

## TODO change the parameter names in this block. what are gamma and recover? 

if (0){
chemo_datagen_fn <- function(inputs, n=150){
    nsim <- nrow(inputs)
    with(inputs, { 
        X.SE1 <- rbinom(nsim, size=n, prob=pi1)
        X.SE2 <- rbinom(nsim, size=n, prob=pi2)
        X.N.hosp <- rbinom(nsim, X.SE1 + X.SE2, gamma.hosp)
        X.N.die <- rbinom(nsim, X.N.hosp, gamma.dead)
        N.amb <- X.SE1 + X.SE2 - X.N.hosp
        X.amb <- rgamma(nsim, N.amb, recover.amb) # sum of N.amb independent exponentials
        N.hosp <- X.N.hosp - X.N.die
        X.hosp <- rgamma(nsim, N.hosp, recover.hosp)
        data.frame(X.amb, X.SE1, X.SE2, X.N.hosp, X.hosp, X.N.die)
    })
}

test_that("EVSI GAM method", { 
    gamf <- "s(X.amb) + s(X.SE1) + s(X.SE2) + s(X.N.hosp) + s(X.hosp) + s(X.N.die)"
    ## method depends on simulating data, so make tests reproducible
    set.seed(1)
    expect_equal(evsi(chemo_nb, chemo_pars, datagen_fn=chemo_datagen_fn,
                      gam_formula=gamf, nsim=1000)$evsi,
                 18.18446, tol=0.01)
})

lik_chemo <- function(Y, inputs){
  loglik <-
      dbinom(Y[,"X.SE1"], size=150, inputs[,pi1], log=TRUE) +
      dbinom(Y[,"X.SE2"], size=150, inputs[,pi2], log=TRUE) +
      dbinom(Y[,"X.N.hosp"], Y[,"X.SE1"] + Y[,"X.SE2"], inputs[,"gamma.hosp"], log=TRUE) +
      dbinom(Y[,"X.N.die"],  Y[,"X.N.hosp"], inputs[,"gamma.dead"], log=TRUE) +
      log(inputs[,"recover.amb"]) * (Y[,"X.SE1"] + Y[,"X.SE2"] - Y[,"X.N.hosp"]) -
      inputs[,"recover.amb"] * Y[,"X.amb"] + 
      log(inputs[,"recover.hosp"]) * (Y[,"X.N.hosp"] - Y[,"X.N.die"]) -
      inputs[,"recover.hosp"] * Y[,"X.hosp"]
  exp(loglik)
}

test_that("EVSI importance sampling method", {
    pars <- c(pi1,pi2,"gamma.hosp","gamma.dead","recover.amb","recover.hosp")
    gamf <- paste0("s(",pars,")", collapse=" + ") 
    set.seed(1)
    expect_equal(
        evsi(chemo_nb, chemo_pars, method="is", pars=pars, nsim=1000, 
             datagen_fn=chemo_datagen_fn, likelihood=lik_chemo, gam_formula=gamf, verbose=FALSE)$evsi,
        9.252767, tol=0.01)
})

}

test_that("EVSI with built-in study designs", {
    set.seed(1)
    expect_equal(
        evsi(chemo_nb, chemo_pars, study="trial_binary", pars=c(pi1, pi2), verbose=FALSE)$evsi, 
        777.1786, tol=0.01)
    set.seed(1)
    e1 <- evsi(chemo_nb, chemo_pars, study="binary", n=100, pars=c(pi2), verbose=FALSE)$evsi
    e2 <- evsi(chemo_nb, chemo_pars, study="binary", n=10000, pars=c(pi2), verbose=FALSE)$evsi
    expect_gt(e2, e1)
    set.seed(1)
    e1 <- evsi(chemo_nb, chemo_pars, study="normal_known", n=100, pars=c(pi2), verbose=FALSE)$evsi
    e2 <- evsi(chemo_nb, chemo_pars, study="normal_known", n=10000, pars=c(pi2), verbose=FALSE)$evsi
    expect_gt(e2, e1)
    set.seed(1)
    e1 <- evsi(chemo_nb, chemo_pars, study="normal_known", n=10000, aux_pars=list(sd=2), pars=c(pi2), verbose=FALSE)$evsi 
    e2 <- evsi(chemo_nb, chemo_pars, study="normal_known", n=10000, aux_pars=list(sd=0.1), pars=c(pi2), verbose=FALSE)$evsi
    ep <- evppi(chemo_nb, chemo_pars, pars=c(pi2), verbose=FALSE)$evppi
    expect_gt(e2, e1)
    expect_gt(ep, e2)
})

test_that("EVSI with built-in study designs: IS method", {
    set.seed(1)
  expect_equal(
    evsi(chemo_nb, chemo_pars, study="trial_binary", pars=c(pi1, pi2), 
         method="is", nsim=1000, verbose=FALSE)$evsi
    , 
    789.1255, tol=0.01)
  
  ## can we include extra cols 
  chemo_pars$xextra <- 1 
  expect_equal(
    evsi(chemo_nb, chemo_pars, study="trial_binary", pars=c(pi1, pi2), 
         method="is", nsim=1000, verbose=FALSE)$evsi
    , 
    789.1255, tol=0.01)
})

test_that("EVSI with multiple sample sizes", { 
  set.seed(1)
  e1 <- evsi(chemo_nb, chemo_pars, study="binary", n=100, pars=c(pi1), verbose=FALSE)
  e2 <- evsi(chemo_nb, chemo_pars, study="binary", n=1000, pars=c(pi1), verbose=FALSE) 
  set.seed(1)
  e3 <- evsi(chemo_nb, chemo_pars, study="binary", n=c(100,1000), pars=c(pi1), verbose=FALSE) 
  expect_equal(e1$evsi, e3$evsi[1])
  expect_equal(e2$evsi, e3$evsi[2])
  set.seed(1)
  e4 <- evsi(chemo_cea, chemo_pars, study="binary", n=c(100,1000), pars=c(pi1), verbose=FALSE) 
  expect_equal(e4$evsi[e4$k==30000 & e4$n==100], e1$evsi, tol=1e-03)
  set.seed(1)
  e1 <- evsi(chemo_cea, chemo_pars, study="binary", n=c(100,1000), pars=c(pi2), 
             verbose=FALSE, method="is", nsim=1000) 
  expect_equal(e1$evsi[2], 617.0609, tol=1e-04)
  set.seed(1)
  e2 <- evsi(chemo_nb, chemo_pars, study="binary", n=c(100,1000), pars=c(pi2), 
             verbose=FALSE, method="is", nsim=1000) 
  expect_equal(e2$evsi[1], 571.2093, tol=1e-04)
})

test_that("EVSI with multiple sample sizes and CEA output", { 
  set.seed(1)
  e1 <- evsi(chemo_cea, chemo_pars, study="binary", n=100, pars=c(pi2), verbose=FALSE)
  e2 <- evsi(chemo_cea, chemo_pars, study="binary", n=1000, pars=c(pi2), verbose=FALSE) 
  set.seed(1)
  e12 <- evsi(chemo_cea, chemo_pars, study="binary", n=c(100,1000), pars=c(pi2), verbose=FALSE) 
  expect_equal(e12$evsi, c(e1$evsi,e2$evsi), tol=1e-03)
  expect_equal(e12$n, c(e1$n,e2$n))
})
