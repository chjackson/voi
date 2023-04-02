pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"
rho <- "logor_side_effects"

if (0){
chemo_datagen_fn <- function(inputs, n=150){
    nsim <- nrow(inputs)
    with(inputs, { 
        X.SE1 <- rbinom(nsim, size=n, prob=p_side_effects_t1)
        X.SE2 <- rbinom(nsim, size=n, prob=p_side_effects_t2)
        data.frame(X.SE1, X.SE2)
    })
}

}

test_that("EVSI with built-in study designs", {
    set.seed(1)
    expect_equal(
        evsi(chemo_nb, chemo_pars, study="trial_binary", pars=c(pi1, pi2), verbose=FALSE)$evsi, 
        242.1268, tol=0.01)
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
    239.1277, tol=0.01)
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
  expect_equal(e4$evsi[e4$k==20000 & e4$n==100], e1$evsi, tol=1e-03)
  set.seed(1)
  e1 <- evsi(chemo_cea, chemo_pars, study="binary", n=c(100,1000), pars=c(pi2), 
             verbose=FALSE, method="is", nsim=1000) 
  expect_equal(e1$evsi[2], 226.7225, tol=1e-04)
  set.seed(1)
  e2 <- evsi(chemo_nb, chemo_pars, study="binary", n=c(100,1000), pars=c(pi2), 
             verbose=FALSE, method="is", nsim=1000) 
  expect_equal(e2$evsi[1], 219.236, tol=1e-04)
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
