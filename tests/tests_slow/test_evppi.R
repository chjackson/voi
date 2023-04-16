pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"
rho <- "logor_side_effects"

if (requireNamespace(INLA,quietly=TRUE)) { 
  test_that("EVPPI with INLA",{
    expect_error(evppi(chemo_nb, chemo_pars, pars=pi2, method="inla", nsim=100), "2 or more parameters")
    pars <- c(pi2,rho)
    set.seed(1)
    expect_equal(
      evppi(chemo_nb, chemo_pars, pars=pars, method="inla", nsim=1000)$evppi, 
      323.7706, tol=1e-02)
    if (interactive()){
      evppi(chemo_nb, chemo_pars, pars=pars, method="inla", nsim=1000, plot_inla_mesh = TRUE)
    }
  })
}

if (requireNamespace(dbarts,quietly=TRUE)) { 
  test_that("EVPPI with BART",{
    pars <- c(pi2,rho)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="bart", nsim=1000)$evppi,
                 324.2451, tol=1e-02)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="bart", nsim=1000, ndpost=2000, se=TRUE)$evppi,
                 324.0516, tol=1e-02)
    
    ## TODO demo of close to EVPI 
    
  })
}

test_that("Standard errors in the earth method",{
  set.seed(1)
  ev1 <- evppi(chemo_nb, chemo_pars, pars=pi2, method="earth", se=TRUE, nsim=100)
  ev2 <- evppi(chemo_nb, chemo_pars, pars=pi2, method="earth", se=TRUE, nsim=500)
  expect_lt(ev2$se[1], ev1$se[1])
})
