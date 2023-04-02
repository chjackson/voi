
if (requireNamespace(INLA,quietly=TRUE)) { 
  ## TODO condition these on the package being installed, and put in the slow section
  test_that("EVPPI with INLA",{
    expect_error(evppi(chemo_nb, chemo_pars, pars=pi2, method="inla", nsim=100), "2 or more parameters")
    pars <- c(pi2,rho)
    set.seed(1)
    expect_equal(
      evppi(chemo_nb, chemo_pars, pars=pars, method="inla", nsim=1000)$evppi, 
      323.7706, tol=1e-02)
  })
}

if (requireNamespace(dbarts,quietly=TRUE)) { 
  test_that("EVPPI with BART",{
    pars <- c(pi2,rho)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="bart", nsim=1000)$evppi,
                 324.2451, tol=1e-02)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="bart", nsim=1000, ndpost=2000, se=TRUE)$evppi,
                 324.0516, tol=1e-02)
  })
}
