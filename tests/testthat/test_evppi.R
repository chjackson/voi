## load_all(".")
#  load_all("voi")

pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"
rho <- "logor_side_effects"

test_that("single-parameter EVPPI",{
    evtest <- evppi(chemo_nb, chemo_pars, pars=pi2)
    expect_equal(evtest$evppi, 262.1109, tol=1e-05)
    evcea <- evppi(chemo_cea, chemo_pars, pars=pi2)
    expect_equal(evcea$evppi[evcea$k==20000], evtest$evppi, tol=1)
    expect_equal(evppi(chemo_nb, chemo_pars[,"p_side_effects_t2"])$evppi,
                 evtest$evppi)
    evcea <- evppi(chemo_cea, chemo_pars, pars=pi2, tidy=FALSE)
    expect_equal(evcea$evppi[evcea$k==20000], evtest$evppi, tol=1)
})

test_that("single-parameter EVPPI, alternative GAM basis",{
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pi2, gam_formula="s(p_side_effects_t2, bs='tp')")$evppi,
                 261.8708, tol=1e-04)
})

test_that("single-parameter EVPPI, GP",{
    expect_equal(
        evppi(chemo_nb, chemo_pars, pars=pi2, method="gp", nsim=100)$evppi,
        285.4876, tol=1e-03)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pi2, method="gp", nsim=100, gp_hyper_n=100)$evppi,
                 284.2202, tol=1e-03)
    gse <- evppi(chemo_nb, chemo_pars, pars=pi2, method="gp", nsim=100, se=TRUE, B=10)$se
    expect_true(is.numeric(gse))
})

test_that("single-parameter EVPPI, earth",{
    expect_equal(
        evppi(chemo_nb, chemo_pars, pars=pi2, method="earth")$evppi, 
        265.1299, tol=1e-03)
})

if (requireNamespace("INLA",quietly=TRUE)) { 
  test_that("EVPPI with INLA",{
    skip_on_cran()
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

if (requireNamespace("dbarts",quietly=TRUE)) { 
  test_that("EVPPI with BART",{
    skip_on_cran()
    pars <- c(pi2,rho)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="bart", nsim=1000)$evppi,
                 324.2451, tol=1e-02)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="bart", nsim=1000, ndpost=2000, se=TRUE)$evppi,
                 324.0516, tol=1e-02)
  })
}

test_that("Standard errors in the earth method",{
  skip_on_cran()
  set.seed(1)
  ev1 <- evppi(chemo_nb, chemo_pars, pars=pi2, method="earth", se=TRUE, nsim=100)
  ev2 <- evppi(chemo_nb, chemo_pars, pars=pi2, method="earth", se=TRUE, nsim=500)
  expect_lt(ev2$se[1], ev1$se[1])
})

test_that("multi-parameter EVPPI, gam",{
    pars <- c(pi2,rho)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="gam")$evppi,
                 333.4516, tol=1e-01)
})

test_that("multi-parameter EVPPI, earth",{
    pars <- c(pi2,rho)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="earth")$evppi,
                 334.0442, tol=1e-02)
})

test_that("multi-parameter EVPPI, s() instead of t() formula",{
    pars <- c(pi2,rho)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pars, method="gam", 
                       gam_formula="s(p_side_effects_t2) + s(logor_side_effects)")$evppi,
                 333.8265, tol=1e-03)
})

test_that("Multiple EVPPI calculations with the same call",{
  evtest <- evppi(chemo_nb, chemo_pars, pars=list(pi1, pi2))
  evtest2 <- evppi(chemo_nb, chemo_pars, pars=pi2)
  expect_equal(evtest$evppi[evtest$pars==pi2], evtest2$evppi)
  evtest <- evppi(chemo_nb, chemo_pars, pars=list(c(pi1, pi2), rho))
  evtest2 <- evppi(chemo_nb, chemo_pars, pars=c(pi1,pi2))
  expect_equal(evtest$evppi[evtest$pars=="p_side_effects_t1,p_side_effects_t2"], 
               evtest2$evppi)
  ## CEA format and multi pars 
  evtest <- evppi(chemo_cea, chemo_pars, pars=list(pi1, pi2))
  evtest2 <- evppi(chemo_nb, chemo_pars, pars=list(pi1, pi2))
  expect_equal(evtest$evppi[evtest$pars==pi2 & evtest$k==20000], 
               evtest2$evppi[evtest2$pars==pi2], tol=1)
  evtest <- evppi(chemo_cea, chemo_pars, pars=list(c(pi1, pi2), rho))
  evcea <- evppi(chemo_cea, chemo_pars, pars=c(pi1, pi2))
  expect_equal(evtest$evppi[evtest$pars=="p_side_effects_t1,p_side_effects_t2" & 
                              evtest$k==20000], evcea$evppi[evcea$k==20000], tol=1e-05)
  
  e1 <- evppi(chemo_cea, chemo_pars, pars=list(pi1, pi2))
  e2 <- evppi(chemo_cea, chemo_pars, pars=list(pi2))
  expect_equal(
    e1$evppi[e1$pars=="p_side_effects_t2" & e1$k==40000],
    e2$evppi[e2$k==40000]
  )
})

test_that("Strong and Oakley single parameter", {
    expect_error(evppi(chemo_nb, chemo_pars, pars=pi2, method="so"),
                 "`n.blocks` is required")
    expect_error(evppi(chemo_nb, chemo_pars, pars=c(pi1,pi2), method="so", n.blocks=20),
                 "only works for single-parameter")
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pi2, method="so", n.blocks=20)$evppi, 261.6701, tol=1e-02)
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pi2, method="so", n.blocks=40)$evppi, 262.3276, tol=1e-02)
    expect_equal(evppi(chemo_cea, chemo_pars, pars=pi2, method="so", n.blocks=20)$evppi[3], 118.0415, tol=1e-02)
    expect_equal(evppi(chemo_cea, chemo_pars, pars=pi2, method="so", n.blocks=40)$evppi[3], 118.3353, tol=1e-02)
})

test_that("Sadatsafavi et al single parameter", {
    expect_error(evppi(chemo_nb, chemo_pars, pars=c(pi1,pi2), method="sal", n.blocks=20),
                 "only works for single-parameter")
    evtest <- evppi(chemo_nb, chemo_pars, pars=pi2, method="sal")
    expect_equal(evtest$evppi, 262.8047, tol=1e-01)
    expect_equal(evppi(chemo_cea, chemo_pars, pars=pi2, method="sal")$evppi[2], 
                 evtest$evppi)
})

test_that("Standard errors for GAM",{
    set.seed(1)
    evtest <- evppi(chemo_nb, chemo_pars, pars=pi2, se=TRUE, B=10)
    expect_equal(evtest$se[1],  2.509941, tol=1e-01)
    evtest <- evppi(chemo_nb, chemo_pars, pars=list(pi1, pi2), se=TRUE, B=10)
    expect_equal(evtest$se[2], 4.274408, tol=1e-01)
})

test_that("Spaces in variable names", {
    cp <- chemo_pars
    names(cp)[names(cp)==pi2] <- "pi 1"
    expect_equal(evppi(chemo_nb, chemo_pars, pars=pi2)$evppi,
                 evppi(chemo_nb, cp, pars="pi 1")$evppi)
    
})

test_that("Variable names matching R built in objects", {
    cp <- chemo_pars
    names(cp)[names(cp)==pi2] <- "pi"
    expect_error(evppi(chemo_nb, cp, pars="pi"), "R internal constant")
    names(cp)[names(cp)=="pi"] <- "letters"
    expect_error(evppi(chemo_nb, cp, pars="letters"), "R internal constant")
})
