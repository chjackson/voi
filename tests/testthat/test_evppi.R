## load_all(".")
#  load_all("voi")

test_that("single-parameter EVPPI",{
    evtest <- evppi(chemo_nb, chemo_pars, pars="pi1")
    expect_equal(evtest, 3.75833452581155, tol=1e-05)
    evcea <- evppi(chemo_cea, chemo_pars, pars="pi1")
    expect_equal(evcea$evppi[evcea$k==30000], evtest)
    expect_equal(evppi(chemo_nb, chemo_pars[,"pi1"]),
                 evtest)
    evcea <- evppi(chemo_cea, chemo_pars, pars="pi1", tidy=FALSE)
    expect_equal(evcea[3], evtest)
})

test_that("single-parameter EVPPI, alternative GAM basis",{
    evppi(chemo_nb, chemo_pars, pars="pi1", gam_formula="s(pi1, bs='tp')")
})


test_that("single-parameter EVPPI, GP",{
    evppi(chemo_nb, chemo_pars, pars="pi1", method="gp", nsim=100)
    evppi(chemo_nb, chemo_pars, pars="pi1", method="gp", nsim=100, gp_hyper_n=100)
})

test_that("single-parameter EVPPI, earth",{
    evppi(chemo_nb, chemo_pars, pars="pi1", method="earth")
})

test_that("multi-parameter EVPPI, gam",{
    pars <- c("pi1","rho")
    evppi(chemo_nb, chemo_pars, pars=pars, method="gam") # 17
})

test_that("multi-parameter EVPPI, earth",{
    pars <- c("pi1","rho")
    evppi(chemo_nb, chemo_pars, pars=pars, method="earth") # 17
})

test_that("EVPPI with INLA",{
    expect_error(evppi(chemo_nb, chemo_pars, pars="pi1", method="inla", nsim=100), "2 or more parameters")
    pars <- c("pi1","rho")
    INLA:::inla.dynload.workaround()
    evppi(chemo_nb, chemo_pars, pars=pars, method="inla", nsim=1000) # 17
})

test_that("4-parameter EVPPI",{
pars <- c("pi1","rho","gamma.hosp","gamma.dead")
evppi(chemo_nb, chemo_pars, pars=pars, method="gam") # 17
})

test_that("multi-parameter EVPPI, s() instead of t() formula",{
pars <- c("pi1","rho")
evppi(chemo_nb, chemo_pars, pars=pars, method="gam", gam_formula="s(pi1) + s(rho)") # 17
})


test_that("Multiple EVPPI calculations with the same call",{
    evtest <- evppi(chemo_nb, chemo_pars, pars=list("pi1", "pi2"))
    expect_equal(evtest$evppi[evtest$pars=="pi1"], 3.75833452581155, tol=1e-05)
    evtest <- evppi(chemo_nb, chemo_pars, pars=list(c("pi1", "pi2"), "rho"))
    expect_equal(evtest$evppi[evtest$pars=="pi1,pi2"], 17.1065526416352, tol=1e-05)
    ## CEA format and multi pars 
    evtest <- evppi(chemo_cea, chemo_pars, pars=list("pi1", "pi2"))
    expect_equal(evtest$evppi[evtest$pars=="pi1" & evtest$k==30000], 3.75833452581155, tol=1e-05)
    evtest <- evppi(chemo_cea, chemo_pars, pars=list(c("pi1", "pi2"), "rho"))
    evcea <- evppi(chemo_cea, chemo_pars, pars=c("pi1", "pi2"))
    expect_equal(evtest$evppi[evtest$pars=="pi1,pi2" & evtest$k==30000], evcea$evppi[evcea$k==30000], tol=1e-05)
})
