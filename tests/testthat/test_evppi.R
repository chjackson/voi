## load_all(".")
#  load_all("voi")

test_that("single-parameter EVPPI",{
    evtest <- evppi(chemo_nb, chemo_pars, pars="pi1")
    expect_equal(evtest, 3.75833452581155, tol=1e-05)
    evcea <- evppi(chemo_cea, chemo_pars, pars="pi1")
    expect_equal(evcea[3], evtest)
    expect_equal(evppi(chemo_nb, chemo_pars[,"pi1"]),
                 evtest)
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
