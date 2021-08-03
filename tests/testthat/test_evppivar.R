## load_all(".")
#  load_all("voi")

test_that("single-parameter variance-based EVPPI",{
    evppivar(chemo_nb[,1], chemo_pars, pars="pi1")
    expect_equal(evppivar(chemo_nb[,1], chemo_pars[,"pi1"])$evppi,
                 evppivar(chemo_nb[,1], chemo_pars, pars="pi1")$evppi)
})

test_that("single-parameter variance-based EVPPI, alternative GAM basis",{
    expect_error({
    evppivar(chemo_nb[,1], chemo_pars, pars="pi1", gam_formula="s(pi1, bs='tp')")
    pi1 <- chemo_pars[,"pi1"]
    evppivar(chemo_nb[,1], pi1, gam_formula="s(pi1, bs='tp')")
    }, NA)
})

test_that("Variance-based EVPPI: list of parameters",{
    e1 <- evppivar(chemo_nb[,1], chemo_pars, pars="pi1")
    e3 <- evppivar(chemo_nb[,1], chemo_pars, pars=list("pi1","pi2"))
    expect_equal(e1$evppi, e3$evppi[e3$pars=="pi1"])
})
