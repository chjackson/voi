## load_all(".")
#  load_all("voi")

test_that("single-parameter variance-based EVPPI",{
    evppivar(chemo_nb[,1], chemo_pars, pars="pi1")
    expect_equal(evppivar(chemo_nb[,1], chemo_pars[,"pi1"]),
                 evppivar(chemo_nb[,1], chemo_pars, pars="pi1"))
})

test_that("single-parameter variance-based EVPPI, alternative GAM basis",{
    expect_error({
    evppivar(chemo_nb[,1], chemo_pars, pars="pi1", gam_formula="s(pi1, bs='tp')")
    pi1 <- chemo_pars[,"pi1"]
    evppivar(chemo_nb[,1], pi1, gam_formula="s(pi1, bs='tp')")
    }, NA)
})

