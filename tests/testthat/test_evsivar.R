
test_that("Variance-based EVSI with built-in study designs", {
    set.seed(1)
    evsivar(chemo_nb[,1], chemo_pars, study="trial_binary", poi=c("pi1", "pi2"), verbose=FALSE)
    evsivar(chemo_nb[,1], chemo_pars, study="binary", poi=c("pi1"), verbose=FALSE)
})
