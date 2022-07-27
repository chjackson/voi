pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"

test_that("Variance-based EVSI with built-in study designs", {
    expect_error({
    set.seed(1)
    evsivar(chemo_nb[,1], chemo_pars, study="trial_binary", pars=c(pi1, pi2), verbose=FALSE)
evsivar(chemo_nb[,1], chemo_pars, study="binary", pars=c(pi1), verbose=FALSE)
    },NA)
})

