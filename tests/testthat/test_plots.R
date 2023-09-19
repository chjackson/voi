pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"

test_that("EVPPI plots work",{
  expect_error({
    evtest <- evppi(chemo_nb, chemo_pars, pars=pi1)
    plot(evtest)
    evtest <- evppi(chemo_nb, chemo_pars, pars=list(pi1,pi2))
    plot(evtest)
    evtest <- evppi(chemo_cea, chemo_pars, pars=c(pi1))
    plot(evtest)
    evtest <- evppi(chemo_cea, chemo_pars, pars=c(pi1,pi2))
    plot(evtest)
    evtest <- evppi(chemo_cea, chemo_pars, pars=as.list(names(chemo_pars)))
    plot(evtest, top=6)
  }, NA)
})
