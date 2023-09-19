pi1 <- "p_side_effects_t1"
pi2 <- "p_side_effects_t2"

test_that("check_regression",{
  evtest <- evppi(chemo_nb, chemo_pars, pars=pi1, check=TRUE)
  et <- check_regression(evtest)
  expect_equal(et$AIC, 165744, tol=1)
  evtest <- evppi(chemo_nb, chemo_pars, pars=list(pi1, c(pi1, pi2)), check=TRUE)
  et2 <- check_regression(evtest,par="p_side_effects_t1,p_side_effects_t2")
  expect_true(et$AIC != et2$AIC)  

  evtest <- evppi(chemo_cea, chemo_pars, pars=c(pi1), check=TRUE)
  etc <- check_regression(evtest, outcome="costs")
  ete <- check_regression(evtest, outcome="effects")
  expect_true(etc$AIC != ete$AIC)
  evtest <- evppi(chemo_cea, chemo_pars, pars=list(pi1, c(pi1, pi2)), check=TRUE)
  etb <- check_regression(evtest, outcome="costs")
  expect_equal(etc$AIC, etb$AIC)
  et <- check_regression(evtest, outcome="costs", 
                         par="p_side_effects_t1,p_side_effects_t2")
  expect_true(et$AIC != etb$AIC)
})
