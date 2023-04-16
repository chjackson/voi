lincomb_pars <- chemo_pars
lincomb_pars$constpar <- 1
lincomb_pars$ldpar <- lincomb_pars$p_side_effects_t1 + lincomb_pars$p_side_effects_t2

test_that("EVPPI with just one constant par",{
  ev1 <- evppi(chemo_nb, lincomb_pars, pars=c("constpar"))
  expect_equal(ev1$evppi, 0)
})

test_that("constant pars get dropped",{
  ev1 <- evppi(chemo_nb, lincomb_pars, pars=c("p_side_effects_t1", "constpar"))
  ev2 <- evppi(chemo_nb, lincomb_pars, pars=c("p_side_effects_t1"))
  expect_equal(ev1$evppi, ev2$evppi)
})

test_that("linearly dependent pars get dropped",{
  ev1 <- evppi(chemo_nb, lincomb_pars, 
               pars=c("p_side_effects_t1", "p_side_effects_t2","ldpar"))
  ev2 <- evppi(chemo_nb, lincomb_pars, 
               pars=c("p_side_effects_t2","ldpar"))
  expect_equal(ev1$evppi, ev2$evppi)
  ## EVPPI for t1, t2 is not the same - due to regression approx
  ev1 <- evppi(chemo_nb, lincomb_pars, 
               pars=c("p_side_effects_t1", "p_side_effects_t2","constpar","ldpar"))
})

