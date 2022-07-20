test_that("Built-in model functions work",{
  nsim <- 10
  pars_used <- names(formals(chemo_model_nb))
  ppars <- chemo_pars_fn(nsim)[pars_used]
  nb <- do.call(chemo_model_nb, ppars[1,])
  cea <- do.call(chemo_model_cea, ppars[1,])
  expect_equivalent(cea["Effects","SoC"]*20000 - cea["Costs","SoC"], nb["SoC"])
})
