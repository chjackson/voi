test_that("errors in model functions", {
  expect_error(
    evsi(outputs=chemo_nb, inputs=chemo_pars, pars="p_side_effects_t2", 
         method="mm", study="binary", n=seq(10, 110, by=10), Q=10, 
         analysis_args = list(a=53, b=60),
         model_fn = chemo_model_cea, par_fn = chemo_pars_fn),
    "output of model_fn should be a vector if `outputs` is in cost-effectiveness format")
})  
