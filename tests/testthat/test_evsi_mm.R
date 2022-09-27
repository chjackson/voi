test_that("moment matching method",{
  evm <- evsi(outputs=chemo_nb, inputs=chemo_pars, 
              pars="p_side_effects_t1",
              method = "mm",
              study =  "binary",
              n = c(100), Q = 5, 
              analysis_args = list(a=53, b=60, n=100),
              model_fn = chemo_model_nb,  par_fn = chemo_pars_fn)
  expect_true(is.data.frame(evm) && nrow(evm)==1)
  
  evm <- evsi(outputs=chemo_cea, inputs=chemo_pars, 
              pars="p_side_effects_t1",
              method = "mm",
              study =  "binary",
              n = c(100), Q = 5, 
              analysis_args = list(a=53, b=60, n=100),
              model_fn = chemo_model_cea,  par_fn = chemo_pars_fn)
  expect_true(is.data.frame(evm) && nrow(evm)==5)
})

test_that("errors in moment matching method",{
  expect_error(
    evsi(outputs=chemo_cea, inputs=chemo_pars, pars="p_side_effects_t1",
         method = "mm", study =  "binary", n = c(100), Q = 5, 
         analysis_fn = analysis_binary, 
         analysis_args = list(a=53, b=60, n=100),
         model_fn = chemo_model_nb,  par_fn = chemo_pars_fn),
    "output of model_fn should have two rows"
  )

  datagen_fn <- function(inputs, n=100){
    nsim <- nrow(inputs)
    data.frame(x1 = rbinom(nsim, size=n, prob=inputs[,"p_side_effects_t1"]))
  }
  
  expect_error(  
    evsi(outputs=chemo_nb, inputs=chemo_pars, 
         pars="p_side_effects_t1",
         method = "mm",
         datagen_fn = datagen_fn, 
         model_fn = chemo_model_nb,  par_fn = chemo_pars_fn),
    "`analysis_fn` should be supplied")
  
  expect_error(  
    evsi(outputs=chemo_nb, inputs=chemo_pars, 
         pars="p_side_effects_t1",
         method = "mm",
         datagen_fn = datagen_fn, 
         model_fn = chemo_model_nb,  par_fn = chemo_pars_fn),
    "`analysis_fn` should be supplied")

})

