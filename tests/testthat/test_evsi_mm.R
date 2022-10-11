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

test_that("errors in analysis function",{

  datagen_fn <- function(inputs, n=100){
    p1 <- inputs[,"p_side_effects_t1"]
    logor <- inputs[,"logor_side_effects"]
    odds1 <- p1 / (1 - p1)
    odds2 <- odds1 * exp(logor)
    p2 <- odds2 / (1 + odds2)
    data.frame(y1 = rbinom(10000, n, p1),
               y2 = rbinom(10000, n, p2))
  }

  expect_error(
    evsi(chemo_nb, chemo_pars, pars="p_side_effects_t1", method="mm",
         datagen_fn = datagen_fn, analysis_fn = function(){}, 
         n = 10, 
         model_fn = chemo_model_nb, par_fn = chemo_pars_fn),
    "`analysis_fn` should have arguments `data`,`args`,`pars` in that order")

  expect_error(
    evsi(chemo_nb, chemo_pars, pars="logor_side_effects", method="mm",
         datagen_fn = datagen_fn, analysis_fn = function(data,args,pars){}, 
         n = 10, 
         model_fn = chemo_model_nb, par_fn = chemo_pars_fn),
    "Parameters logor_side_effects returned by `analysis_fn` are not arguments to model_fn")

  analysis_fn_wrongpars <- function(data,args,pars){
    data.frame(wrongpars = rep(0, 100))
  }

  chemo_model_logor <- function(p_side_effects_t1, 
                                logor_side_effects,
                                p_home_home, p_home_hospital, p_home_recover,
                                p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                                c_home_care, c_hospital, c_death,
                                u_recovery, u_home_care, u_hospital){
    odds1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
    odds2 <- exp(logor_side_effects) * odds1
    p_side_effects_t2 <- odds2 / (1 + odds2)
    chemo_model_nb(p_side_effects_t1, 
                   p_side_effects_t2,
                   p_home_home, p_home_hospital, p_home_recover,
                   p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                   c_home_care, c_hospital, c_death,

                   u_recovery, u_home_care, u_hospital)
  }

  expect_error(
    evsi(chemo_nb, chemo_pars, pars="logor_side_effects", method="mm",
         datagen_fn = datagen_fn, analysis_fn = analysis_fn_wrongpars, 
         n = 10, 
         model_fn = chemo_model_logor, par_fn = chemo_pars_fn),
    "Parameters logor_side_effects not found in data frame returned by `analysis_fn`")
})
