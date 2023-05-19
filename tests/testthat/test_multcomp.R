## Create a fake decision model output with three interventions 

pi2 <- "p_side_effects_t2"

set.seed(1)
nsim <- nrow(chemo_cea$e)
rr3_sim <- rgamma(nsim, 100, 100)
multcomp_cea <- chemo_cea
multcomp_cea$c <- as.data.frame(multcomp_cea$c)
multcomp_cea$e <- as.data.frame(multcomp_cea$e)
multcomp_cea$c$trt3 <- multcomp_cea$c$Novel * 1.2 * rr3_sim
multcomp_cea$e$trt3 <- multcomp_cea$e$Novel * 0.8 * rr3_sim      # new trt is worse 
multcomp_nb <- multcomp_cea$e*20000 - multcomp_cea$c
multcomp_pars <- cbind(chemo_pars, rr3_sim)

test_that("EVPPI for decision models with three decision options",{
  expect_equal(evpi(multcomp_nb), 417, tol=10)
  expect_equal(evppi(multcomp_nb, multcomp_pars, par=pi2)$evppi, 262, tol=1)
  expect_equal(evpi(multcomp_cea)$evpi[2], 417, tol=10)
  expect_equal(evppi(multcomp_cea, multcomp_pars, par=pi2)$evppi[2], 262, tol=1)
})


## Test the moment matching EVSI method
## Adapt the model function to three interventions

multcomp_model_cea <- function(p_side_effects_t1, p_side_effects_t2,
                                p_hospitalised_total, p_died,
                                lambda_home, lambda_hosp,
                                c_home_care, c_hospital, c_death,
                                u_recovery, u_home_care, u_hospital,
                                rate_longterm, rr3){
  ce <- chemo_model_cea(p_side_effects_t1 = p_side_effects_t1,
                        p_side_effects_t2 = p_side_effects_t2,
                        p_hospitalised_total = p_hospitalised_total,
                        p_died = p_died,
                        lambda_home = lambda_home,
                        lambda_hosp = lambda_hosp,
                        c_home_care = c_home_care,
                        c_hospital = c_hospital,
                        c_death = c_death,
                        u_recovery = u_recovery,
                        u_home_care = u_home_care,
                        u_hospital = u_hospital,
                        rate_longterm = rate_longterm)
  
  trt3 <- ce[,2] * rr3
  cbind(ce, trt3)
}

multcomp_pars_fn <- function(n){
  cbind(chemo_pars_fn(n), 
        rr3 = 1.2*rgamma(n, 100, 100))
}


test_that("moment matching method",{
  set.seed(100)
  expect_error(
    evsi(outputs=multcomp_nb, inputs=chemo_pars, 
              pars="p_side_effects_t1",
              method = "mm",  study =  "binary", n = c(100), Q = 5, 
              analysis_args = list(a=53, b=60, n=100),
              model_fn = chemo_model_nb,  par_fn = chemo_pars_fn),
    "Number of decision options")
  
  expect_error(
    evm <- evsi(outputs=chemo_cea, inputs=chemo_pars, 
              pars="p_side_effects_t1",
              method = "mm",
              study =  "binary",
              n = c(100), Q = 5, 
              analysis_args = list(a=53, b=60, n=100),
              model_fn = multcomp_model_cea,  par_fn = multcomp_pars_fn),
    "Number of decision options")

  set.seed(1)
  evm2 <- evsi(outputs=multcomp_cea, inputs=multcomp_pars, 
               pars="p_side_effects_t2",
               method = "mm", study =  "binary",
               n = c(100), Q = 5, 
               analysis_args = list(a=53, b=60, n=100),
               model_fn = multcomp_model_cea,  par_fn = multcomp_pars_fn)
  evm2
  expect_equal(evm2$evsi[2], 167, tol=1) 

  evm3 <- evsi(outputs=multcomp_cea, inputs=multcomp_pars, 
               pars="p_side_effects_t2",
               method = "gam", study =  "binary",
               analysis_args = list(a=53, b=60, n=100))
  expect_equal(evm3$evsi[2], 221, tol=1)
})
