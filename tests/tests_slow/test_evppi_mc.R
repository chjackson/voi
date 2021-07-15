
## TODO feedback
## Shousn't need n_population  [ is there a draft on te chemo model in the book ? ]
## model function should separate costs and effects 
## was the model fn designed to return more than one evaluation? 
## Where are n_side_effects, n_patients, n_hospitalised, n_died, rr_side_effects_mu, rr_side_effects_sd from?  fill in for the moment 

## TODO coding
## Parallel 

## not sure bother about
## bias and uncertainty estimates, required number of samples
## Box 2 in Brennan 
## depends on what write about in book. 

## Latest version 
model_fn_chemo <- function(
  n_side_effects_pred_t1, 
  n_side_effects_pred_t2,
  p_home_home, p_home_hospital, p_home_recover,
  p_hospital_hospital, p_hospital_recover, p_hospital_dead,
  c_home_care, c_hospital, c_death,
  u_recovery, u_home_care, u_hospital,
  time_horizon,
  n_population) {
  wtp <- 20000
  ce <- cea_fn_chemo(  n_side_effects_pred_t1, 
                       n_side_effects_pred_t2,
                       p_home_home, p_home_hospital, p_home_recover,
                       p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                       c_home_care, c_hospital, c_death,
                       u_recovery, u_home_care, u_hospital,
                       time_horizon,
                       n_population)
  nb <- wtp * ce["e",] - ce["c",]
}

model_fn_chemo <- chemo_nb_fn
par_fn <- chemo_prior_pars

if (0) {
  evppi_mc(model_fn_chemo, par_fn, nouter=1000, ninner=100, k=20000, 
           pars="pi1" 
           #,mfargs = list(n_population=1000, time_horizon=15)
  )


}
