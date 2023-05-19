## Generate chemotherapy PSA inputs and outputs
## using functions built into package

## Work out data constants used in the model 

setwd("inst/Chemotherapy_Book")
source("02_data/01_Data_Inputs.R")
source("02_data/02_Assumption_Inputs.R")
setwd("../..")
vars <- c("n_patients", "n_side_effects", "n_hospitalised", "n_died",
  "c_death_mu", "c_death_sd", "c_home_care_mu", "c_home_care_sd", 
  "c_hospital_mu", "c_hospital_sd", "c_treatment_1", "c_treatment_2", 
  "dat_side_effs_SC", "logor_side_effects_mu", "logor_side_effects_sd", 
  "n_died", "n_hospitalised", "n_patients", "n_side_effects", "p_recovery_home_mu", 
  "p_recovery_home_sd", "p_recovery_hosp_mu", "p_recovery_hosp_sd", 
  "rate_longterm_mu", "rate_longterm_sd", "time_horizon", "time_horizon_longterm", 
  "u_home_care_mu", "u_home_care_sd", "u_hospital_mu", "u_hospital_sd", 
  "u_recovery_mu", "u_recovery_sd")
chemo_constants <- setNames(lapply(vars, get), vars)

set.seed(1)
nsim <- 10000
pars_used <- names(formals(chemo_model_nb))

## Generate PSA inputs
chemo_pars <- chemo_pars_fn(nsim)

## Generate model outputs at these inputs
ppars <- chemo_pars[pars_used]
effs <- costs <- matrix(nrow=nsim, ncol=2)
for (i in 1:nsim){
  ce <- do.call(chemo_model_cea, ppars[i,])
  effs[i,] <- ce[1,]
  costs[i,] <- ce[2,]
}
colnames(effs) <- colnames(costs) <- c("SoC","Novel")
chemo_nb <- as.data.frame(effs*20000 - costs)
chemo_cea <- list(e = effs, c = costs, k = 10000 * 1:5)

library(usethis)
use_data(chemo_constants, overwrite=TRUE)
use_data(chemo_pars, overwrite=TRUE)
use_data(chemo_cea, overwrite=TRUE)
use_data(chemo_nb, overwrite=TRUE)
