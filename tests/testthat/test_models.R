## Test built-in example decision model functions 

# source("~/uncertainty/evi/EVSI_Chemotherapy_Model-master/chemotherapy_model.R")
## checked dists of all inputs and outputs match Anna's version 

if (0) { # TODO 
nsim <- 100
ppars <- chemo_prior_pars(nsim)
cesim <- matrix(nrow=nsim, ncol=4)
for (i in 1:nsim){
    cesim[i,] <- do.call(chemo_costeff_fn, ppars[i,ce_pars])
}

}
