## Generate chemotherapy PSA inputs and outputs
## using functions built into package

library(devtools)
load_all(".")

set.seed(1)

nsim <- 10000
chemo_pars <- chemo_prior_pars(nsim)

cesim <- matrix(nrow=nsim, ncol=4)
ce_pars <- names(formals(chemo_costeff_fn))
for (i in 1:nsim){
    cesim[i,] <- as.numeric(t(do.call(chemo_costeff_fn, chemo_pars[i,ce_pars])))
}

chemo_cea <- list(e = cesim[,1:2], c = cesim[3:4], k = 10000 * 1:5)

chemo_nb <- chemo_cea$e*30000 - chemo_cea$c

use_data(chemo_pars, overwrite=TRUE)
use_data(chemo_cea, overwrite=TRUE)
use_data(chemo_nb, overwrite=TRUE)


## TODO needs replacing with the current functions 
## and all the tests redoing 
## Make a start to work out those constants 

dat_side_effs_SC <- read.csv("~/work/voibook/Chemotherapy_Book/01_data_raw/Side_Effects_SC.csv",
                             header = TRUE)

n_patients <- nrow(dat_side_effs_SC) # 111
n_side_effects <- sum(dat_side_effs_SC$Side_Effects) ## 52
n_hospitalised <- sum(dat_side_effs_SC$Hosp) ## 43
n_died <- sum(dat_side_effs_SC$Death) # 8

