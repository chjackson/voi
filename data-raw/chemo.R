## Generate chemotherapy PSA inputs and outputs
## using functions built into package

## Work out data constants used in the model 
dat_side_effs_SC <- read.csv("inst/Chemotherapy_Book/01_data_raw/Side_Effects_SC.csv",
                             header = TRUE)
n_patients <- nrow(dat_side_effs_SC) # 111
n_side_effects <- sum(dat_side_effs_SC$Side_Effects) ## 52
n_hospitalised <- sum(dat_side_effs_SC$Hosp) ## 43
n_died <- sum(dat_side_effs_SC$Death) # 8
# These and other constants get saved inside the functions in R/model_chemo_source.R

set.seed(1)
nsim <- 10000
pars_used <- names(formals(chemo_model_nb))

## Generate PSA inputs
chemo_pars <- chemo_pars_fn(nsim)

## Generate model outputs at these inputs
ppars <- chemo_pars[pars_used]
nbs <- effs <- costs <- matrix(nrow=nsim, ncol=2)
for (i in 1:nsim){
  nbs[i,] <- do.call(chemo_model_nb, ppars[i,]) # wtp=20000
  ce <- do.call(chemo_model_cea, ppars[i,])
  effs[i,] <- ce["Effects",]
  costs[i,] <- ce["Costs",]
}
colnames(nbs) <- colnames(effs) <- colnames(costs) <- colnames(ce)

chemo_nb <- nbs
chemo_cea <- list(e = effs, c = costs, k = 10000 * 1:5)

use_data(chemo_pars, overwrite=TRUE)
use_data(chemo_cea, overwrite=TRUE)
use_data(chemo_nb, overwrite=TRUE)
