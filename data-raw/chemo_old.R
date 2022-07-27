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
