### BUILT-IN STUDY DESIGNS 

evppi(outputs=chemo_nb, inputs=chemo_pars, 
      pars=c("p_side_effects_t1","p_side_effects_t2"), method="gam")

# EVSI calculation using GAM regression.
evsi_builtin_rb <- evsi(outputs = chemo_nb,
                        inputs = chemo_pars,
                        study = "trial_binary",
                        pars = c("p_side_effects_t1", 
                                 "p_side_effects_t2"),
                        n = seq(50, 500, by = 50),
                        method = "gam")
evsi_builtin_rb

# EVSI calculation using Importance Sampling
evsi_builtin_is <- evsi(outputs = chemo_nb,
                        inputs = chemo_pars,
                        study = "trial_binary",
                        pars = c("p_side_effects_t1", 
                                 "p_side_effects_t2"),
                        n = 1000,
                        method = "is")
evsi_builtin_is 

# Beta prior for standard care is set using the number of events
beta_params_t1 <- c(1 + chemo_constants$n_side_effects,
                    1 + chemo_constants$n_patients - chemo_constants$n_side_effects)
# Beta prior for the novel intervention is approximated from the 
# mean and standard deviation of the PA distribution for the 
# probability of side effects.
beta_params_t2 <- betaPar(mean(chemo_pars$p_side_effects_t2),
                          sd(chemo_pars$p_side_effects_t2))
# EVSI calculation with moment matching method
evsi_builtin_mm <- evsi(outputs = chemo_nb,
                        inputs = chemo_pars,
                        study = "trial_binary",
                        pars = c("p_side_effects_t1", "p_side_effects_t2"),
                        n = seq(50, 500, by = 50),
                        method = "mm",
                        model_fn = chemo_model_nb,
                        analysis_args = list(a1 = beta_params_t1[1],
                                             b1 = beta_params_t1[2],
                                             a2 = beta_params_t2$alpha,
                                             b2 = beta_params_t2$beta),
                        par_fn = chemo_pars_fn)




### BESPOKE STUDY DESIGNS 

if (requireNamespace(truncnorm, quietly=TRUE)){

# Data generation function - aggregate data (for regression method)
utility_datagen_fn_agg <- function(inputs, n = 500){
  dat_indiv <- utility_datagen_fn_indiv(inputs, n = n)
  X_hospital_mean <- rowMeans(dat_indiv)
  data_save_dat <- data.frame(X_hospital_mean = X_hospital_mean)
  return(data_save_dat)
}

# Data generation function - individual data (for other EVSI methods)
# Vectorised version: does not run any faster than the version in Chemotherapy_Book 
utility_datagen_fn_indiv <- function(inputs, n = 500){
  m_hospital <- rep(inputs[, "u_hospital"], n[1])
  sd_hospital <- rep(inputs[, "sd_iid_hospital"], n[1])
  X_hospital <- truncnorm::rtruncnorm(n[1], mean=m_hospital, 
                                      sd=sd_hospital, a=-Inf, b=1)
  X_hospital <- matrix(X_hospital, nrow=nrow(inputs), ncol=n[1])
  data_save_dat <- as.data.frame(cbind(X_hospital = X_hospital))
  return(data_save_dat)
}

chemo_pars$sd_iid_hospital <- runif(nrow(chemo_pars), 0.00001, 0.4)

## Regression Based Method  ( 12 sec )
evsi_utility_rb <- evsi(outputs = chemo_nb, inputs = chemo_pars,
                        pars = c("u_hospital"),  n=1000, 
                        n = seq(50, 1050, by = 500),
                        method = "gam",
                        datagen_fn = utility_datagen_fn_agg)


## Moment Matching Method
# Analysis function based on JAGS
utility_analysis_fn <- function(data, args, pars){
  # Create the data list for JAGS
  data_jags <- list(X_hospital = as.vector(data),
                    n = args$n,
                    alpha_hospital = betaPar(
                      args$u_hospital_mu,
                      args$u_hospital_sd
                    )$alpha,
                    beta_hospital = betaPar(
                      args$u_hospital_mu,
                      args$u_hospital_sd
                    )$beta)
  
  trial <- "
  model { 
    for(i in 1:n){
      X_hospital[i] ~ dnorm(u_hospital, tau_hospital)T(, 1)
    }
    u_hospital ~ dbeta(alpha_hospital, beta_hospital)
    sd_hospital ~ dunif(0.00001, 0.4)
    tau_hospital <- 1 / sd_hospital ^ 2
  }
  "

  jagsmod <- rjags::jags.model(textConnection(trial), 
                               data=data_jags, quiet=TRUE)
  update(jagsmod, 250, progress.bar="none")
  print(args$n.iter)
  sam <- rjags::coda.samples(jagsmod, variable.names="u_hospital", 
                             n.iter=args$n.iter, progress.bar="none")
  u_hospital <- as.data.frame(sam[[1]])[,1]
  return(data.frame(u_hospital = u_hospital))
}

analysis_args <- list(n = 30,
                      u_hospital_mu = chemo_constants$u_hospital_mu,
                      u_hospital_sd = chemo_constants$u_hospital_sd,
                      n.iter = 500)

set.seed(1)
evsi_utility_mm <- evsi(outputs = chemo_nb,
                        inputs = chemo_pars,
                        pars = c("u_hospital"),
                        pars_datagen = c("u_hospital","sd_iid_hospital"),
                        n = 500, 
                        method = "mm",  nsim = 1000, 
                        datagen_fn = utility_datagen_fn_indiv,
                        model_fn = chemo_model_nb,
                        analysis_args = analysis_args,
                        analysis_fn = utility_analysis_fn, 
                        par_fn = chemo_pars_fn,
                        Q = 10, verbose=TRUE)

## Importance Sampling
# Likelihood function - vectorised - faster 
utility_likelihood <- function(Y, inputs){
  data_vec <- unlist(Y)
  ndat <- length(data_vec)
  data_vec <- rep(data_vec, each = nrow(inputs))
  m_hospital <- rep(inputs[, "u_hospital"], ndat)
  sd_hospital <- rep(inputs[, "sd_iid_hospital"], ndat)
  ll <- log(truncnorm::dtruncnorm(data_vec, mean=m_hospital, sd=sd_hospital,
                              a=-Inf, b=1))
  ll_mat <- matrix(ll, nrow=nrow(inputs), ncol=ndat)
  ll <- exp(rowSums(ll_mat))
  return(ll)
}

set.seed(1)
evsi_utility <- evsi(outputs = chemo_nb,
                     inputs = chemo_pars,
                     pars = c("u_hospital"),
                     pars_datagen = c("u_hospital", "sd_iid_hospital"),
                     # n = seq(50, 1000, by = 200),
                     n = 50, 
                     method = "is",
                     nsim = 500, 
                     datagen_fn = utility_datagen_fn_indiv,
                     likelihood = utility_likelihood, verbose=TRUE)

}
