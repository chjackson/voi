
setwd("inst/Chemotherapy_Book")
source("04_analysis/02_baseline_model_output.R")
source("06_figs/01_plotting_functions.R")
setwd("../..")

chemotherapy_output <- list(e = m_costs_effects[, "Effects", ],
                            c = m_costs_effects[, "Costs", ],
                            k = seq(0, 50000, length.out = 501))

chemo_cea_501 <- chemotherapy_output
use_data(chemo_cea_501, overwrite=TRUE)

# Data generation function
OR_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  logor_side_effects <- inputs[, "logor_side_effects"]
  
  # Odds for side effects for treatment 1
  odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  # Odds for side effects on treatment 2
  odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
  # Probability of side effects under treatment 2
  p_side_effects_t2    <- odds_side_effects_t2 / (1 + odds_side_effects_t2)
  
  # Data generation
  X1 <- rbinom(length(p_side_effects_t1), n, p_side_effects_t1)
  X2 <- rbinom(length(p_side_effects_t2), n, p_side_effects_t2)
  data_save <- data.frame(X1 = X1, X2 = X2)
  return(data_save)
}

# Analysis function based on JAGS
OR_analysis_fn <- function(data, args, pars){
  X1 <- data$X1
  X2 <- data$X2
  
  data_jags <- list(X1 = X1,
                    X2 = X2,
                    n = args$n,
                    n_side_effects = args$n_side_effects,
                    n_patients = args$n_patients,
                    logor_side_effects_mu = args$logor_side_effects_mu,
                    logor_side_effects_sd = args$logor_side_effects_sd)
  
  LogOR_trial <- function(){
    # Probability of side effects under treatment 1
    p_side_effects_t1 ~ dbeta(1 + n_side_effects, 
                              1 + n_patients - n_side_effects)
    
    # Log odds of side effects on treatment 2
    logor_side_effects ~ dnorm(logor_side_effects_mu, logor_side_effects_sd)
    # Odds of side effects on treatment 1
    odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
    # Odds for side effects on treatment 2
    odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
    
    # Probability of side effects under treatment 2
    p_side_effects_t2    <- odds_side_effects_t2 / (1 + odds_side_effects_t2)
    
    X1 ~ dbin(p_side_effects_t1, n)
    X2 ~ dbin(p_side_effects_t2, n)
  }
  
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  R2OpenBUGS::write.model(LogOR_trial,filein)
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  bugs.data <- R2jags::jags(
    data =  data_jags,
    parameters.to.save = pars,
    model.file = filein, 
    n.chains = 1, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  
  return(data.frame(logor_side_effects = bugs.data$BUGSoutput$sims.matrix[, pars[1]]))
}

chemo_evsi_or <- evsi(outputs = chemotherapy_output,
                      inputs = m_params,
                      pars = c("logor_side_effects"),
                      pars_datagen = c("p_side_effects_t1", "logor_side_effects"),
                      n = seq(50, 1500, by = 50),
                      method = "mm",
                      datagen_fn = OR_datagen_fn,
                      model_fn = calculate_costs_effects,
                      analysis_args = list(n_side_effects = n_side_effects,
                                           n_patients = n_patients,
                                           n = 500,
                                           logor_side_effects_mu = logor_side_effects_mu,
                                           logor_side_effects_sd = logor_side_effects_sd,
                                           n.iter = 5250),
                      analysis_fn = OR_analysis_fn, 
                      par_fn = generate_psa_parameters)

use_data(chemo_evsi_or, overwrite=TRUE)
