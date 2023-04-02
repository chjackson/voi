### TODO Check results match Anna's


###### (a) TRIAL THAT UPDATES LOG OR ####### 

# Data generation function
OR_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  p_side_effects_t2 <- inputs[, "p_side_effects_t2"]
  X1 <- rbinom(length(p_side_effects_t1), n, p_side_effects_t1)
  X2 <- rbinom(length(p_side_effects_t2), n, p_side_effects_t2)
  # Create odds ratio as summary statistic
  OR <- (n - X2) / X2 / ((n - X1) / X1)
  data_save <- data.frame(OR = OR)
  return(data_save)
}

# EVSI calculation using GAM regression.
evsi_OR <- evsi(outputs = chemo_nb,
                inputs = chemo_pars,
                pars = c("logor_side_effects"),
                pars_datagen = c("p_side_effects_t1", 
                                 "p_side_effects_t2"),
                n = seq(50, 500, by = 100),
                method = "gam",
                datagen_fn = OR_datagen_fn,
                par_fn = chemo_pars_fn)
evsi_OR

trial_binary_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  p_side_effects_t2 <- inputs[, "p_side_effects_t2"]
  X1 <- rbinom(length(p_side_effects_t1), n, p_side_effects_t1)
  X2 <- rbinom(length(p_side_effects_t2), n, p_side_effects_t2)
  data_save <- data.frame(X1 = X1, X2 = X2)
  return(data_save)
}

# Analysis function based on JAGS
OR_analysis_fn <- function(data, args, pars){
  X1 <- data$X1
  X2 <- data$X2
  
  data_jags <- list(X1 = X1,  X2 = X2, n = args$n,
                    n_side_effects = args$n_side_effects,
                    n_patients = args$n_patients,
                    logor_side_effects_mu = args$logor_side_effects_mu,
                    logor_side_effects_sd = args$logor_side_effects_sd)
  
  LogOR_trial <- "
  model { 
    # Probability of side effects under treatment 1
    p_side_effects_t1 ~ dbeta(1 + n_side_effects, 
                              1 + n_patients - n_side_effects)
    
    # Log odds of side effects on treatment 2
    logor_side_effects ~ dnorm(logor_side_effects_mu, 
                               logor_side_effects_sd)
    # Odds of side effects on treatment 1
    odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
    # Odds for side effects on treatment 2
    odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
    
    # Probability of side effects under treatment 2
    p_side_effects_t2 <- 
      odds_side_effects_t2 / (1 + odds_side_effects_t2)
    
    X1 ~ dbin(p_side_effects_t1, n)
    X2 ~ dbin(p_side_effects_t2, n)
  }
  "
  filein <- file.path(tempdir(), fileext="datmodel.txt")
  cat(LogOR_trial,file=filein)
  
  jagsmod <- rjags::jags.model(textConnection(LogOR_trial), 
                               data=data_jags, quiet=TRUE)
  update(jagsmod, 250, progress.bar="none")
  sam <- rjags::coda.samples(jagsmod, variable.names="logor_side_effects", 
                             n.iter=args$n.iter, progress.bar="none")
  logor_side_effects <- as.data.frame(sam[[1]])[,1]
  
  return(data.frame(logor_side_effects = logor_side_effects))
}

### NOTE analysis function should only return log OR, to match the model_fn 
### TODO delete upstream as not in book 

# EVSI calculation using the moment matching method.
analysis_args <- list(n_side_effects = chemo_constants$n_side_effects,
                      n_patients = chemo_constants$n_patients,
                      n = 500,
                      logor_side_effects_mu = chemo_constants$logor_side_effects_mu,
                      logor_side_effects_sd = chemo_constants$logor_side_effects_sd,
                      n.iter = 7500)
evsi_OR <- evsi(outputs = chemo_nb,
                inputs = chemo_pars,
                pars = c("logor_side_effects"),
                pars_datagen = c("p_side_effects_t1", "p_side_effects_t2"),
                n = seq(50, 500, by = 10),
                method = "mm",
                datagen_fn = OR_datagen_fn,
                model_fn = chemo_model_lor_nb,
                analysis_args = analysis_args, 
                analysis_fn = OR_analysis_fn, 
                par_fn = chemo_pars_fn, Q=20)
evsi_OR




###### (a) TRIAL THAT UPDATES LOG OR ####### 
