################################################################################
#### Functions for the Chemotherapy Model
################################################################################

## Function to generate the PSA parameters
pars_fn_chemo <- function(n,n_population){
  # n: The number of PSA simulations to be drawn

    ## CJ FILL IN THESE
    n_side_effects <- 27
    n_patients <- 111
    n_hospitalised <- 17
    n_died <- 1 
    rr_side_effects_mu <- 0.65 
    rr_side_effects_sd <- 0.1
    time_horizon <- 15
    c_home_care_mu <- 2300
    c_home_care_sd <- 90
    c_hospital_mu <- 6500
    c_hospital_sd <- 980
    c_death_mu <- 4200
    c_death_sd <- 560
    u_recovery_mu <- 0.98
    u_recovery_sd <- 0.01 
    u_home_care_mu <- 0.5 
    u_home_care_sd <- 0.02 
    u_hospital_mu <- 0.2 
    u_hospital_sd <- 0.03 
    n_population <- 1000 
  
  # Probability of side effects under treatment 1
  p_side_effects_t1 <- rbeta(n, 
                             1 + n_side_effects, 
                             1 + n_patients - n_side_effects)
  
  # Relative risk of side effects on treatment 2
  rr_side_effects <- rnorm(n, rr_side_effects_mu, rr_side_effects_sd)
  # Probability of side effects under treatment 2
  p_side_effects_t2 <- rr_side_effects * p_side_effects_t1

  # Predictive distribution of the number of side effects
  n_side_effects_pred_t1 <- rbinom(n, n_population, p_side_effects_t1)
  n_side_effects_pred_t2 <- rbinom(n, n_population, p_side_effects_t2)

  
  ## Variables to define transition probabilities
  # Probability that a patient is hospitalised over the time horizon
  p_hospitalised_total <- rbeta(n, 
                                1 + n_hospitalised, 
                                1 + n_side_effects - n_hospitalised)
  # Probability that a patient dies over the time horizon given they were 
  # hostpialised
  p_died <- rbeta(n, 1 + n_died, 4 + p_hospitalised_total - n_died)
  # Lambda_home: Conditional probability that a patient recovers considering 
  # that they are not hospitalised
  betapars <- betaPar(0.45, 0.02)
  lambda_home <- rbeta(n, betapars$alpha, betapars$beta)
  # Lambda_hosp: Conditional probability that a patient recovers considering 
  # that they do not die
  betapars <- betaPar(0.35, 0.02)
  lambda_hosp <- rbeta(n, betapars$alpha, betapars$beta)
  
  ## Transition Probabilities
  p_home_hospital <- 1 - (1 - p_hospitalised_total) ^ (1 / time_horizon)
  p_home_home <- (1 - lambda_home) * (1 - p_home_hospital)
  p_home_recover <- lambda_home * (1 - p_home_hospital)
  p_hospital_dead <- 1 - (1 - p_died) ^ (1 / time_horizon)
  p_hospital_hospital <- (1 - lambda_hosp) * (1 - p_hospital_dead)
  p_hospital_recover <- lambda_hosp * (1 - p_hospital_dead)
        
  ## Health State Costs
  lnpars <- lognPar(c_home_care_mu, c_home_care_sd)
  c_home_care <- rlnorm(n, lnpars$meanlog, lnpars$sdlog)
  lnpars <- lognPar(c_hospital_mu, c_hospital_sd)
  c_hospital <- rlnorm(n, lnpars$meanlog, lnpars$sdlog)
  lnpars <- lognPar(c_death_mu, c_death_sd)
  c_death <- rlnorm(n, lnpars$meanlog, lnpars$sdlog)
  
  ## Health Utilities
  betapars <- betaPar(u_recovery_mu, u_recovery_sd)
  u_recovery <- rbeta(n, betapars$alpha, betapars$beta)
  betapars <- betaPar(u_home_care_mu, u_home_care_sd)
  u_home_care <- rbeta(n, betapars$alpha, betapars$beta)
  betapars <- betaPar(u_hospital_mu, u_hospital_sd)
  u_hospital <- rbeta(n, betapars$alpha, betapars$beta)
  
  # Specify a matrix containing all the parameters
  params_matrix <- data.frame(
    n_side_effects_pred_t1,
    n_side_effects_pred_t2,
    p_home_hospital, p_home_home, p_home_recover,
    p_hospital_hospital, p_hospital_recover, p_hospital_dead, 
    c_home_care, c_hospital, c_death,
    u_recovery, u_home_care, u_hospital,
    p_side_effects_t1, 
    p_side_effects_t2,
    rr_side_effects,
    p_hospitalised_total, p_died,
    lambda_home, lambda_hosp)
  
  return(params_matrix)
}

## Function to calculate average time each patient with adverse events spends
## in the health states of the Markov model
calculate_state_occupancy_markov_model <- function(
  n_side_effects_pred_t1, 
  n_side_effects_pred_t2,
  p_home_home, p_home_hospital, p_home_recover,
  p_hospital_hospital, p_hospital_recover, p_hospital_dead,
  time_horizon)
  # All function arguments come from the generate_psa_parameters function except
  # time_horizon which is a model assumption
{ 
  ## Markov transition probability matrix 
  ## States: Home care, Hospital care, Recovery, Death
  MM.mat <- matrix(c(p_home_home, p_home_hospital, p_home_recover, 0,
                     0, p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                     0, 0, 1, 0,
                     0, 0, 0, 1),
                   nrow = 4, ncol = 4, byrow = TRUE)
  
  ## Number of patients in each state for each time point 
  ## 3 dimensions: number of states, number of time points, 
  ## number of treatment options
  trace <- array(0, dim = c(4, time_horizon + 1, 2)) 
  # Initialise with the predicted number of side effects in the population
  trace[1, 1, ] <- c(n_side_effects_pred_t1,n_side_effects_pred_t2)
  
  # Run the markove model over the time horizon
  for(i in 2:(time_horizon + 1)){   
    trace[, i, 1] <- trace[, i - 1, 1] %*% MM.mat 
    trace[, i, 2] <- trace[, i - 1, 2] %*% MM.mat
  }
  
  return(trace) # 4*16*2 array
}

## Function to calculate the costs and effects from our model
cea_fn_chemo <- function(
  n_side_effects_pred_t1, 
  n_side_effects_pred_t2,
  p_home_home, p_home_hospital, p_home_recover,
  p_hospital_hospital, p_hospital_recover, p_hospital_dead,
  c_home_care, c_hospital, c_death,
  u_recovery, u_home_care, u_hospital,
  time_horizon,
  n_population)
  # All function arguments come from the generate_psa_parameters function except
  # time_horizon which is a model assumption
{

    # CJ ADDED THESE
    c_treatment_1 <- 110 
    c_treatment_2 <- 420 

  # Calculate the trace matrix from the markov model function
  m_markov_trace <- calculate_state_occupancy_markov_model(
    n_side_effects_pred_t1,
    n_side_effects_pred_t2,
    p_home_home, p_home_hospital, p_home_recover,
    p_hospital_hospital, p_hospital_recover, p_hospital_dead,
    time_horizon)
  
  ## costs and effectiveness for four states
  c_state_vector <- c(c_home_care, c_hospital, 0, 0)
  u_state_vector <- c(u_recovery, u_home_care, u_hospital, 0)
  
  ## Estimate the cost of side effects from the Markov model
  c_side_effects <- array(NA, dim = 2)
  
  ## Average cost for both Soc and novel treatment per person
  ## (The cost includes one-off cost of death for all patients who died)
  c_side_effects[1] <- (sum(c_state_vector %*% m_markov_trace[, , 1]) + 
                          c_death * m_markov_trace[4, time_horizon + 1, 1])/
    n_population 
  c_side_effects[2] <- (sum(c_state_vector %*% m_markov_trace[, , 2]) + 
                c_death * m_markov_trace[4, time_horizon + 1, 2])/
    n_population 
  c_drug <- c(c_treatment_1, c_treatment_2)
  c_overall <- c(c_drug + c_side_effects)
  
  ## Total QALY of side effects for both Soc and novel treatment
  u_side_effects <- array(NA, dim = 2)
  u_side_effects[1] <- sum(u_state_vector %*% m_markov_trace[,,1])
  u_side_effects[2] <- sum(u_state_vector %*% m_markov_trace[,,2])
  ## QALY of total number of patients who do not experience adverse events for 15 days
  n_no_side_effects <- n_population - 
                          c(n_side_effects_pred_t1,
                            n_side_effects_pred_t2)
  u_no_side_effects <-  n_no_side_effects * u_recovery * (time_horizon + 1)
  
  ## Average effect for both Soc and novel treatment per person
  u_overall <- c(u_side_effects + u_no_side_effects)/
    n_population           
  
  names(c_overall) <- paste0("cost",seq_along(c_overall))
  names(u_overall) <- paste0("eff",seq_along(u_overall))
  
  return(rbind("c"=c_overall, "e"=u_overall))
}

calculate_net_benefit <- function(
  costs_effects,
  wtp=20000)
{
    nb <- wtp * costs_effects["e",] - 
      costs_effects["c",]
    return(nb)
}


betaPar <- function(m, s) {
  # m:  Mean of the Beta distribution
  # m: Standard deviation of the Beta distribution
  
  var <- s ^ 2
  alpha <- ((1 - m) / var - 1 / m) * m ^ 2
  beta <- alpha * (1 / m - 1)
  
  return(
    list(alpha = alpha, beta = beta)
  )
}

## Function to transform values for mean and standard deviation into parameters 
## for a Log-Normal distribution

lognPar <- function(m,s) {
  # m: Mean of Log-Normal distribution
  # s: Standard deiviation of Log-Normal distribution
  
  var <- s^2
  meanlog <- log(m) - 0.5 * log(1 + var/m^2)
  varlog <- log(1 + (var/m^2))
  sdlog <- sqrt(varlog)
  
  return(
    list(meanlog = meanlog, sdlog = sdlog)
  )
}
