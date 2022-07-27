## Copied from
## https://github.com/convoigroup/Chemotherapy_Book/blob/main/02_data/01_Data_Inputs.R
## https://github.com/convoigroup/Chemotherapy_Book/blob/main/02_data/02_Assumption_Inputs.R
## https://github.com/convoigroup/Chemotherapy_Book/blob/main/03_R/02_model_functions.R
## https://github.com/convoigroup/Chemotherapy_Book/blob/main/03_R/02_misc_functions.R
## on 20/07/2022 
## but with the constants pasted inside the functions 
## CJ TODO also change relative risk to log odds ratio 

################################################################################
#### Misc Functions for the Chemotherapy Model
################################################################################

## Function to transform values for mean and standard deviation into parameters 
## for a Beta distribution

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


################################################################################
#### Functions for the Chemotherapy Model
################################################################################


## Function to generate the PSA parameters
## Function to generate the PSA parameters
generate_psa_parameters <- function(n){
  ## n: The number of PSA simulations to be drawn
  time_horizon <- 50
  n_patients <- 111
  n_side_effects <- 52
  n_hospitalised <- 43
  n_died <- 8

  ## Parameters for the PSA distribution of the risk reduction of side effects
  logor_side_effects_mu <- log(0.54)
  logor_side_effects_sd <- 0.3

  ## Parameters for the PSA distribution of the recovery time for patients who are
  ## not hospitalised
  p_recovery_home_mu <- 0.21
  p_recovery_home_sd <- 0.03

  ## Parameters for the PSA distribution of the recovery time for patients who are
  ## hospitalised
  p_recovery_hosp_mu <- 0.03
  p_recovery_hosp_sd <- 0.0065

  ## Parameters for the PSA distribution of the costs of treating patients at home
  c_home_care_mu <- 830
  c_home_care_sd <- sqrt(150)

  ## Parameters for the PSA distribution of the costs of treating patients in
  ## hospital
  c_hospital_mu <- 2400
  c_hospital_sd <- sqrt(1880)

  ## Parameters for the PSA distribution of the one-off cost of death
  c_death_mu <- 1710
  c_death_sd <- sqrt(760)

  ## Parameters for the PSA distribution of the utility for recovered patients
  u_recovery_mu <- 0.98
  u_recovery_sd <- sqrt(0.001)

  ## Parameters for the PSA distribution of the utility of patients who are treated
  ## at home. 
  u_home_care_mu <- 0.7
  u_home_care_sd <- sqrt(0.02)

  ## Parameters for the PSA distribution of the utility of treating patients in 
  ## hospital
  u_hospital_mu <- 0.3
  u_hospital_sd <- sqrt(0.03)

  # Probability of side effects under treatment 1
  p_side_effects_t1 <- rbeta(n, 
                             1 + n_side_effects, 
                             1 + n_patients - n_side_effects)
  
  # Log odds of side effects on treatment 2
  logor_side_effects <- rnorm(n, logor_side_effects_mu, logor_side_effects_sd)
  # Odds of side effects on treatment 1
  odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  # Odds for side effects on treatment 2
  odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)

  # Probability of side effects under treatment 2
  p_side_effects_t2    <- odds_side_effects_t2 / (1 + odds_side_effects_t2)

  ## Variables to define transition probabilities
  # Probability that a patient is hospitalised over the time horizon
  p_hospitalised_total <- rbeta(n, 
                                1 + n_hospitalised, 
                                1 + n_side_effects - n_hospitalised)
  # Probability that a patient dies over the time horizon given they were 
  # hospitalised
  p_died <- rbeta(n, 1 + n_died, 1 + n_hospitalised - n_died)
  # Lambda_home: Conditional probability that a patient recovers considering 
  # that they are not hospitalised
  betapars <- betaPar(p_recovery_home_mu, p_recovery_home_sd)
  lambda_home <- rbeta(n, betapars$alpha, betapars$beta)
  # Lambda_hosp: Conditional probability that a patient recovers considering 
  # that they do not die
  betapars <- betaPar(p_recovery_hosp_mu, p_recovery_hosp_sd)
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
    p_side_effects_t1, 
    p_side_effects_t2,
    p_home_hospital, p_home_home, p_home_recover,
    p_hospital_hospital, p_hospital_recover, p_hospital_dead, 
    c_home_care, c_hospital, c_death,
    u_recovery, u_home_care, u_hospital,
    logor_side_effects,
    p_hospitalised_total, p_died,
    lambda_home, lambda_hosp)
  
  return(params_matrix)
}

## Function to calculate average time each patient with adverse events spends
## in the health states of the Markov model
calculate_state_occupancy_markov_model <- function(
  p_side_effects_t1, 
  p_side_effects_t2,
  p_home_home, p_home_hospital, p_home_recover,
  p_hospital_hospital, p_hospital_recover, p_hospital_dead)
  # All function arguments come from the generate_psa_parameters function except
  # time_horizon which is a model assumption
{ 
    time_horizon <- 50

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
  trace[1, 1, ] <- c(p_side_effects_t1, p_side_effects_t2)
  
  # Run the markove model over the time horizon
  for(i in 2:(time_horizon + 1)){   
    trace[, i, 1] <- trace[, i - 1, 1] %*% MM.mat 
    trace[, i, 2] <- trace[, i - 1, 2] %*% MM.mat
  }
  
  return(trace) # 4*16*2 array
}

## Function to calculate the costs and effects from our model
calculate_costs_effects <- function(
  p_side_effects_t1, 
  p_side_effects_t2,
  p_home_home, p_home_hospital, p_home_recover,
  p_hospital_hospital, p_hospital_recover, p_hospital_dead,
  c_home_care, c_hospital, c_death,
  u_recovery, u_home_care, u_hospital)
  # All function arguments come from the generate_psa_parameters function
{
    time_horizon <- 50

    ## Drug costs
    c_treatment_1 <- 120
    c_treatment_2 <- 10300

  # Calculate the trace matrix from the markov model function
  m_markov_trace <- calculate_state_occupancy_markov_model(
    p_side_effects_t1,
    p_side_effects_t2,
    p_home_home, p_home_hospital, p_home_recover,
    p_hospital_hospital, p_hospital_recover, p_hospital_dead)
  
  ## costs and effectiveness for four states
  c_state_vector <- c(c_home_care, c_hospital, 0, 0)
  u_state_vector <- c(u_home_care, u_hospital, u_recovery, 0)
  
  ## Estimate the cost of side effects from the Markov model
  c_side_effects <- array(NA, dim = 2)
  
  ## Average cost for both Soc and novel treatment per person
  ## (The cost includes one-off cost of death for all patients who died)
  c_side_effects[1] <- (sum(c_state_vector %*% m_markov_trace[, , 1]) + 
                          c_death * m_markov_trace[4, time_horizon + 1, 1])
  c_side_effects[2] <- (sum(c_state_vector %*% m_markov_trace[, , 2]) + 
                c_death * m_markov_trace[4, time_horizon + 1, 2])
  c_drug <- c(c_treatment_1, c_treatment_2)
  c_overall <- c(c_drug + c_side_effects)
  
  ## Total QALY of side effects for both Soc and novel treatment
  u_side_effects <- array(NA, dim = 2)
  u_side_effects[1] <- sum(u_state_vector %*% m_markov_trace[,,1])
  u_side_effects[2] <- sum(u_state_vector %*% m_markov_trace[,,2])
  ## QALY of total number of patients who do not experience adverse events for 15 days
  p_no_side_effects <- 1 - 
                          c(p_side_effects_t1,
                            p_side_effects_t2)
  u_no_side_effects <-  p_no_side_effects * u_recovery * (time_horizon + 1)
  
  ## Average effect for both Soc and novel treatment per person
  u_overall <- c(u_side_effects + u_no_side_effects)
  
  names(c_overall) <- paste0("cost",seq_along(c_overall))
  names(u_overall) <- paste0("eff",seq_along(u_overall))
  output <- array(NA, dim = c(1, length(u_overall), 2),
                  dimnames = list(NA, c("SoC", "Novel"),
                                  c("Effects", "Costs")))
  output[1, , 1] <- u_overall
  output[1, , 2] <- c_overall
  
  return(output)
}

calculate_net_benefit <- function(
  costs_effects,
  wtp)
{
    if(!is.null(dim(costs_effects))){
      nb <- wtp * costs_effects[, , 1] - 
        costs_effects[, , 2]
    }

  return(nb)
}

