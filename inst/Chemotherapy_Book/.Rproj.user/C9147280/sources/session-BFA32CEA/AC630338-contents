################################################################################
#### Analysis for the Chemotherapy Example
################################################################################
rm(list = ls())
set.seed(123)

### Load the Inputs
source("02_data/01_data_inputs.R")
source("02_data/02_Assumption_Inputs.R")

### Load the Functions
source("03_R/01_misc_functions.R")
source("03_R/02_model_functions.R")

### Generate the parameters for the PSA
n_psa_size <- 5000
m_params <- generate_psa_parameters(n_psa_size)

### Run the model
m_costs_effects <- array(NA, dim = c(n_psa_size, 2, 2))
for(s in 1:n_psa_size){
  v_params <- m_params[s, ]
  m_costs_effects[s, , ] <- do.call(calculate_costs_effects,
                                  v_params[names(formals(calculate_costs_effects))])
}

# Set the column names for the output
dimnames(m_costs_effects)[2:3] <- dimnames(
  do.call(calculate_costs_effects,
          v_params[names(formals(calculate_costs_effects))])
  )[2:3]

# Estimate net benefit for different willingness-to-pay values
n_wtp <- 51
wtp_seq <- seq(0, 50000, length.out = n_wtp)

# Initialise the net benefit matrix
m_net_benefit <- array(NA, dim = c(n_psa_size, 2, n_wtp))
w <- 1
for(wtp in wtp_seq){
  m_net_benefit[ , , w] <- calculate_net_benefit(m_costs_effects, wtp = wtp)
  w <- w + 1
}

