###########################################################
### Expected Value of Perfect Information Analysis      ###
###########################################################

### Packages
library(voi)
library(BCEA)

## Run the model
source("04_analysis/01_model_run.R")

## Expected Value of Perfect Information - single WTP
# Specify willingness to pay
wtp_fix = 20000
# Calculate EVPI from net benefit
evpi(m_net_benefit[ , , wtp_seq == wtp_fix])

## Baseline Cost-Effectiveness Formatting
# The output from the cost-effectiveness model should be formatted for 
# the voi package.
# Use BCEA package to create a cost-effectiveness object.
chemotherapy_output <- bcea(e = m_costs_effects[, , "Effects"],
                            c = m_costs_effects[, , "Costs"])

## Expected Value of Perfect Information
# Calculate
EVPI <- evpi(chemotherapy_output)
# WTP = 20000
EVPI$evpi[EVPI$k == wtp_fix]

# Plot
plot(EVPI,
     xlab = "Willingness-to-Pay",
     ylab = "EVPI",
     main = "Expected Value of Perfect Information",
     type = "l")

