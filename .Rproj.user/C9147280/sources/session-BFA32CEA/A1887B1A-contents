###############################################################
### Expected Value of Partial Perfect Information Analysis  ###
###############################################################
### Packages
library(voi)
library(ggplot2)
library(dplyr)

## Run the model
source("04_analysis/02_baseline_model_output.R")

## Baseline Cost-Effectiveness Formatting
# The output from the cost-effectiveness model should be formatted for 
# the voi package.
# Use BCEA package to create a cost-effectiveness object.
chemotherapy_output <- list(e = m_costs_effects[, , "Effects"],
                            c = m_costs_effects[, , "Costs"],
                            k = seq(0, 50000, length.out = 501))

# EVPI
EVPI <- evpi(chemotherapy_output)
evpi.plot <- ggplot(EVPI, aes(x=k, y=evpi)) +
  geom_line()

# Select the WTP on the grid that is closest to the ICER
# Represents the maximum decision uncertainty (but not necessarily the max VOI)
wtp.max <- chemotherapy_output$k[
  which.min(abs(chemotherapy_output$k - ICER))
]

## Single parameter EVPPI
# Perform the calculations to determine the EVPPI for all parameters at all WTP
pars_all <- as.list(names(m_params))
ev_single <- evppi(outputs=chemotherapy_output, inputs=m_params, pars=pars_all)

ev_single <- ev_single %>%
  mutate(pars = rep(unlist(pars_all), each = length(chemotherapy_output$k)))

# Explore EVPPI at maximum uncertainty
ev_single %>%
  filter(k == wtp.max) %>%
  arrange(desc(evppi)) %>%
  mutate(evppi = round(evppi, 2))

# Plot key parameters
ev_single %>%
  filter(k == wtp.max) %>%
  plot(order = TRUE)

ggsave("06_figs/EVPPI_ICER.pdf")


## Calculate standard errors
pars_small <- list("u_home_care", "u_hospital", "p_home_recover", "lambda_home")
nb <- wtp.max * chemotherapy_output$e - chemotherapy_output$c
ev_small_se <- evppi(outputs=nb, inputs=m_params, pars=pars_small,
                     se = TRUE)
ev_small_se

# Plot key parameters across WTP
EVPI <- EVPI %>%
  mutate(pars = "all_parameters")

ev_single %>%
  filter(pars %in%c("logor_side_effects",
                    "p_home_recover",
                    "u_home_care",
                    "u_hospital",
                    "lambda_home",
                    "p_home_home",
                    "p_side_effects_t1")) %>%
  ggplot(aes(x=k, y=evppi, group=pars, col = pars)) +
  geom_line() + 
  geom_line(data = EVPI, aes(x = k, y = evpi), linetype = "dashed") + 
  theme_bw()

ggsave("06_figs/EVPPI_WTP.pdf")

## EVPPI Groups
# Randomised Trial
par_RCT <- list(
  "side_effects" = c("logor_side_effects"),
  "side_effects_and_follow_up" = c("logor_side_effects",
                                   "p_hospitalised_total","p_died",
                                   "lambda_home","lambda_hosp")
)
ev_RCT <- evppi(outputs=chemotherapy_output, inputs=m_params, pars=par_RCT)

# Explore at maximum willingness to pay
ev_RCT %>%
  filter(k == wtp.max) 

# Explore across willingness to pay
ev_RCT %>%
  ggplot(aes(x=k, y=evppi, group=pars, col = pars)) +
  geom_line() + 
  geom_line(data = EVPI, aes(x = k, y = evpi), linetype = "dashed") + 
  theme_bw()

ggsave("06_figs/EVPPI_RCT.pdf")

# All Studies
par_groups <- list(
  "side_effects" = c("logor_side_effects"),
  "hosp_trans_probs" = c("p_died","lambda_hosp"),
  "trans_probs" = c("p_side_effects_t1", "p_hospitalised_total", 
                    "lambda_home", "p_died","lambda_hosp"),
  "side_effects_and_follow_up" = c("logor_side_effects", "p_hospitalised_total",
                                   "p_died","lambda_home","lambda_hosp"),
  "costs" = c("c_home_care","c_hospital","c_death"),
  "utilities" = c("u_recovery","u_home_care","u_hospital")
)
ev_groups <- evppi(outputs=chemotherapy_output, inputs=m_params, pars=par_groups)

# Explore EVPPI at maximum uncertainty
ev_groups %>%
  filter(k == wtp.max) %>%
  arrange(desc(evppi)) %>%
  mutate(evppi = round(evppi, 2))

# Plot key parameters
ev_groups %>%
  filter(k == wtp.max) %>%
  plot(order = TRUE)

# Plot all groups across WTP
ev_groups %>%
  ggplot(aes(x=k, y=evppi, group=pars, col = pars)) +
  geom_line() + 
  geom_line(data = EVPI, aes(x = k, y = evpi), linetype = "dashed") + 
  theme_bw()


## Population Level EVPPI
pop_size <- 46000 * c(sum(1 / (1 + 0.035)^(0:5)),
                      sum(1 / (1 + 0.035)^(0:10)),
                      sum(1 / (1 + 0.035)^(0:15)))

ev_groups  %>%
  filter(k == 20000) %>%
  arrange(desc(evppi)) %>%
  mutate(evppi_5Y = evppi * pop_size[1],
         evppi_10Y = evppi * pop_size[2],
         evppi_15Y = evppi * pop_size[3])

ggsave("06_figs/EVPPI_groups.pdf")
