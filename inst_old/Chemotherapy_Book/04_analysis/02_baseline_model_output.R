################################################################################
#### Baseline Results for Chemotherapy Model
################################################################################

### Packages
library(colorspace)

### Run the Model
source("04_analysis/01_model_run.R")

### Mean costs and effects
mean_costs_effects <- apply(m_costs_effects, c(2, 3), mean)
write.csv(mean_costs_effects, "07_tables/Mean_Costs_Effects.csv")

### Incremental Cost-Effectiveness Ratio (ICER)
ICER <- mean(m_costs_effects[, 1, 2] - m_costs_effects[, 2, 2]) / 
  mean(m_costs_effects[, 1, 1] - m_costs_effects[, 2, 1])
ICER
write.csv(ICER, "07_tables/ICER.csv")

### Cost-Effectiveness Acceptability Curve
colours <- rainbow_hcl(12)
m_max <- apply(m_net_benefit, c(1,3), which.max)
pdf("06_figs/CEAC.pdf")
plot(wtp_seq, apply(m_max == 2, 2, mean), type = "l", lwd = 3,
     ylim = c(0,1),
     col = colours[1],
     xlab = "Willingness-to-pay",
     ylab = "Probability of Cost Effectiveness for Treatment 2",
     main = "Cost Effectiveness Acceptability Curve")
dev.off()

### Cost-Effectiveness Acceptability Frontier
m_ceaf <- cbind(apply(m_max == 1, 2, mean), apply(m_max == 2, 2, mean))
pdf("06_figs/CEAF.pdf")
plot(wtp_seq, m_ceaf[, 2], type = "l", lwd = 3,
     ylim = c(0,1),
     col = colours[1],
     xlab = "Willingness-to-pay",
     ylab = "Probability of Cost Effectiveness",
     main = "Cost Effectiveness Acceptability Frontier")
points(wtp_seq, m_ceaf[, 1], type = "l", col = colours[7], lwd = 3)

v_optimal_treatment <- apply(apply(m_net_benefit, c(2,3), mean), 2, which.max)
v_ceaf <- vector(length = n_wtp)
for(i in 1:n_wtp){
  v_ceaf[i] <- m_ceaf[i, v_optimal_treatment[i]]
}

points(wtp_seq, v_ceaf, pch = 19, cex = 0.5)
points(wtp_seq, v_ceaf, pch = 0)
legend("bottomright", c("Treatment 1", "Treatment 2", "Frontier"),
       lwd = c(2, 2, NA), pch = c(NA, NA, 0),
       col = c(colours[1], colours[7], "black"))
dev.off()

### Baseline willingness to pay = 20000
baseline_wtp <- which(wtp_seq == 20000)

# Average net benefit
mean_net_benefit <- apply(m_net_benefit[, , baseline_wtp], 2, mean)
mean_net_benefit
names(mean_net_benefit) <- c("Treatment 1", "Treatment 2")
mean_net_benefit / 20000
write.csv(mean_net_benefit, "07_tables/Mean_Net_Benefit.csv")

# INB
INB <- diff(mean_net_benefit)
INB
q_INB <- quantile(apply(m_net_benefit[, , baseline_wtp], 1, diff), prob = c(0.025, 0.975))
q_INB

# Optimal treatment
d_star <-  which.max(mean_net_benefit)
d_star

