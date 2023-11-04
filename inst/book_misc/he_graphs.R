library(tidyverse)
library(viridis)

bet <- SHELF::fitdist(vals=c(0.2, 0.4, 0.6),
                      probs=c(0.025, 0.5, 0.975),
                      lower=0, upper=1)$Beta

## https://en.wikipedia.org/wiki/Logit-normal_distribution
dlogitnorm <- function(x, mu, sd){
  dnorm(qlogis(x), mu, sd) / (x*(1-x))
}

pdat <- data.frame(x = seq(0, 1, by=0.01)) %>%
  mutate(dbet = dbeta(x, bet[["shape1"]], bet[["shape2"]]),
         dln = dlogitnorm(x, qlogis(0.4), 
                     (qlogis(0.6) - qlogis(0.2))/(2*qnorm(0.975)))) %>%
  pivot_longer(cols = c("dbet", "dln"), names_to="dist", values_to = "dens") %>%
  mutate(dist = fct_recode(dist, "Beta(8.9, 13.1)"="dbet", "Logit-normal(-0.4,0.46)"="dln"))

green <- viridis(5)[4]
purple <- viridis(5)[2]
pdf("~/work/voibook/voibook/chapters/chapter01_modelling/elic_compare.pdf", width=6, height=4)
ggplot(pdat, aes(x=x, y=dens, col=dist, linetype=dist)) + 
  geom_line(lwd=2) + 
  scale_color_manual(name = "", values = c(green, purple)) +
  scale_linetype_manual(name = "", values=c(1, 3)) +
  theme_bw() + 
  theme(legend.position = c(0.8, 0.8) ) +
  xlab(expression(italic(p))) + ylab("Probability density") + 
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
dev.off()
