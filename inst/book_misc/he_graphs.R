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
         dln = plogitnorm(x, qlogis(0.4), 
                     (qlogis(0.6) - qlogis(0.2))/(2*qnorm(0.975))))

green <- viridis(5)[4]
purple <- viridis(5)[2]
pdf("~/work/voibook/voibook/Figures/02-healthecon/elic_compare.pdf", width=6, height=4)
ggplot(pdat, aes(x=x, y=dbet)) + 
  geom_line(lwd=2, col=green) + 
  geom_line(aes(y=dln), col=purple, lwd=2) + 
  geom_text(aes(x=0.55, y=3), label="Beta(8.9, 13.1)", col=green, hjust=0) + 
  geom_text(aes(x=0.55, y=2.6), label="Logit-normal(-0.4,0.46)", col=purple, hjust=0) + 
  theme_bw() + 
  xlab(expression(italic(p))) + ylab("Probability density") + 
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
dev.off()
