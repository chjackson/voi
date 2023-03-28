## Illustrations of GP regression 

library(ggplot2)
set.seed(1)
y <- sample(chemo_nb, 40) / 1000000
x <- sample(chemo_pars[,"logor_side_effects"], 40)
qplot(y, x)

## So can we fit a GP regression using inbuilt tools? 
## can we use gpFunc using hyperpar

# need to understand what gpFunc does and modularise it 
# sets: names of POI?  why not do this outside ? it does some error checking. 
# it removes any parameters that are constant.  Needs external warning 

inputs <- remove_constant_linear_cols(inputs)

xpred <- seq(min(x), max(x), by=0.001)
mod <- gp(y, x, Xpred=xpred)

# lower the first par to make less smooth
# nu controls the extent of measurement error. 
# Note this does not do any prediction outside the data points so the fit is jagged. 

mod <- gp(y, x, hyper = c(1, 0.5637912), Xpred=xpred)
mod <- gp(y, x, hyper = c(0.5, 0.5637912), Xpred=xpred)

mod1 <- gp(y, x, hyper = c(1, 0.001), Xpred=xpred)
mod2 <- gp(y, x, hyper = c(0.5, 0.001), Xpred=xpred)
mod3 <- gp(y, x, hyper = c(0.1, 0.001), Xpred=xpred)
mod4 <- gp(y, x, hyper = c(0.01, 0.001), Xpred=xpred)

dat <- data.frame(y, x)
datpred1 <- data.frame(x=xpred, fitted=mod1$pred, delta=1)
datpred2 <- data.frame(x=xpred, fitted=mod2$pred, delta=0.5)
datpred3 <- data.frame(x=xpred, fitted=mod3$pred, delta=0.1)
datpred4 <- data.frame(x=xpred, fitted=mod4$pred, delta=0.01)
datpred <- rbind(datpred1, datpred2, datpred3)
datpred$delta <- factor(datpred$delta)

cols <- RColorBrewer::brewer.pal(4, "Greys")[2:4]

pdf("Figures/04-evppi/gp_regression.pdf", width=6, height=3)  
ggplot(dat, aes(x=x, y=y)) +
  geom_line(data=datpred, aes(x=x, y=fitted, col=delta), lwd=1, alpha=0.8) + 
  scale_color_manual(breaks=levels(datpred$delta), values=cols) + 
  coord_cartesian(ylim=c(0.9, 1.01)) +
  geom_point(size=2) + 
  xlab("Predictor") + ylab("Outcome") + 
  labs(col=expression(delta))
dev.off()
