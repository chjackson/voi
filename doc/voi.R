## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(1) 
nsam <- 10000
inputs <- data.frame(p1 = rnorm(nsam, 1, 1), 
                     p2 = rnorm(nsam, 0, 2))

## -----------------------------------------------------------------------------
outputs_nb <- data.frame(t1 = 0, 
                         t2 = inputs$p1 - inputs$p2)

## -----------------------------------------------------------------------------
outputs_cea <- list( 
  e = data.frame(t1 = 0, t2 = inputs$p1), 
  c = data.frame(t1 = 0, t2 = inputs$p2), 
  k = c(1, 2, 3)
)

## -----------------------------------------------------------------------------
decision_current <- 2
nb_current <- 1
decision_perfect <- ifelse(outputs_nb$t2 < 0, 1, 2)
nb_perfect <- ifelse(decision_perfect == 1, 0, outputs_nb$t2)
(evpi <- mean(nb_perfect) - nb_current)

## -----------------------------------------------------------------------------
opp_loss <- nb_perfect - nb_current
mean(opp_loss)

## -----------------------------------------------------------------------------
library(voi)
evpi(outputs_nb)
evpi(outputs_cea)

## -----------------------------------------------------------------------------
prob_correct <- 1 - pnorm(0, 1, sqrt(5))

## -----------------------------------------------------------------------------
mean_truncnorm <- function(mu, sig, lower=-Inf, upper=Inf){ 
  a <- (lower-mu)/sig
  b <- (upper-mu)/sig
  mu + sig * (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
}
enb_correct <- mean_truncnorm(1, sqrt(5), lower=0) 
mean_nb_perfect <- enb_correct * prob_correct
(evpi <- mean_nb_perfect - nb_current)

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1")
evppi(outputs_nb, inputs, pars=c("p1","p2"))

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars=list("p1","p2"))
evppi(outputs_nb, inputs, pars=list("p1",c("p1","p2")))

## -----------------------------------------------------------------------------
evppi(outputs_cea, inputs, pars=list("p1",c("p1","p2")))

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", method="gp", nsim=1000)

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", method="earth")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')
#  install.packages('splancs')
#  devtools::install_version("ldr", version = "1.3.3", repos = "http://cran.uk.r-project.org")

## ----eval=FALSE---------------------------------------------------------------
#  evppi(outputs_nb, inputs, pars=c("p1","p2"), method="inla")

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars=c("p1","p2"), method="gam", gam_formula="s(p1) + s(p2)")

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", n.blocks=20, method="so")

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", method="sal")

## -----------------------------------------------------------------------------
model_fn_nb <- function(p1, p2){ 
  c(0, p1 - p2) 
}

## -----------------------------------------------------------------------------
model_fn_cea <- function(p1, p2){ 
  rbind(e = c(0, p1), 
        c = c(0, p2)) 
}

## -----------------------------------------------------------------------------
par_fn <- function(n){
  data.frame(p1 = rnorm(n, 1, 1),
             p2 = rnorm(n, 0, 2))
}

## ----eval=FALSE---------------------------------------------------------------
#  evppi_mc(model_fn_nb, par_fn, pars="p1", ninner=1000, nouter=100)

## -----------------------------------------------------------------------------
datagen_normal <- function(inputs, n=100, sigma=1){
  data.frame(xbar = rnorm(n = nrow(inputs),
                          mean = inputs[,"p1"],
                          sd = sigma / sqrt(n)))
}

evsi(outputs_nb, inputs, datagen_fn = datagen_normal, n=10)
set.seed(1)
evsi(outputs_nb, inputs, datagen_fn = datagen_normal, n=100)
evsi(outputs_nb, inputs, datagen_fn = datagen_normal, n=1000)


## -----------------------------------------------------------------------------
evsi(outputs_nb, inputs, study = "normal_known", n=100, pars = "p1")
evsi(outputs_nb, inputs, study = "normal_known", n=1000, pars = "p1")

## -----------------------------------------------------------------------------
likelihood_normal <- function(Y, inputs, n=100, sig=1){
  mu <- inputs[,"p1"]
  dnorm(Y[,"xbar"], mu, sig/sqrt(n))
}

evsi(outputs_nb, inputs, datagen_fn = datagen_normal, likelihood = likelihood_normal, 
     n=100, pars = "p1", method="is", nsim=1000)

## -----------------------------------------------------------------------------
evsi(outputs_nb, inputs, study = "normal_known", n=100, pars = "p1", method="is", nsim=1000)

