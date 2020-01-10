
load_all(".")

#evppi(chemo_nb, chemo_pars, poi="pi1")
#evppi(chemo_cea, chemo_pars, poi="pi1")

##' Input is data frame of parameter simulations. One row per simulation, one column per parameter
##' Output is data frame of summarised future data simulations.  One row per simulation, one column per dimension of the summarised data. 
rfn <- function(inputs, n=150){
    nsim <- nrow(inputs)
    with(inputs, { 
        X.SE1 <- rbinom(nsim, size=n, prob=pi1)
        X.SE2 <- rbinom(nsim, size=n, prob=pi2)
        X.N.hosp <- rbinom(nsim, X.SE1 + X.SE2, gamma)
        X.N.die <- rbinom(nsim, X.N.hosp, gamma2)
        N.amb <- X.SE1 + X.SE2 - X.N.hosp
        X.amb <- rgamma(nsim, N.amb, recover.amb) # sum of N.amb independent exponentials
        N.hosp <- X.N.hosp - X.N.die
        X.hosp <- rgamma(nsim, N.hosp, recover.hosp)
        data.frame(X.amb, X.SE1, X.SE2, X.N.hosp, X.hosp, X.N.die)
    })
}

Tdata <- rfn(chemo_pars)
prepost.s <- gam(INB ~ te(X.amb, X.SE1, X.SE2, X.N.hosp, X.hosp, X.N.die, k=3), data=Tdata)

gamf <- "te(X.amb, X.SE1, X.SE2, X.N.hosp, X.hosp, X.N.die, k=3)"
evsi(chemo_nb, chemo_pars, rfn=rfn, gam_formula=gamf, nsim=1000) # 
evsi(chemo_cea, chemo_pars, rfn=rfn, gam_formula=gamf, nsim=1000) # 
evsi(chemo_nb, chemo_pars, rfn=rfn, gam_formula=gamf) # 21

poi <- c("pi1","pi2","gamma","gamma2","recover.amb","recover.hosp") 
gamf <- sprintf("te(%s, k=3)", paste(poi, collapse=", "))
gamf <- paste0("s(",poi,")", collapse=" + ")

lik_chemo <- function(Y, inputs){
  loglik <-
      dbinom(Y[,"X.SE1"], size=150, inputs[,"pi1"], log=TRUE) +
      dbinom(Y[,"X.SE2"], size=150, inputs[,"pi2"], log=TRUE) +
      dbinom(Y[,"X.N.hosp"], Y[,"X.SE1"] + Y[,"X.SE2"], inputs[,"gamma"], log=TRUE) +
      dbinom(Y[,"X.N.die"],  Y[,"X.N.hosp"], inputs[,"gamma2"], log=TRUE) +
      log(inputs[,"recover.amb"]) * (Y[,"X.SE1"] + Y[,"X.SE2"] - Y[,"X.N.hosp"]) -
      inputs[,"recover.amb"] * Y[,"X.amb"] + 
      log(inputs[,"recover.hosp"]) * (Y[,"X.N.hosp"] - Y[,"X.N.die"]) -
      inputs[,"recover.hosp"] * Y[,"X.hosp"]
  exp(loglik)
}

load_all(".")

nsam <- 1000

simdat <- generate_data(chemo_pars[1:nsam,], rfn=rfn)
lik_chemo(simdat[1,], chemo_pars[1:nsam,])

evsi(chemo_nb[1:nsam,], chemo_pars[1:nsam,], method="is", poi=poi, rfn=rfn, likelihood=lik_chemo, gam_formula=gamf)
