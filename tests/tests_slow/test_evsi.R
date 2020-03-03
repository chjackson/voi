# load_all(".")

#evppi(chemo_nb, chemo_pars, poi="pi1")
#evppi(chemo_cea, chemo_pars, poi="pi1")

##' Input is data frame of parameter simulations. One row per simulation, one column per parameter
##' Output is data frame of summarised future data simulations.  One row per simulation, one column per dimension of the summarised data. 
datagen_fn <- function(inputs, n=150){
    nsim <- nrow(inputs)
    with(inputs, { 
        X.SE1 <- rbinom(nsim, size=n, prob=pi1)
        X.SE2 <- rbinom(nsim, size=n, prob=pi2)
        X.N.hosp <- rbinom(nsim, X.SE1 + X.SE2, gamma.hosp)
        X.N.die <- rbinom(nsim, X.N.hosp, gamma.dead)
        N.amb <- X.SE1 + X.SE2 - X.N.hosp
        X.amb <- rgamma(nsim, N.amb, recover.amb) # sum of N.amb independent exponentials
        N.hosp <- X.N.hosp - X.N.die
        X.hosp <- rgamma(nsim, N.hosp, recover.hosp)
        data.frame(X.amb, X.SE1, X.SE2, X.N.hosp, X.hosp, X.N.die)
    })
}

##' Check the function for simulating summarised data works
Tdata <- datagen_fn(chemo_pars, n=150)
head(Tdata)

##' Calculate EVSI using GAM method (Strong et al)

set.seed(1) # method depends on simulating data, so make this test reproducible

##' Using default GAM formula te(everything, k=4) this is too computationally intensive. 
# evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn, nsim=1000) 

##' So we simplify the GAM formula 
##' Haven't thought about what an appropriate GAM formula might be in this example.  Can we give users guidance - related to how they expect their parameters to affect the outcome, given their model structure?   Similarly - how many simulations to use? 

##' Fewer parameters in the basis
gamf <- "te(X.amb, X.SE1, X.SE2, X.N.hosp, X.hosp, X.N.die, k=3, bs='cr')"
evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn, gam_formula=gamf, nsim=1000) # 19, runs in a few minutes

##' All smooths and no interactions 
gamf <- "s(X.amb) + s(X.SE1) + s(X.SE2) + s(X.N.hosp) + s(X.hosp) + s(X.N.die)"
evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn, gam_formula=gamf, nsim=1000) # 16, runs instantly 

##' Check works for different outcome formats (note results for same WTP will differ due to random seed)
evsi(chemo_nb, chemo_pars, datagen_fn=datagen_fn, gam_formula=gamf, nsim=1000) # 
evsi(chemo_cea, chemo_pars, datagen_fn=datagen_fn, gam_formula=gamf, nsim=1000) # 



##' Calculate EVSI using importance sampling method (Menzies)

##' Function to evaluate the likelihood for parameters given simulated data. 
##' Returns one row 
lik_chemo <- function(Y, inputs){
  loglik <-
      dbinom(Y[,"X.SE1"], size=150, inputs[,"pi1"], log=TRUE) +
      dbinom(Y[,"X.SE2"], size=150, inputs[,"pi2"], log=TRUE) +
      dbinom(Y[,"X.N.hosp"], Y[,"X.SE1"] + Y[,"X.SE2"], inputs[,"gamma.hosp"], log=TRUE) +
      dbinom(Y[,"X.N.die"],  Y[,"X.N.hosp"], inputs[,"gamma.dead"], log=TRUE) +
      log(inputs[,"recover.amb"]) * (Y[,"X.SE1"] + Y[,"X.SE2"] - Y[,"X.N.hosp"]) -
      inputs[,"recover.amb"] * Y[,"X.amb"] + 
      log(inputs[,"recover.hosp"]) * (Y[,"X.N.hosp"] - Y[,"X.N.die"]) -
      inputs[,"recover.hosp"] * Y[,"X.hosp"]
  exp(loglik)
}

##' Check the likelihood function works 
nsam <- 10
simdat <- datagen_fn(chemo_pars[1:nsam,])
simdat
lik_chemo(simdat[1,], chemo_pars[1:nsam,])

##' A GAM is needed to calculate the EVPPI in this method, though the default formula is intensive in this application and we override it. 
poi <- c("pi1","pi2","gamma.hosp","gamma.dead","recover.amb","recover.hosp") 
gamf <- sprintf("te(%s, k=3)", paste(poi, collapse=", "))  # runs very slowly with this.  9 
gamf <- paste0("s(",poi,")", collapse=" + ") # runs fast with this.  10 

evsi(chemo_nb, chemo_pars, method="is", poi=deadpoi, nsim=1000, 
     datagen_fn=datagen_fn, likelihood=lik_chemo, gam_formula=gamf)
