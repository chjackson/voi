## Chemotherapy health economic model
##
## Functions to generate model inputs and outputs based on current information


##' Chemotherapy health economic model parameters
##'
##' Generate a sample from the uncertainty distribution of the
##' parameters.
##'
##' This is the state of knowledge under current information, prior to
##' any further research.  The distributions are all of standard
##' analytic forms.  The distributions of pi1, gamma.hosp and
##' gamma.dead are obtained from simple beta/binomial conjugate
##' Bayesian inference based on current data.
##'
##' @param n Number of samples to draw 
##' 
chemo_prior_pars <- function(n){
    num.pat <- 111 # Number of patients (observed data)  
    num.se <- 27   # Number of patients with side effects
    num.hosp <- 17  # Number of patients require hospital care
    num.dead <- 1   # Number of Deaths
    pi1 <- rbeta(n, 1 + num.se, 1 + num.pat - num.se)
    rho <- rnorm(n, 0.65, 0.1)
    pi2 <- rho * pi1
    n.pred <- 1000
    TH <- 15
    SE1 <- rbinom(n, n.pred, pi1)
    SE2 <- rbinom(n, n.pred, pi2)
    gamma.hosp <- rbeta(n, 1 + num.hosp, 1 + num.se - num.hosp)
    gamma.home <- 1 - gamma.hosp
    gamma.dead <- rbeta(n, 1 + num.dead, 4 + num.hosp - num.dead)
    
    ## Convert mean and variance for a probability to beta parameters
    betaPar <- function(mu, var) {
        alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
        beta <- alpha * (1 / mu - 1)
        list(alpha = alpha, beta = beta)
    }

    betapars <- betaPar(0.45, 0.02)
    lambda.home.rec.TH <- rbeta(n, betapars$alpha, betapars$beta)
    betapars <- betaPar(0.35, 0.02)
    lambda.hosp.rec.TH <- rbeta(n, betapars$alpha, betapars$beta)
    lambda.home.hosp <- gamma.hosp / TH
    lambda.home.home <- (1 - lambda.home.rec.TH)*(1 - lambda.home.hosp)
    lambda.home.rec <- (1 - lambda.home.hosp)*lambda.home.rec.TH
    lambda.hosp.dead <- gamma.dead/TH
    lambda.hosp.hosp <- (1 - lambda.hosp.rec.TH)*(1-lambda.hosp.dead)
    lambda.hosp.rec <- (1 - lambda.hosp.dead)*lambda.hosp.rec.TH 

    ## Convert mean and SD on natural scale to log normal parameters
    lognPar <- function(m,s) {
        s2 <- s^2
        meanlog <- log(m) - 0.5 * log(1 + s2/m^2)
        s2log <- log(1 + (s2/m^2))
        sdlog <- sqrt(s2log)
        list(meanlog = meanlog, sdlog = sdlog)
    }

    lnpars <- lognPar(2300, 90)
    c.home <- rlnorm(n, lnpars$meanlog, lnpars$sdlog)
    lnpars <- lognPar(6500, 980)
    c.hosp <- rlnorm(n, lnpars$meanlog, lnpars$sdlog)
    lnpars <- lognPar(4200, 560)
    c.dead <- rlnorm(n, lnpars$meanlog, lnpars$sdlog)

#### TODO verify second pars are all really variances, not SDs
    betapars <- betaPar(0.98, 0.001)
    e.chemo <- rbeta(n, betapars$alpha, betapars$beta)
    betapars <- betaPar(0.5, 0.02)
    e.home <- rbeta(n, betapars$alpha, betapars$beta)
    betapars <- betaPar(0.2, 0.03)
    e.hosp <- rbeta(n, betapars$alpha, betapars$beta)

    recover.amb <- -log(1 - lambda.home.rec)
    recover.hosp <- -log(1 - lambda.hosp.rec)

    data.frame(
#### variables needed explicitly to evaluate decision model functions below
        SE1, SE2,
        lambda.home.home, lambda.home.hosp, lambda.home.rec,
        lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead, 
        c.home, c.hosp, c.dead,
        e.chemo, e.home, e.hosp,

### other variables that we could calculate EVPPI for
        pi1, pi2, rho,
        gamma.hosp, gamma.dead,
        recover.amb, recover.hosp,
        lambda.home.rec.TH, lambda.hosp.rec.TH)
}

## Calculate average time each patient with adverse events spends in the health states of the Markov model
chemo_markov_model <- function(SE1, SE2,
                               lambda.home.home, lambda.home.hosp, lambda.home.rec,
                               lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead)
    ## SE: predictive distribution of number of patients in 1000 that would experience adverse events
    ## a two-vector containing the number of adverse events for Soc and novel treatment
    
{ 
    ## Markov transition probability matrix 
    ## States: Home care, Hospital care, Recovery, Death
    TH <- 15 
    MM.mat <- matrix(c(lambda.home.home, lambda.home.hosp, lambda.home.rec, 0,
                       0, lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead,
                       0, 0, 1, 0,
                       0, 0, 0, 1),
                     nrow = 4, ncol = 4, byrow = TRUE)
    
    ## Number of patients in each state for each time point 
    ## 3 dimensions: number of states, number of time points, length of side effects
    trace <- array(0, dim = c(4, TH + 1, 2)) 
    trace[1, 1, ] <- c(SE1,SE2)
    
    for(i in 2:(TH + 1)){   
        trace[, i, 1] <- trace[, i - 1, 1] %*% MM.mat 
        trace[, i, 2] <- trace[, i - 1, 2] %*% MM.mat
    }
    
    trace # 4*16*2 array
}


##' Chemotherapy health economic model
##'
##' Calculate expected costs and effects for two treatments as a function of the input parameters.   Each input parameter is a scalar.
##'
##' TODO refer to example of looping to generate PSA output matrix (currently in data-raw)
##'
##' @param SE1 TODO 
##'
##' @param SE2 TODO 
##'
##' @param lambda.home.home Transition probability from 
##'
##' @param lambda.home.hosp Transition probability from 
##'
##' @param lambda.home.rec Transition probability from 
##'
##' @param lambda.hosp.hosp Transition probability from
##'
##' @param lambda.hosp.rec Transition probability from
##' 
##' @param lambda.hosp.dead Transition probability from
##'
##' @param c.home Cost TODO 
##'
##' @param c.hosp Cost TODO 
##'
##' @param c.dead Cost TODO 
##'
##' @param e.chemo Effect TODO 
##
##' @param e.home Effect TODO 
##'
##' @param e.hosp Effect TODO 
##
##' 
##' 
chemo_costeff_fn <- function(SE1, SE2,
                             lambda.home.home, lambda.home.hosp, lambda.home.rec,
                             lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead,
                             c.home, c.hosp, c.dead,
                             e.chemo, e.home, e.hosp){

### auxilary arguments.  TODO make into formal args when we're sure of the format we'll need when calling this function 
    TH <- 15
    N <- 1000

    trace <- chemo_markov_model(SE1, SE2,
                                lambda.home.home, lambda.home.hosp, lambda.home.rec,
                                lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead)
    
    ## costs and effectiveness for four states
    c.states <- c(c.home, c.hosp, 0, 0)
    e.states <- c(e.home, e.hosp, e.chemo, 0)
    
    c.se <- array(NA, dim = 2)

    ## Average cost for both Soc and novel treatment per person per time point
    ## (The cost includes one-off cost of death for patients who died at the end of 15 days)
    c.se[1] <- (sum(c.states %*% trace[,,1]) + c.dead * trace[4, TH + 1, 1])/
      (N * (TH + 1)) 
    c.se[2] <- (sum(c.states %*% trace[,,2]) + c.dead * trace[4, TH + 1, 2])/
      (N * (TH + 1)) 
    c.drug <- c(110, 420)
    cost <- c(c.drug + c.se)

    ## Total QALY of side effects for both Soc and novel treatment
    q.se <- array(NA, dim = 2)
    q.se[1] <- sum(e.states %*% trace[,,1])
    q.se[2] <- sum(e.states %*% trace[,,2])
    ## QALY of total number of patients who do not experience adverse events for 15 days
    q.chemo <- (N - c(SE1,SE2)) * e.chemo * (TH + 1)
    
    ## Average effect for both Soc and novel treatment per person per time point
    effs <- c(q.se + q.chemo)/
      (N * (TH + 1))           

    names(cost) <- paste0("cost",seq_along(cost))
    names(effs) <- paste0("eff",seq_along(effs))
    
    c(cost, effs)
}

chemo_nb_fn <- function(SE1, SE2,
                        lambda.home.home, lambda.home.hosp, lambda.home.rec,
                        lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead,
                        c.home, c.hosp, c.dead,
                        e.chemo, e.home, e.hosp, wtp=20000)
{
    ce <- chemo_costeff_fn(SE1, SE2,
                           lambda.home.home, lambda.home.hosp, lambda.home.rec,
                           lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead,
                           c.home, c.hosp, c.dead,
                           e.chemo, e.home, e.hosp)
    nb <- ce[c("eff1","eff2")]*wtp - ce[c("cost1","cost2")]
    nb
}
