evsi_mm <- function(outputs, 
                    inputs,
                    output_type,
                    pars,
                    datagen_fn,
                    analysis_model,
                    analysis_options,
                    decision_model,
                    n=100, likelihood,
                    Q,
                    npreg_method="gam",
                    verbose, ...){
    ## might be neater to have just one function that does both nb and cea formats 
    if (output_type=="nb"){
        evsi_mm_nb(nb=outputs, inputs=inputs, pars=pars, datagen_fn=datagen_fn, n=n,
                   Q=Q, 
                   analysis_model=analysis_model,
                   analysis_options=analysis_options,
                   decision_model=decision_model, 
                   npreg_method=npreg_method,
                   verbose=verbose, ...)
    }
    else if (output_type=="cea"){
        evsi_mm_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                    inputs=inputs, pars=pars, datagen_fn=datagen_fn, n=n,
                    Q=Q, 
                    analysis_model=analysis_model,
                    analysis_options=analysis_options,
                    decision_model=decision_model, 
                    npreg_method=npreg_method,
                    verbose=verbose, ...)
    }
}

evsi_mm_nb <- function(nb, inputs, pars, datagen_fn, n, Q=Q, 
                       analysis_model,
                       analysis_options,
                       decision_model, 
                       npreg_method,
                       verbose, ...){

    ## TODO adapt from EVSI package

    ## Get quantiles of input parameters from PSA sample 
    ## TODO multiple sample sizes argument to this 
    quants <- mm_gen_quantiles(pars, inputs, Q)

    ## Instead of a monolithic "mm.post.var" function
    ## break it down into logical steps
    ## generate data / fit bayesian model / run HE model 

    ## Generate future data given Q quantiles
    var_sim <- vector(Q, mode="list")
    for(i in 1:Q){
        simdata <- datagen_fn(quants[i,], quants[i,"N"])

        ## Fit Bayesian model to future data to get posterior
        ## Skeleton example of an "analysis_model" function is analysis_model_jags below. 
        postpars <- analysis_model(simdata, analysis_options)

        ## TODO append simulations from the prior for all parameters in the decision model
        ## Users could specify these in analysis_model through extra pars in their JAGS model that are simulated from their priors
        ## Any not supplied there might be filled in through an extra user-supplied function to simulate parameters from the model - a built-in example is chemo_prior_pars 
        
        ## Run HE model again with posterior sample (TODO get exact arg formats correct )
        inb <- decision_model(postpars)
        var_sim[[q]] <- var(inb)
    }

    ## Calculate EVPPI
    ## version that works on net benefit 
    evppi_fit <- fitted_npreg(nb, inputs=inputs, pars=pars, method=npreg_method, verbose=verbose, ...)

    ## version that splits costs, effects and WTP 
#    evppi_cfit <- fitted_npreg(costs, inputs=inputs, pars=pars, method=npreg_method, verbose=verbose, ...)
#    evppi_efit <- fitted_npreg(effects, inputs=inputs, pars=pars, method=npreg_method, verbose=verbose, ...)

    ## Now rewrite evsi.calc.   Instead of accepting a monolithic mm.var object,
    ## break that object into logical components that might come from inputs to evsi_mm_nb: 
    ## * simulated variances
    ## * EVPPIs 
    ## * original simulations from decision model 
    ## * WTPs [ will need to handle in evsi_mm_cea method, but not in evsi_mm_nb ] 
    ## * sample sizes 

    evsi <- evsi.calc(var_sim,    # simulated variances 
                      evppi_res,  # EVPPI results  [ todo handle nb and cea formats ] 
                      nb,         # 
                      wtp=NULL,
                      N=NULL,
                      CI=NULL)

    ## Perhaps split into two functions, one to calculate EVSI the standard moment-matching way, and a second function to interpolate EVSIs for different sample sizes by regression 
    ## For multiple sample sizes, I reckon to just use our favourite MCMC updater to do the nonlinear regression.
    ## It's a different situation from when the user specifies the Bayesian model governing their proposed study.
    ## Here the Bayesian model is pre-specified, so we can pre-specify our favourite software to fit it [ though will it always follow the same functional shape ? ]

    ## Return object of similar format so we can do same output analysis / graph things as EVSI package
}

## Example of "analysis_model" function that user could supply 
## Main thing user needs to supply is a textfile of JAGS model code
## All data needed by this model should be returned by the user-supplied datagen_fn. 
## User would call evsi(..., analysis_model = analysis_model_jags, analysis_options = list(jags_model_file = "my_model.txt", nburn=1000, niter=10000))

analysis_model_jags <- function(simdata, jags_options){
    if (is.null(jags_options$model_file)) stop("JAGS model file must be specified in `analysis_options`") 
    if (is.null(jags_options$nburn)) nburn <- 1000
    if (is.null(jags_options$niter)) niter <- 10000  # or such like 
    ## then convert simdata from data frame to JAGS format

    ## build inits

    ## run JAGS

    ## extract sample from posterior, return as data frame 
}


## This may have some code in common with evsi_mm_nb
## Let's see how it works out
## Or just have one function that does both if it turns out tidier. 

evsi_mm_cea <- function(costs, effects, wtp, inputs, pars, datagen_fn, n, Q=Q, 
                    analysis_model,
                    decision_model, 
                    npreg_method,
                    verbose, ...){
}


mm_gen_quantiles <- function(pars, inputs, Q, N.size = NULL){
    quantiles.parameters <- array(NA, dim = c(Q, length(pars)))
    colnames(quantiles.parameters) <- pars
    for(i in 1:length(pars)){
        quantiles.parameters[,i] <- sample(quantile(inputs[,pars[i]],
                                                    probs = 1:Q / (Q + 1), type = 4))
    }

# TODO handle multiple sample sizes    
#    if (!is.null(N.size)) {
#        N.size <- round(exp(seq(log(min(N.size)),log(max(N.size)),length.out = Q)))
#        quantiles.aug <- cbind(quantiles.parameters, N = N.size)
#        quantiles.parameters <- quantiles.aug
#    }
    as.data.frame(quantiles.parameters)
}
