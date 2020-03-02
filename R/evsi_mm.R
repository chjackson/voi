evsi_mm <- function(outputs, 
                    inputs,
                    output_type,
                    poi,
                    datagen_fn,
                    analysis_model,
                    decision_model,
                    n=100, likelihood,
                    Q,
                    npreg_method="gam",
                    verbose, ...){
    ## might be neater to have just one function that does both nb and cea formats 
    if (output_type=="nb"){
        evsi_mm_nb(nb=outputs, inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                   Q=Q, 
                   analysis_model=analysis_model,
                   decision_model=decision_model, 
                   npreg_method=npreg_method,
                   verbose=verbose, ...)
    }
    else if (output_type=="cea"){
        evsi_mm_cea(costs=outputs$c, effects=outputs$e, wtp=outputs$k,
                    inputs=inputs, poi=poi, datagen_fn=datagen_fn, n=n,
                    Q=Q, 
                    analysis_model=analysis_model,
                    decision_model=decision_model, 
                    npreg_method=npreg_method,
                    verbose=verbose, ...)
    }
}

evsi_mm_nb <- function(nb, inputs, poi, datagen_fn, n, Q=Q, 
                       analysis_model,
                       decision_model, 
                       npreg_method,
                       verbose, ...){

    ## TODO adapt from EVSI package

    ## Get quantiles of input parameters from PSA sample 
    ## TODO multiple sample sizes argument to this 
    quants <- mm_gen_quantiles(poi, inputs, Q)

    ## Instead of a monolithic "mm.post.var" function
    ## break it down into logical steps
    ## generate data / fit bayesian model / run HE model 

    ## Generate future data given Q quantiles
    var_sim <- vector(Q, mode="list")
    for(i in 1:Q){
        simdata <- datagen_fn(quants[i,], quants[i,"N"])

        ## Fit Bayesian model to future data to get posterior
        ## TODO return data frame of simulations to be passed to model function 
        postpars <- analysis_model(simdata, analysis_options)

        ## TODO append simulations from the prior (decide most sensible place to do this, could be within analysis_model, or outside)

        ## Run HE model again with posterior sample (TODO get exact arg formats correct )
        inb <- decision_model(postpars)
        var_sim[[q]] <- var(inb)
    }

    ## Calculate EVPPI
    ## version that works on net benefit 
    evppi_fit <- fitted_npreg(nb, inputs=inputs, poi=poi, method=npreg_method, verbose=verbose, ...)

    ## version that splits costs, effects and WTP 
    evppi_cfit <- fitted_npreg(costs, inputs=inputs, poi=poi, method=npreg_method, verbose=verbose, ...)
    evppi_efit <- fitted_npreg(effects, inputs=inputs, poi=poi, method=npreg_method, verbose=verbose, ...)

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

    ## For multiple sample sizes, I reckon to just use our favourite MCMC updater to do the nonlinear regression.
    ## It's a different situation from when the user specifies the Bayesian model governing their proposed study.
    ## Here the Bayesian model is pre-specified, so we can pre-specify our favourite software to fit it 

    ## Return object of similar format so we can do same output analysis / graph things
}


## This may have some code in common with evsi_mm_nb
## Let's see how it works out
## Or just have one function that does both if it turns out tidier. 

evsi_mm_cea <- function(costs, effects, wtp, inputs, poi, datagen_fn, n, Q=Q, 
                    analysis_model,
                    decision_model, 
                    npreg_method,
                    verbose, ...){
}


mm_gen_quantiles <- function(poi, inputs, Q, N.size = NULL){
    quantiles.parameters <- array(NA, dim = c(Q, length(poi)))
    colnames(quantiles.parameters) <- poi
    for(i in 1:length(poi)){
        quantiles.parameters[,i] <- sample(quantile(param.mat[,poi[i]],
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
