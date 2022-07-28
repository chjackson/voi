evsi_mm <- function(outputs, 
                    inputs,
                    pars,
                    datagen_fn,
                    study, 
                    analysis_fn,
                    analysis_args,
                    model_fn,
                    par_fn,
                    n=100, 
                    Q,
                    npreg_method="gam",
                    verbose, ...){
  model_fn <- check_model_fn(model_fn, par_fn, mfargs=NULL, class(outputs)[1], verbose=verbose)

  ## TODO check analysis_fn has required arguments (data, args, pars)
  analysis_fn <- form_analysis_fn(study, analysis_fn) 
  analysis_args <- form_analysis_args(analysis_args, study, n)
    
  if (length(n) > 1)
    stop("Only one sample size `n` at a time is currently handled in method=`mm`")
  ## Get quantiles of input parameters from PSA sample 
  quants <- mm_gen_quantiles(pars, inputs, Q)
  if (inherits(outputs,"nb")){
    evsi_mm_nb(outputs, inputs=inputs, pars=pars, datagen_fn=datagen_fn, n=n,
               quants=quants, Q=Q, 
               analysis_fn=analysis_fn,
               analysis_args=analysis_args,
               model_fn=model_fn, 
               par_fn=par_fn,
               npreg_method=npreg_method,
               verbose=verbose, ...)
  }
  else if (inherits(outputs,"cea")) {
    evsi_mm_cea(outputs, inputs=inputs, pars=pars, datagen_fn=datagen_fn, n=n,
                quants=quants, Q=Q, 
                analysis_fn=analysis_fn,
                analysis_args=analysis_args,
                model_fn=model_fn, 
                par_fn=par_fn,
                npreg_method=npreg_method,
                verbose=verbose, ...)
  }
}

evsi_mm_nb <- function(outputs, inputs, pars, datagen_fn, n,
                       quants, Q, 
                       analysis_fn,
                       analysis_args,
                       model_fn, 
                       par_fn,
                       npreg_method="gam",
                       verbose=FALSE, ...){
  fits <- evsi_mm_core(nb=outputs, inputs, pars, datagen_fn, n, quants, Q,
                       analysis_fn, analysis_args,
                       model_fn, par_fn, output_row=NULL, npreg_method, verbose, ...)  
  res <- data.frame(
    n = n,
    evsi = calc_evppi(fits$fit_rescaled),
    evppi = calc_evppi(fits$fit)
  )
  if (res$evppi < res$evsi) warning("EVSI > EVPPI may result from approximation error")
  res
}

evsi_mm_cea <- function(outputs, inputs, pars, datagen_fn, n,
                        quants, Q, 
                        analysis_fn,
                        analysis_args,
                        model_fn, 
                        par_fn,
                        npreg_method="gam",
                        verbose=FALSE, ...){
  cfits <- evsi_mm_core(outputs$c, inputs, pars, datagen_fn, n, quants, Q,
                        analysis_fn, analysis_args,
                        model_fn, par_fn,
                        output_row = "c",
                        npreg_method, verbose, ...)
  efits <- evsi_mm_core(outputs$e, inputs, pars, datagen_fn, n, quants, Q,
                        analysis_fn, analysis_args,
                        model_fn, par_fn,
                        output_row = "e",
                        npreg_method, verbose, ...)
  evsi <- calc_evppi_ce(cfits$fit_rescaled, efits$fit_rescaled, outputs$k, verbose=verbose)$evppi
  evppi <- calc_evppi_ce(cfits$fit, efits$fit, outputs$k, verbose=verbose)
  res <- cbind(evppi, n=n, evsi=evsi)
  res
}

evsi_mm_core <- function(nb, # could actually be nb, or c, or e
                         inputs, pars, datagen_fn, n,
                         quants, Q, 
                         analysis_fn,
                         analysis_args,
                         model_fn, 
                         par_fn,
                         output_row=NULL,
                         npreg_method="gam",
                         verbose=FALSE, ...){
  ## Generate future data given Q quantiles
  var_sim <- matrix(nrow=Q, ncol=ncol(nb)-1)
  pb <- progress::progress_bar$new(total = Q) 
  for(i in 1:Q){
    ## Generate one dataset given parameters equal to a specific prior quantile
    simdata <- datagen_fn(inputs = quants[i,,drop=FALSE], n = n, pars=pars)

    ## Fit Bayesian model to future data to get a sample from posterior(pars|simdata)
    postpars <- analysis_fn(simdata, analysis_args, pars)
    niter <- nrow(postpars)
    
    ## Combine with samples from the prior for remaining parameters of the decision model
    modelpars <- par_fn(niter)
    modelpars <- modelpars[names(formals(model_fn))]
    modelpars[names(postpars)] <- postpars
    
    ## Run the decision model, giving a sample from posterior(INB|simdata)
    inbpost <- matrix(nrow=niter, ncol=ncol(nb)-1)
    for (j in 1:niter){
      ## TODO check handles error for nb 1 col or otherwise wrong format 
      nbpost <- do.call(model_fn, modelpars[j,,drop=FALSE])
      if (!is.null(output_row)) nbpost <- nbpost[mfi(nbpost)[[output_row]],]
      inbpost[j,] <- nbpost[-1] - nbpost[1]
    }
    var_sim[i,] <- apply(inbpost, 2, var) # TODO CHECK works with >2 decision options.  TODO should we account for the covariance
    pb$tick()
  }
  mean_prep_var <- apply(var_sim, 2, mean)
  inbprior <- nb[,-1,drop=FALSE] - nb[,1]
  prior_var <- apply(inbprior, 2, var)
  var_prep_mean <- max(0, prior_var - mean_prep_var)  # Why is the 0 needed here? Monte Carlo error?
  
  ## Calculate fitted values for EVPPI
  fit <- fitted_npreg(nb, inputs=inputs, pars=pars, method=npreg_method, verbose=verbose, ...)
  fitn1 <- fit[,-1,drop=FALSE]
  var_fit <- apply(fitn1, 2, var)
  mean_fit <- apply(fitn1, 2, mean)
  
  ## s1/s2 is the prop of variance explained by new data.
  ## if 1 then EVSI=EVPPI.  if 0 then EVSI=0. 
  s1 <- sqrt(var_prep_mean)
  s2 <- sqrt(var_fit)
  fit_rescaled <- (fitn1 - mean_fit) * s1/s2 + mean_fit
  fit_rescaled <- cbind(0, fit_rescaled)

  list(fit=fit, fit_rescaled=fit_rescaled, p_shrink=s1/s2)
}

#' @param pars Character vector of parameters of interest 
#'
#' @param inputs Data frame with sampled values from current distribution of `pars`
#'
#' @param Q Number of equally-spaced quantiles to generate 
#'
#' @return Data frame with one column for each parameter in `pars` and one row per quantile
#' For each variable, the quantiles are randomly permuted.
#'
#' A future version of this function should perhaps use Sobol sequences (randtoolbox package).
#' 
#' @noRd
mm_gen_quantiles <- function(pars,
                             inputs,
                             Q,
                             N.size = NULL){
    quantiles.parameters <- array(NA, dim = c(Q, length(pars)))
    colnames(quantiles.parameters) <- pars
    for(i in 1:length(pars)){
        quantiles.parameters[,i] <- sample(quantile(inputs[,pars[i]],
                                                    probs = 1:Q / (Q + 1), type = 4))
    }
    as.data.frame(quantiles.parameters)
}
