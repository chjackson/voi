evsi_mm <- function(outputs, 
                    inputs,
                    pars,
                    pars_datagen,
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
  check_pars_in_modelfn(pars, model_fn)
  analysis_args <- form_analysis_args(analysis_args, study, n)
  analysis_fn <- form_analysis_fn(study, analysis_fn, analysis_args, datagen_fn, inputs, n, pars, pars_datagen) 
    
  ## Get quantiles of input parameters from PSA sample 
  quants <- mm_gen_quantiles(pars_datagen, inputs, Q)
  if (inherits(outputs,"nb")){
    evsi_mm_nb(outputs, inputs=inputs,
               pars=pars, pars_datagen=pars_datagen,
               datagen_fn=datagen_fn, n=n,
               quants=quants, Q=Q, 
               analysis_fn=analysis_fn,
               analysis_args=analysis_args,
               model_fn=model_fn, 
               par_fn=par_fn,
               npreg_method=npreg_method,
               verbose=verbose, ...)
  }
  else if (inherits(outputs,"cea")) {
    evsi_mm_cea(outputs, inputs=inputs,
                pars=pars, pars_datagen=pars_datagen,
                datagen_fn=datagen_fn, n=n,
                quants=quants, Q=Q, 
                analysis_fn=analysis_fn,
                analysis_args=analysis_args,
                model_fn=model_fn, 
                par_fn=par_fn,
                npreg_method=npreg_method,
                verbose=verbose, ...)
  }
}

evsi_mm_nb <- function(outputs, inputs, pars, pars_datagen, datagen_fn, n,
                       quants, Q, 
                       analysis_fn,
                       analysis_args,
                       model_fn, 
                       par_fn,
                       npreg_method="gam",
                       verbose=FALSE, ...){
  fits <- evsi_mm_core(nb=outputs, inputs, pars, pars_datagen, datagen_fn, n, quants, Q,
                       analysis_fn, analysis_args,
                       model_fn, par_fn, output_row=NULL, npreg_method, verbose, ...)

  evsis <- evppis <- numeric(length(n))
  for (i in 1:length(n)){
    evsis[i] <- calc_evppi(fits$fit_rescaled[,,i])
  }
  res <- data.frame(n=n, evsi=evsis)
  attr(res, "evppi") <- data.frame(evppi = calc_evppi(fits$fit))
  if (any(attr(res,"evppi")$evppi < res$evsi))
    message("EVSI > EVPPI may result from approximation error")
  res
}

evsi_mm_cea <- function(outputs, inputs, pars,
                        pars_datagen,
                        datagen_fn, n,
                        quants, Q, 
                        analysis_fn,
                        analysis_args,
                        model_fn, 
                        par_fn,
                        npreg_method="gam",
                        verbose=FALSE, ...){
  cfits <- evsi_mm_core(outputs$c, inputs, pars, pars_datagen, datagen_fn, n, quants, Q,
                        analysis_fn, analysis_args,
                        model_fn, par_fn,
                        output_row = "c",
                        npreg_method, verbose, ...)
  efits <- evsi_mm_core(outputs$e, inputs, pars, pars_datagen, datagen_fn, n, quants, Q,
                        analysis_fn, analysis_args,
                        model_fn, par_fn,
                        output_row = "e",
                        npreg_method, verbose, ...)

  evsis <- evppis <- vector(length(n), mode="list")
  for (i in 1:length(n)){
    evsis[[i]] <- calc_evppi_ce(cfits$fit_rescaled[,,i], efits$fit_rescaled[,,i],
                                outputs$k, verbose=verbose)$evppi
    evsis[[i]] <- cbind(n=n[i], k=outputs$k, evsi=evsis[[i]])
  }
  res <- as.data.frame(do.call("rbind", evsis))
  attr(res, "evppi") <- calc_evppi_ce(cfits$fit, efits$fit, outputs$k, verbose=verbose)
  res
}

evsi_mm_core <- function(nb, # could actually be nb, or c, or e
                         inputs, 
                         pars, 
                         pars_datagen, 
                         datagen_fn, 
                         n,
                         quants, Q, 
                         analysis_fn,
                         analysis_args,
                         model_fn, 
                         par_fn,
                         output_row=NULL,
                         npreg_method="gam",
                         verbose=FALSE, ...){
  ## Generate future data given Q quantiles
  ncomp <- ncol(nb) - 1 # number of comparisons, not number of decision options
  var_sim <- matrix(nrow=Q, ncol=ncomp)
  pb <- progress::progress_bar$new(total = Q) 
  if (length(n) > 1) {
    nfit <- unique(round(seq(sqrt(min(n)), sqrt(max(n)), length=Q)^2))
  } else nfit <- rep(n, Q)

  for(i in 1:Q){
    ## Generate one dataset given parameters equal to a specific prior quantile
    simdata <- datagen_fn(inputs = quants[i,,drop=FALSE], n = nfit[i], pars=pars_datagen)
    
    ## Fit Bayesian model to future data to get a sample from posterior(pars|simdata)
    analysis_args$n <- nfit[i]
    if (verbose) message("Running Bayesian analysis function...")
    postpars <- analysis_fn(simdata, analysis_args, pars)
    niter <- nrow(postpars)
    if (i==1)
      priorpars <- par_fn(niter)
    
    ## Combine with samples from the prior for remaining parameters of the decision model
    ## (newly-generated samples, in case number of posterior samples (niter) desired is more than nrow(inputs))
    modelpars <- priorpars[names(formals(model_fn))]
    modelpars[names(postpars)] <- postpars
    
    ## Run the decision model, giving a sample from posterior(INB|simdata)
    inbpost <- matrix(nrow=niter, ncol=ncomp)
    if (verbose) message("Evaluating decision model at updated parameters...")
    for (j in 1:niter){
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
  
  if (verbose) message("Calculating fitted values for EVPPI...")
  fit <- fitted_npreg(nb, inputs=inputs, pars=pars, method=npreg_method, verbose=verbose, ...)
  fitn1 <- fit[,-1,drop=FALSE]
  var_fit <- apply(fitn1, 2, var)
  mean_fit <- apply(fitn1, 2, mean)

  var_prep_mean <- matrix(nrow=length(n), ncol=ncomp)
  if (length(n) > 1){
    ## do regression to relate variance reduction to sample size 
    quants <- c(0.025, 0.125, 0.5, 0.875, 0.975)
    var_red_n <- matrix(nrow=length(n), ncol=5, dimnames=list(NULL, quants))
    if (verbose) message("Running regression on sample size...")
    for (d in 1:ncomp){
      var_red <- prior_var[d] - var_sim[,d] 
      if (sd(var_red)==0)
        beta <- 0
      else 
        beta <- regression_on_sample_size(var_red, nfit, var_fit[d])
      for (j in 1:length(n)){
        var_red_n[j,] <- quantile(var_fit[d] * n[j] / (n[j] + beta), quants)
        var_prep_mean[j,d] <- var_red_n[j,"0.5"] # TODO return uncertainty in these
      }
    }
  } else {
    var_prep_mean[1,] <- pmax(0, prior_var - mean_prep_var)  # Why is the 0 needed here? Monte Carlo error?
  }
  
  fit_rescaled <- array(dim = c(dim(fit), length(n)))
  p_shrink <- numeric(length(n))
  for (j in 1:length(n)){  
    s1 <- sqrt(var_prep_mean[j,])
    s2 <- sqrt(var_fit)
    p_shrink[j] <- s1/s2 # prop of unc explained by new data. if 1 then EVSI=EVPPI, if 0 then EVSI=0. 
    fit_rescaled[,,j] <- cbind(0,
    (fitn1 - mean_fit) * p_shrink[j] + mean_fit
    )
  }
  
  list(fit=fit, fit_rescaled=fit_rescaled, p_shrink=p_shrink)
}

#' @param pars Character vector of parameters required to generate the study data
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
    quants <- array(NA, dim = c(Q, length(pars)))
    colnames(quants) <- pars
    for(i in 1:length(pars)){
        quants[,i] <- sample(quantile(inputs[,pars[i]],
                                                    probs = 1:Q / (Q + 1), type = 4))
    }
    as.data.frame(quants)
}


regression_on_sample_size <- function(var_reduction,
                                      sample_size,
                                      var_base) {
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("JAGS and the R package \"rjags\" must be installed to use this function.", call. = FALSE)
  }
  dat <- list(
    N = length(sample_size),
    y = as.vector(var_reduction),
    x = as.vector(sample_size),
    sigma_mu = sd(var_reduction)/2,
    sigma_tau = 1/sd(var_reduction),
    Nmax_par = max(sample_size)/2,
    shape_Nmax = 0.05 / max(sample_size),
    var_base = var_base
  )
  ini <- list(sigma=sd(var_reduction)/2, beta=1)
  mod <- "model {
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- var_base * (x[i]/(x[i] + beta))
    }
    tau <- 1/(sigma*sigma)
    sigma ~ dt(sigma_mu, sigma_tau, 3)I(0, )
    beta ~ dnorm(Nmax_par, shape_Nmax)I(0, )
  }
  "
  jagsmod <- rjags::jags.model(textConnection(mod), data=dat, inits=ini, quiet=TRUE)
  update(jagsmod, 1000, progress.bar="none")
  sam <- rjags::coda.samples(jagsmod, variable.names="beta", n.iter=3000, progress.bar="none")
  beta <- as.data.frame(sam[[1]])[,1]
  beta
}
