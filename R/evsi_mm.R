#' Moment matching method for calculating EVSI
#'
#' @inheritParams evsi
#'
#' @param ... Options passed to nonparametric regression functions
#'
#' @return Data frame with EVSI estimates, as documented in \code{\link{EVSI}}.
#'
#' @noRd
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
                    Q=50,
                    npreg_method="gam",
                    verbose, ...){
  check_pars(pars, inputs, evppi=FALSE)
  model_fn <- check_model_fn(model_fn, par_fn, mfargs=NULL, outputs, verbose=verbose)
  check_pars_in_modelfn(pars, model_fn)
  analysis_args <- form_analysis_args(analysis_args, study, n)
  analysis_fn <- form_analysis_fn(study, analysis_fn, analysis_args, datagen_fn, inputs, n, pars, pars_datagen) 
    
  if (inherits(outputs,"nb")){
    evsi_mm_nb(outputs, inputs=inputs,
               pars=pars, pars_datagen=pars_datagen,
               datagen_fn=datagen_fn, n=n,
               Q=Q, 
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
                Q=Q, 
                analysis_fn=analysis_fn,
                analysis_args=analysis_args,
                model_fn=model_fn, 
                par_fn=par_fn,
                npreg_method=npreg_method,
                verbose=verbose, ...)
  }
}

#' Moment matching method for calculating EVSI
#'
#' @inheritParams evsi
#'
#' Specific to `outputs` of net benefit form
#'
#' @noRd 
evsi_mm_nb <- function(outputs, inputs, pars, pars_datagen, datagen_fn, n,
                       Q, 
                       analysis_fn,
                       analysis_args,
                       model_fn, 
                       par_fn,
                       npreg_method="gam",
                       verbose=FALSE, ...){
  fits <- evsi_mm_core(nb=outputs, inputs, pars, pars_datagen, datagen_fn, n, Q,
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

#' Moment matching method for calculating EVSI
#'
#' @inheritParams evsi
#'
#' Specific to `outputs` of "cost-effectiveness analysis" form
#'
#' @noRd 
evsi_mm_cea <- function(outputs, inputs, pars,
                        pars_datagen,
                        datagen_fn, n,
                        Q, 
                        analysis_fn,
                        analysis_args,
                        model_fn, 
                        par_fn,
                        npreg_method="gam",
                        verbose=FALSE, ...){
  cfits <- evsi_mm_core(outputs$c, inputs, pars, pars_datagen, datagen_fn, n, Q,
                        analysis_fn, analysis_args,
                        model_fn, par_fn,
                        output_row = "c",
                        npreg_method, verbose, ...)
  efits <- evsi_mm_core(outputs$e, inputs, pars, pars_datagen, datagen_fn, n, Q,
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

#' Moment-matching method for EVSI calculation - all the hard work is
#' done in this core function.
#'
#' @inheritParams evsi
#' 
#' @param nb Matrix of outputs, with columns indicating decision
#'   options, and rows indicating samples. Could be either net
#'   benefits, costs or effects.
#'
#' @param output_row Name of the row of the decision model output
#' that is used.  Only used if `nb` is costs or effects.
#'
#' @return A list with components:
#'
#' `fit` Fitted values for EVPPI calculation.  Matrix with `nsim` rows
#' (number of samples from the parameter uncertainty distribution) and
#' `d` columns (number of decision options).  "Fitted values" means
#' the expected net benefit given further information, conditionally
#' on specific parameter values.
#' 
#' `fits` Array of fitted values for EVSI calculation, with dimensions
#' `nsim`  x  `d`  x  `length(n)`, where `length(n)` is the number of
#' distinct sample sizes that the EVSI calculation was requested for.
#'
#' `p_shrink` is the ratio of standard deviations, describing the
#' proportion of uncertainty explained by the further information
#' gained from the proposed study.
#'
#' @noRd
evsi_mm_core <- function(nb, # could actually be nb, or c, or e
                         inputs, 
                         pars, 
                         pars_datagen, 
                         datagen_fn, 
                         n,
                         Q, 
                         analysis_fn,
                         analysis_args,
                         model_fn, 
                         par_fn,
                         output_row=NULL,
                         npreg_method="gam",
                         verbose=FALSE, ...){
  ## Determine grid of Q parameters from quantiles
  quants <- mm_gen_quantiles(pars_datagen, inputs, Q)
  ## Generate future data given these parameters
  ncomp <- ncol(nb) - 1 # number of comparisons, not number of decision options
  ncov <- ncomp*(ncomp+1)/2
  var_sim <- matrix(nrow=Q, ncol=ncov)
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
    var_sim[i,] <- covvec(inbpost)

    pb$tick()
  }
  
  mean_prep_var <- apply(var_sim, 2, mean)
  inbprior <- nb[,-1,drop=FALSE] - nb[,1]
  prior_var <- covvec(inbprior)
  
  if (verbose) message("Calculating fitted values for EVPPI...")
  fit <- fitted_npreg(nb, inputs=inputs, pars=pars, method=npreg_method, verbose=verbose, ...)
  fitn1 <- fit[,-1,drop=FALSE]
  var_fit <- covvec(fitn1)
  mean_fit <- apply(fitn1, 2, mean)

  var_prep_mean <- matrix(nrow=length(n), ncol=ncov)
  if (length(n) > 1){
    ## do regression to relate variance reduction to sample size 
    quants <- c(0.025, 0.125, 0.5, 0.875, 0.975)
    var_red_n <- matrix(nrow=length(n), ncol=5, dimnames=list(NULL, quants))
    if (verbose) message("Running regression on sample size...")
    for (d in 1:ncov){
      var_red <- prior_var[d] - var_sim[,d] 
      if (sd(var_red)==0)
        beta <- 0
      else 
        beta <- regression_on_sample_size(var_red, nfit, var_fit[d])
      for (j in 1:length(n)){
        var_red_n[j,] <- quantile(var_fit[d] * n[j] / (n[j] + beta), quants)
        var_prep_mean[j,d] <- var_red_n[j,"0.5"]
      }
    }
  } else {
    var_prep_mean[1,] <- pmax(0, prior_var - mean_prep_var) # Correct for Monte Carlo error
  }
  
  fit_rescaled <- array(dim = c(dim(fit), length(n)))
  if (ncomp==1){
    p_shrink <- numeric(length(n))
    s2 <- sqrt(var_fit)
  } else {
    p_shrink <- vector(length(n), mode="list")
    s2inv <- MatrixSqrt(cov(fitn1), inverse=TRUE)
  }
  for (j in 1:length(n)){  
    fit_rescaled[,1,j] <- 0 
    if (ncomp == 1){
      s1 <- sqrt(var_prep_mean[j,])
      p_shrink[j] <- s1/s2 # prop of unc explained by new data. if 1 then EVSI=EVPPI, if 0 then EVSI=0.
      fit_rescaled[,2,j] <- (fitn1[,1] - mean_fit[1]) * p_shrink[j] + mean_fit[1]
    }
    else {
      s1 <- MatrixSqrt(covvec2mat(var_prep_mean[j,]))
      p_shrink[[j]] <- s2inv %*% s1
      mean_mat <- matrix(mean_fit, nrow=nrow(fitn1), ncol = ncomp, byrow = TRUE)
      fit_rescaled[,-1,j] <- (fitn1 - mean_mat) %*% p_shrink[[j]] + mean_mat
    }  
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


#' Bayesian nonlinear regression of the variance reduction in terms of
#' the proposed study sample size
#'
#' @return A sample from the posterior distribution of the parameter
#' `beta` governing this dependence.
#' 
#' @noRd
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


## Covariances as a vector rather than a matrix
covvec <- function(x){ 
  covmat2vec(cov(x))
}

## Convert a covariance matrix to a vector (excluding the above-diagonals)
covmat2vec <- function(x){
  c(diag(x), x[lower.tri(x)])
}

## Convert a vectorised covariance matrix back to a matrix
covvec2mat <- function(x){
  n <- trunc(sqrt(length(x)*2))
  mat <- diag(x[1:n])
  covs <- x[(n+1):length(x)]
  mat[lower.tri(mat)] <- covs
  mat[upper.tri(mat)] <- t(mat)[upper.tri(t(mat))]
  mat
}

#' Matrix square root.  Used to implement (co)variance reduction when
#' using the moment matching method for models with two or more
#' treatment comparisons (i.e. three or more decision options).
#'
#' @noRd
MatrixSqrt <- function(x, inverse=FALSE){
  e <- eigen(x)
  if (any(e$values<0)){ 
    e <- eigen(Matrix::nearPD(x)$mat)
  }
  sq <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
  if (inverse) chol2inv(chol(sq)) else sq
}
