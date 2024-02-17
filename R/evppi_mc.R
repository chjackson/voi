##' Traditional two-level Monte Carlo estimator of EVPPI.
##'
##' Traditional two-level Monte Carlo estimator of the expected value of partial
##' perfect information from a decision-analytic model.  Only useful in the
##' simplest of examples.  For realistically complex examples, the methods
##' implemented in the \code{\link{evppi}} function, based on regression,
##' will usually be much more computationally efficient.
##'
##' See the \href{https://chjackson.github.io/voi/articles/voi.html#evppimc}{package overview / Get Started vignette} for an example of using this function. 
##'
##'
##' @param pars A character vector giving the parameters of interest, for which
##'   the EVPPI is required.   This should correspond to an explicit argument to
##'   \code{model_fn}.
##'
##'   The parameters of interest are assumed to have uncertainty distributions
##'   that are independent of those of the other parameters.
##'
##' @param model_fn A function to evaluate a decision-analytic model at a given
##'   set of parameters. This should have one argument per parameter, and return either:
##'
##'   (net benefit format) a vector giving the net benefit for each decision
##'   option, or
##'
##'   (cost-effectiveness analysis format) a matrix or data frame with two rows,
##'   and one column for each decision option.  If the rows have names
##'   \code{"e"} and \code{"c"} then these are assumed to be the effects and
##'   costs respectively.
##'
##'   Otherwise, the first row is assumed to be the effects, and the second the
##'   costs.
##'
##' @param par_fn  A function to generate a random sample of values for the
##'   parameters of \code{model_fn}. This should return a matrix or a data frame
##'   with named columns matching the arguments of \code{model_fn}.
##'
##'   If any required arguments to \code{model_fn} are not supplied in this
##'   return value, then \code{evppi_mc} looks for them in the list supplied as
##'   the \code{mfargs} argument.
##'
##'   If any required arguments are not found in the results of \code{par_fn} or
##'   \code{mfargs}, and if \code{model_fn} defines default values for those
##'   arguments, then those default values are used.
##'
##'   The first argument of \code{par_fn} should be an integer \code{n} denoting
##'   the number of random values to draw for each parameter.  The object
##'   returned by \code{par_fn} should then have \code{n} rows, and one column
##'   for each parameter. If one value is drawn, then \code{par_fn} is also
##'   allowed to return a vector, but this should still be named.
##'
##'   The parameters may be correlated.  If we wish to compute the EVPPI for a
##'   parameter which is correlated with a different parameter q, then `par_fn`
##'   must have an argument with the name of that parameter.  If that argument
##'   is set to a fixed value, then `par_fn` should return a sample drawn
##'   conditionally on that value.  If that argument is not supplied, then
##'   `par_fn` must return a sample drawn from the marginal distribution. See
##'   the vignette for an example.
##'
##' @param nouter Number of outer samples
##'
##' @param ninner Number of inner samples
##'
##' @param k Vector of willingness-to-pay values.  Only used if
##'   \code{model_fn} is in cost-effectiveness analyis format.
##'
##' @param mfargs Named list of additional arguments to supply to
##'   \code{model_fn}.
##'
##' @param verbose Set to \code{TRUE} to print some additional messages to
##' help with debugging.
##' 
##' @return A data frame with a column \code{pars}, indicating the parameter(s),
##'   and a column \code{evppi}, giving the corresponding EVPPI. 
##'
##'   If \code{outputs} is of "cost-effectiveness analysis" form, so that there is
##'   one EVPPI per willingness-to-pay value, then a column \code{k} identifies the 
##'   willingness-to-pay.
##'
##' @export
evppi_mc <- function(model_fn, par_fn, pars, nouter, ninner,
                     k=NULL, mfargs=NULL, verbose=FALSE){
    model_fn <- check_model_fn(model_fn, par_fn, mfargs, verbose=verbose)
    nopt <- attr(model_fn, "nopt")
    check_parfnn(par_fn, model_fn)
    check_int(nouter, "nouter")
    check_int(ninner, "ninner")
    pars_rep <- par_fn(n=nouter)
    check_evppimc_pars(pars, model_fn, pars_rep)
    pars_rep <- pars_rep[,pars,drop=FALSE]
    rese <- evppimc(model_fn=model_fn, par_fn=par_fn, pars=pars, pars_rep=pars_rep,
                    nouter=nouter, ninner=ninner, nopt=nopt, wtp=k, mfargs=mfargs) 
    if (inherits(model_fn, "cea")){
        res <- data.frame(k=k, evppi=rese) 
    } else res <- data.frame(evppi=rese)
    res
}

evppimc <- function(model_fn, ...){
    UseMethod("evppimc", model_fn)
}
    
##' @noRd
evppimc.nb <- function(model_fn, par_fn, pars, pars_rep, nouter, ninner, nopt, mfargs, ...) {
    nb_current <- matrix(nrow=nouter, ncol=nopt)
    nb_ppi <- numeric(nouter)
    nf <- names(formals(model_fn))
    defaults <- get_default_args(model_fn, pars_rep)   
    pars_corr <- intersect(names(formals(par_fn)), nf)
    pb <- progress::progress_bar$new(total = nouter) 
    for (i in 1:nouter){
        nbi <- matrix(nrow=ninner, ncol=nopt)
        for (j in 1:ninner){
            parsfix <- sample_conditional(par_fn, pars, pars_corr, vals=pars_rep[i,,drop=FALSE])
            args <- c(parsfix, mfargs, defaults)[nf]
            nbi[j,] <- do.call(model_fn, args)
        }
        pb$tick()
        nb_ppi[i] <- max(colMeans(nbi))
        args <- c(par_fn(1), mfargs, defaults)[nf]
        nb_current[i,] <- do.call(model_fn, args)
    }
    mean(nb_ppi) - max(colMeans(nb_current))
}

##' @noRd
evppimc.cea <- function(model_fn, par_fn, pars, pars_rep, nouter, ninner, nopt, wtp, mfargs, ...) {
    nwtp <- length(wtp)
    if (nwtp < 1)
        stop("If `model_fn` is in cost-effectiveness format, at least one willingness-to-pay should be supplied in `k`")
    res <- numeric(nwtp)
    ce_current <- array(dim=c(nouter, 2, nopt))
    cost_ppi <- eff_ppi <- matrix(nrow=nouter, ncol=nopt)
    nf <- names(formals(model_fn))
    defaults <- get_default_args(model_fn, pars_rep)   
    pars_corr <- intersect(names(formals(par_fn)), nf)
    for (i in 1:nouter){
        cei <- array(dim=c(ninner, 2, nopt))
        for (j in 1:ninner){
            parsfix <- sample_conditional(par_fn, pars, pars_corr, vals=pars_rep[i,,drop=FALSE])
            args <- c(parsfix, mfargs, defaults)[nf]
            resj <- do.call(model_fn, args)
            cei[j,,] <- resj
        }
        inds <- mfi(res)
        eff_ppi[i,] <- colMeans(cei[,inds$c,])
        cost_ppi[i,] <- colMeans(cei[,inds$e,])
        args <- c(par_fn(1), mfargs, defaults)[nf]
        ce_current[i,,] <- do.call(model_fn, args)
    }
    for (k in 1:nwtp){
        nb_ppi <- apply(eff_ppi*wtp[k] - cost_ppi, 1, max)
        nb_current <- ce_current[,inds$e,]*wtp[k] - ce_current[,inds$c,]
        res[k] <- mean(nb_ppi) - max(colMeans(nb_current))
    }
    res
}

## Sample from the joint distribution defined by `par_fn`, given fixed values
## `vals` for parameters `pars`.   If `par_fn` has any arguments other than `n`
## (named in `pars_corr`), these are assumed to be parameters that are
## correlated with other parameters, so that fixing their value will change the
## conditional distribution of the remaining parameters.  The values `vals` of
## the parameters fixed in the inner EVPPI loop are then supplied for these
## arguments, to allow par_fn to compute the appropriate conditional
## distribution

sample_conditional <- function(par_fn, pars, pars_corr, vals) {
    args_fixed <- if (length(pars_corr) == 0) NULL else as.list(vals[,pars_corr])
    parsfix <- do.call("par_fn", c(list(n=1), args_fixed))
    parsfix[,pars] <- vals
    parsfix
}

check_int <- function(n, name){
    if (!is.numeric(n) || (n < 2))
        stop(sprintf("%s is `%s`, this should be a number greater than 1", name, n))
}

check_evppimc_pars <- function(pars, model_fn, pars_rep){
    badpars <- setdiff(pars, names(formals(model_fn)))
    if (length(badpars)>0){
        stop(sprintf("parameters of interest `%s` not found in arguments of `model_fn`", paste(badpars,collapse=",")))
    }    
    badpars <- setdiff(pars, colnames(pars_rep))
    if (length(badpars)>0){
        stop(sprintf("parameters of interest `%s` not found in columns of object returned by `par_fn`",
                     paste(badpars,collapse=",")))
    }
}
