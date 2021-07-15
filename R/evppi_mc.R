##' Traditional two-level Monte Carlo estimator of EVPPI.
##'
##' Traditional two-level Monte Carlo estimator of the expected value of partial
##' perfect information from a decision-analytic model.  Only useful in the
##' simplest of examples.  For realistically complex examples, the methods
##' implemented in the \code{\link{evppi}} function will usually be preferred.
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
##'   set of parameters. This should either return:
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
##' @param nouter Number of outer samples
##'
##' @param ninner Number of inner samples
##'
##' @param wtp Vector of willingness-to-pay values.  Only used if \code{model_fn} is in cost-effectiveness analyis format.
##'
##' @param mfargs Named list of additional arguments to supply to \code{model_fn}. 
##' 
##' @export
evppi_mc <- function(model_fn, par_fn, pars, nouter, ninner,
                     wtp=NULL, mfargs=NULL){
    model_fn <- check_model_fn(model_fn, par_fn, mfargs)
    nopt <- attr(model_fn, "nopt")
    check_parfnn(par_fn, model_fn)
    check_int(nouter, "nouter")
    check_int(ninner, "ninner")
    pars_rep <- par_fn(n=nouter)
    check_evppimc_pars(pars, model_fn, pars_rep)
    pars_rep <- pars_rep[,pars,drop=FALSE]
    rese <- evppimc(model_fn, par_fn, pars, pars_rep, nouter, ninner, nopt, wtp, mfargs) 
    if (inherits(model_fn, "cea")){
        res <- data.frame(wtp=wtp, evppi=rese) 
    } else res <- data.frame(evppi=rese)
    res
}

evppimc <- function(model_fn, ...){
    UseMethod("evppimc", model_fn)
}
    
evppimc.nb <- function(model_fn, par_fn, pars, pars_rep, nouter, ninner, nopt, wtp, mfargs) {
    nb_current <- matrix(nrow=nouter, ncol=nopt)
    nb_ppi <- numeric(nouter)
    nf <- names(formals(model_fn))
    defaults <- get_default_args(model_fn, pars_rep)   
#    registerDoParallel() # TODO 
    pb <- progress::progress_bar$new(total = nouter) 
    for (i in 1:nouter){
        nbi <- matrix(nrow=ninner, ncol=nopt)
        for (j in 1:ninner){
            parsfix <- par_fn(1)
            parsfix[,pars] <- pars_rep[i,]
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

evppimc.cea <- function(model_fn, par_fn, pars, pars_rep, nouter, ninner, nopt, wtp, mfargs) {
    nwtp <- length(wtp)
    if (nwtp < 1)
        stop("If `model_fn` is in cost-effectiveness format, at least one willingness-to-pay should be supplied in `wtp`")
    res <- numeric(nwtp)
    ce_current <- array(dim=c(nouter, 2, nopt))
    cost_ppi <- eff_ppi <- matrix(nrow=nouter, ncol=nopt)
    nf <- names(formals(model_fn))
    defaults <- get_default_args(model_fn, pars_rep)   
    for (i in 1:nouter){
        cei <- array(dim=c(ninner, 2, nopt))
        for (j in 1:ninner){
            parsfix <- par_fn(1)
            parsfix[,pars] <- pars_rep[i,]
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
