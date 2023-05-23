##' Expected net benefit of sampling
##'
##' Calculates the expected net benefit of sampling for a typical study to inform
##' a health economic evaluation, given estimates of the per-person expected
##' value of sample information, decision population size and study setup and
##' per-participant costs.   The optimal sample size for each willingness-to-pay,
##' population size and time horizon is also determined.
##'
##' \code{pop},\code{time} and \code{dis} may be supplied as vectors
##' of different lengths.  In that case, the ENBS is calculated for all
##' possible combinations of the values in these vectors. 
##'
##' @param evsi Data frame giving estimates of the expected value of sample
##'   information, as returned by \code{\link{evsi}}.  This may contain
##'   multiple estimates, one for each sample size and willingness to pay.
##'
##' @param costs_setup Setup costs of the study.  This can either be a constant,
##'   or a vector of two elements giving a 95% credible interval (with mean
##'   defined by the midpoint), or a vector of three elements assumed to define
##'   the mean and 95% credible interval.
##'
##' @param costs_pp Per-participant costs of the study, supplied in the same
##'   format as \code{cost_setup}.
##'
##' @param pcut Cut-off probability which defines a "near-optimal" sample size.
##'   The minimum and maximum sample size for which the ENBS is within
##'   \code{pcut} (by default 5%) of its maximum value will be determined.
##'
##' @param smooth If \code{TRUE}, then the maximum ENBS is determined after
##'   fitting a nonparametric regression to the data frame \code{x}, which
##'   estimates and smooths the ENBS for every integer sample size in the range
##'   of \code{x$n}.  The regression is done using the default settings of
##'   \code{\link[mgcv]{gam}} from the \pkg{mgcv} package. 
##'
##'   If this is \code{FALSE}, then no smoothing or interpolation is done, and
##'   the maximum is determined by searching over the values supplied in
##'   \code{x}.
##'
##' @param smooth_df Basis dimension for the smooth regression. Passed as the
##' `k` argument to the `s()` term in \code{\link[mgcv]{gam}}.   Defaults to
##' 6, or the number of unique sample sizes minus 1 if this is lower.  Set
##' to a higher number if you think the smoother does not capture the
##' relation of ENBS to sample size accurately enough. 
##'
##' @inheritParams pop_voi
##'
##' @return Data frame with components \code{enbs} giving the ENBS, and
##'   \code{sd} giving the corresponding standard deviation.  The rows of the
##'   data frame correspond to the rows of \code{evsi}, and any \code{n} and
##'   \code{k} are inherited from \code{evsi}.  Additional columns include:
##'
##'   \code{pce}: the probability that the study is cost-effective, i.e. that
##'   the ENBS is positive, obtained from a normal distribution defined by the
##'   estimate and standard deviation.
##'   
##'   \code{enbsmax}: The maximum ENBS for each willingness-to-pay \code{k}.
##'   
##'   \code{nmax}: The sample size \code{n} at which this maximum is achieved.
##'
##' A second data frame is returned as the \code{"enbsmax"} attribute.
##' This has one row per willingness-to-pay (`k`), giving the optimal
##' ENBS (`enbsmax`) the optimal sample size (`nmax`) and an interval
##' estimate for the optimal sample size (`nlower` to `nupper`).
##'
##' If \code{pop}, \code{time} or \code{dis} were supplied as vectors
##' of more than one element, then additional columns will be returned
##' in these data frames to identify the population, time or discount
##' rate for each ENBS calculation.  An index \code{ind} is also returned
##' to identify the unique combination that each row refers to.
##'   
##' @references Value of Information for Healthcare Decision Making
##' (CRC Press, eds. Heath, Kunst and Jackson: forthcoming)
##'
##' @export
enbs <- function(evsi, costs_setup, costs_pp, pop, time, dis=0.035,
                 smooth=FALSE, smooth_df=NULL, pcut=0.05){
  costs_setup <- costs_elic(costs_setup)
  costs_pp <- costs_elic(costs_pp)

  ## Handle vectorised population size, time horizon and discount
  pdt <- expand.grid(pop=pop, time=time, dis=dis)
  npdt <- nrow(pdt)
  pdt$ind <- 1:nrow(pdt)
  pdt <- pdt[rep(pdt$ind, each=nrow(evsi)),]
  evsi <- evsi[rep(1:nrow(evsi), npdt),]
  for (i in c("pop","time","dis"))
    evsi[[i]] <- pdt[[i]]
    
  ## Calculate population EVSI, costs, ENBS, and their SEs
  pop_evsi <- pop_voi(evsi$evsi, evsi$pop, evsi$time, dis)
  costs <- costs_setup["mean"] + evsi$n * costs_pp["mean"]
  enbs <- pop_evsi - costs
  evsi_sd <- if (is.null(evsi$sd)) 0 else evsi$sd
  costs_sd <- sqrt(costs_setup["sd"]^2 + evsi$n^2  *  costs_pp["sd"]^2)
  pop_mult <- (pop_evsi / evsi$evsi)^2
  enbs_sd <- sqrt((pop_mult*evsi_sd)^2 + costs_sd^2)
  enbs <- data.frame(n=evsi$n, k=evsi$k, enbs=enbs, sd=enbs_sd,
                     pce = pnorm(0, enbs, enbs_sd, lower.tail=FALSE))

  ## Unique combinations of WTP, population, time and discount to optimise for
  evsi$ind <- interaction(evsi$k, evsi$pop, evsi$time, evsi$dis)
  evsi$ind <- enbs$ind <- match(evsi$ind, unique(evsi$ind))
  ind_lookup <- evsi[,c("ind","k","pop","time","dis")]
  ind_lookup <- ind_lookup[!duplicated(ind_lookup$ind),]

  ## Determine the optimal sample size for each of these
  maxlist <- lapply(split(enbs, enbs$ind), enbs_opt, pcut=pcut,
                    smooth=smooth, smooth_df=smooth_df)
  enbsmax <- do.call(rbind, maxlist)
  rownames(enbsmax) <- NULL
  has_combs <- 0
  for (i in c("k", "pop", "time", "dis")){
    if (length(unique(ind_lookup[[i]])) > 1){
      enbsmax[[i]] <- ind_lookup[[i]][match(enbsmax$ind, ind_lookup$ind)]
      enbs[[i]] <- ind_lookup[[i]][match(enbs$ind, ind_lookup$ind)]
      has_combs <- has_combs + 1
    }
  }
  for (i in c("enbsmax","nmax","nlower","nupper"))
    enbs[[i]] <- enbsmax[[i]][match(enbs$ind, enbsmax$ind)]
  if (has_combs > 1) enbsmax$ind <- enbs$ind <- NULL
  attr(enbs,"enbsmax") <- enbsmax
  enbs
}

##' Determine the optimum sample size in an analysis of the expected net benefit
##' of sampling
##'
##' The optimum sample size for a given willingness to pay is determined either
##' by a simple search over the supplied ENBS estimates for different sample
##' sizes, or by a regression and interpolation method.
##'
##' @param x Data frame containing a set of ENBS estimates for
##'   different sample sizes, which will be optimised over.  Usually
##'   this is for a common willingness-to-pay. The required components
##'   are \code{enbs} and \code{n}.
##'
##' @param keep_preds If \code{TRUE} and \code{smooth=TRUE} then the data frame of
##' predictions from the smooth regression model is stored in the \code{"preds"}
##' attribute of the result.
##'
##' @inheritParams enbs
##' 
##' @return A data frame with one row, and the following columns:
##' 
##' \code{ind}: An integer index identifying, e.g. the willingness to pay and other common characteristics of the ENBS estimates (e.g. incident population size, decision time horizon). This is copied from \code{x$ind}.
##' 
##' \code{enbsmax}: the maximum ENBS
##' 
##' \code{nmax}: the sample size at which this maximum is achieved
##' 
##' \code{nlower}: the lowest sample size for which the ENBS is within
##' 
##' \code{pcut} (default 5%) of its maximum value
##' 
##' \code{nupper}: the corresponding highest ENBS
##'
##' @export
enbs_opt <- function(x, pcut=0.05, smooth=FALSE, smooth_df=NULL, keep_preds=FALSE){
  if (smooth) { 
    nrange <- seq(min(x$n), max(x$n), by=1)
    if (is.null(smooth_df)) smooth_df <- min(6, length(unique(x$n)) - 1)
    mod <- mgcv::gam(enbs~s(n, k=smooth_df), data=x)
    enbs_smooth <- predict(mod, newdata=list(n=nrange))
    x <- data.frame(n=nrange, enbs=enbs_smooth,
                    ind = if (is.null(x$ind)) 1 else x$ind[1]
                    )
  } 
  maxind <- which.max(x$enbs)
  x$enbsmax <- x$enbs[maxind]
  x$nmax <- x$n[maxind]
  near_max <- x$n[x$enbs > x$enbsmax - abs(pcut*x$enbsmax)]
  x$nlower <- min(near_max)
  x$nupper <- max(near_max)
  res <- x[maxind,,drop=FALSE]
  res <- res[,c("ind","enbsmax","nmax","nlower","nupper"),drop=FALSE]
  if (keep_preds) attr(res, "preds") <- x
  res
}

##' Population expected value of information
##'
##' Convert per-person expected value of information to the population
##' expected value of information, given a discount rate over some
##' time horizon.
##'
##' Calculated as \code{voi*pop/dis*(1 - exp(-dis*time))}, or \code{voi*pop}
##' if the discount rate is zero.  This is a continuous-time variant
##' of the typical discrete-time discounting formula.
##'
##' Any arguments may be supplied as vectors, in which case, all
##' arguments are replicated to the length of the longest argument.
##'
##' @param voi Vector of estimates of any per-person value of information
##' measure, e.g. the \code{evsi} column of the data frame returned by
##' \code{\link{evsi}} or the correspondingly-named columns of the
##' data frames returned by \code{\link{evppi}} or \code{\link{evpi}}. 
##' 
##' @param pop Size of the population who would be affected by the decision.
##'
##' @param time Time horizon over which discounting will be applied.
##'
##' @param dis Discount rate used when converting per-person to population EVSI.
##'
##' @return A vector of population VoI estimates.
##' 
##' @export
pop_voi <- function(voi, pop, time, dis=0.035){

  ## vectorising...
  nmax <- max(length(voi), length(pop), length(time), length(dis))
  voi <- rep(voi, length.out=nmax)
  pop <- rep(pop, length.out=nmax)
  time <- rep(time, length.out=nmax)
  dis <- rep(dis, length.out=nmax)

  ifelse(dis==0,
         voi*pop, 
         voi*pop/dis*(1 - exp(-dis*time)))
}

costs_elic <- function(costs){
  if (!is.numeric(costs)) stop("costs should be numeric")
  if (!length(costs) %in% c(1,2,3))
    stop("length of costs argument should be 1, 2 or 3")
  costs <- sort(costs)
  if (length(costs)==1) {
    cmean <- costs
    csd <- clsd <- 0
  }
  if (length(costs)==2) {
    cmean <- mean(costs)
    csd <- (costs[2] - costs[1]) / 4
    clsd <- (log(costs[2]) - log(costs[1])) / 4
  }
  if (length(costs)==3) {
    cmean <- costs[2]
    csd <- (costs[3] - costs[1]) / 4
    clsd <- (log(costs[3]) - log(costs[1])) / 4
  }
  ## log versions currently unused.
  c(mean=cmean, sd=csd, lmean=log(cmean), lsd=clsd) 
}
