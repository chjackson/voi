##' Import results of probabilistic analysis from heemod
##'
##' [heemod](https://CRAN.R-project.org/package=heemod) is a package
##' for constructing common forms of health economic decision models.
##' The outputs from probabilistic analysis of these models can be
##' imported using these functions, to allow Value of Information
##' measures to be calculated for them using the \pkg{voi} package.
##'
##' @aliases import_heemod_outputs import_heemod_inputs
##'
##' @param obj Object returned by the \code{\link[heemod]{run_psa}}
##'   function in \pkg{heemod}, containing samples from probabilistic
##'   analysis of a decision model.
##'
##' @param k Vector of willingness-to-pay values.  The default is
##' inherited from the \code{\link[BCEA]{bcea}} function from the \pkg{BCEA}
##' package.
##'
##' @return \code{import_heemod_outputs} produces a list of model
##'   outputs in "cost-effectiveness analysis" format, that can be
##'   supplied as the \code{outputs} argument to \code{\link{evppi}}
##'   and similar functions in the \pkg{voi} package.  Both the
##'   \pkg{heemod} and \pkg{BCEA} packages need to be installed to use
##'   this.
##'
##' \code{import_heemod_inputs} produces a data frame with samples of
##' parameter values under uncertainty, that can be supplied as the
##' \code{inputs} argument to \code{\link{evppi}} and similar functions
##' in \pkg{voi}.
##'
##' @name import_heemod
NULL

##' @rdname import_heemod
##' @export
import_heemod_outputs <- function(obj, k=NULL){
  if (!requireNamespace("heemod",quietly=TRUE)) {
    stop("The `heemod` package is required")
  }
  if (!requireNamespace("BCEA",quietly=TRUE)) {
    stop("The 'BCEA' package is required")
  }
  heemod::run_bcea(obj,k=k,ref=1)[c("c","e","k")]
}

##' @rdname import_heemod
##' @export
import_heemod_inputs <- function(obj){
  obj$psa[!duplicated(obj$psa$.index), obj$resamp_par]
}
