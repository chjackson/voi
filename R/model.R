## Facilities to deal with functions that evaluate decision-analytic models, and functions that
## generate parameters for those models. 

##' Check that a decision-analytic model function is of the appropriate form.  Detect if it returns net benefit or CEA format
##' note check_outputs adds a class 
##'
##' @return A modified copy of `model_fn` with a \code{class} attribute indicating whether it is
##' in net benefit \code{"nb"} or cost-effectiveness \code{"cea"} format, and an attribute \code{"nopt"} giving the
##' number of decision options. 
##' 
##' @keywords internal 
check_model_fn <- function(model_fn, par_fn, mfargs){
    ## Test by evaluating the model function at a single set of parameters/arguments
    pars <- check_parfn1(par_fn, model_fn, mfargs)
    defaults <- get_default_args(model_fn, pars)
    res <- do.call(model_fn, c(pars, mfargs, defaults)[names(formals(model_fn))])

    if (is.vector(res)){
        class(model_fn) <- c("nb", attr(model_fn, "class"))
        attr(model_fn, "nopt") <- length(res)
    }  else if (is.matrix(res) || is.data.frame(res)) {
        class(model_fn) <- c("cea", attr(model_fn, "class"))
        attr(model_fn, "nopt") <- ncol(res)
        if (nrow(res) != 2) {
            stop("If `model_fn` returns a matrix or data frame it should have two rows, one for effects and one for costs")
        }
    } else stop("`model_fn` should return a vector, matrix or data frame")
    describe_modelfn(model_fn)
    if (attr(model_fn, "nopt")==1)
        stop("model_fn should describe more than one decision option")
    model_fn
}

describe_modelfn <- function(model_fn, ...){
    UseMethod("describe_modelfn", model_fn)
}
describe_modelfn.nb <- function(model_fn, ...){
    plural <- if (attr(model_fn, "nopt") > 1) "s" else ""
    message(sprintf("model_fn returns net benefit for %s decision option%s", attr(model_fn, "nopt"), plural))
}
describe_modelfn.cea <- function(model_fn, ...){
    plural <- if (attr(model_fn, "nopt") > 1) "s" else ""
    message(sprintf("model_fn returns effects and costs for %s decision option%s", attr(model_fn, "nopt"), plural))
}

check_parfn1 <- function(par_fn, model_fn, mfargs){
    fn_try <- try(pars <- par_fn(1), silent=TRUE)
    if (inherits(fn_try, "try-error")){
        stop("Evaluating `par_fn` returned an error")
    }
    if (is.vector(pars)) {
        if (is.null(names(pars)))
            stop("pars(1) should return a named vector, matrix or data frame")
        pars <- as.data.frame(as.list(pars))
    } else if (is.matrix(pars)){
            if (is.null(colnames(pars)))
                stop("If pars(1) returns a matrix, the columns should be named")
            pars <- as.data.frame(pars)
    } else if (!is.data.frame(pars)) {
        stop("pars(1) should return a named vector, matrix or data frame")
    }
    model_pars <- names(formals(model_fn))
    defaults <- get_default_args(model_fn, pars)
    missing_pars <- setdiff(model_pars, c(names(pars),names(mfargs),names(defaults)))
    if (length(missing_pars) > 0)
        stop("The following parameters of `model_fn` were not found in the components of pars(1) or in the `...` argument: ",paste(missing_pars,collapse=","))
    pars 
}

# get the arguments to function fn that are not supplied in `supplied`

get_default_args <- function(fn, supplied=NULL){
    fs <- formals(fn)
    has_default <- sapply(fs, function(x) { if (is.name(x) && !nzchar(x)) FALSE else TRUE } )
    default_args <- as.list(fs[has_default])
    if (!is.null(supplied))
        default_args <- default_args[!(names(default_args) %in% names(supplied))]
    default_args
}

check_parfnn <- function(par_fn, model_fn){
    fn_try <- try(pars <- par_fn(2), silent=TRUE)
    if (inherits(fn_try, "try-error")){
        stop("Evaluating `par_fn` returned an error") # TODO print the actual error? 
    }
    if (!(is.matrix(pars) || is.data.frame(pars)))
        stop("par_fn(n) for n>1 should return a matrix or data frame")
    if (nrow(pars) != 2) {
        stop("par_fn(n) should have n rows")
    }
    as.data.frame(pars) # return value currently unused
}

## todo handle more indices
mfi <- function(res){
    ci <- match("c", rownames(res))
    if (is.na(ci)) ci <- 2
    ei <- match("e", rownames(res))
    if (is.na(ei)) ei <- 1
    list(c=ci, e=ei)
}
