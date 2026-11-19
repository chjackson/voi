# Calculate the expected value of partial perfect information from a decision-analytic model

Calculate the expected value of partial perfect information from a
decision-analytic model

## Usage

``` r
evppi(
  outputs,
  inputs,
  pars = NULL,
  method = NULL,
  se = FALSE,
  B = 1000,
  nsim = NULL,
  verbose = FALSE,
  check = FALSE,
  ...
)
```

## Arguments

- outputs:

  This could take one of two forms

  "net benefit" form: a matrix or data frame of samples from the
  uncertainty distribution of the expected net benefit. The number of
  rows should equal the number of samples, and the number of columns
  should equal the number of decision options.

  "cost-effectiveness analysis" form: a list with the following named
  components:

  `"c"`: a matrix or data frame of samples from the distribution of
  costs. There should be one column for each decision option.

  `"e"`: a matrix or data frame of samples from the distribution of
  effects, likewise.

  `"k"`: a vector of willingness-to-pay values.

  Objects of class `"bcea"`, as created by the BCEA package, are in this
  "cost-effectiveness analysis" format, therefore they may be supplied
  as the `outputs` argument.

  Users of heemod can create an object of this form, given an object
  produced by `run_psa` (`obj`, say), with
  [`import_heemod_outputs`](https://chjackson.github.io/voi/reference/import_heemod.md).

  If `outputs` is a matrix or data frame, it is assumed to be of "net
  benefit" form. Otherwise if it is a list, it is assumed to be of "cost
  effectiveness analysis" form.

- inputs:

  Matrix or data frame of samples from the uncertainty distribution of
  the input parameters of the decision model. The number of columns
  should equal the number of parameters, and the columns should be
  named. This should have the same number of rows as there are samples
  in `outputs`, and each row of the samples in `outputs` should give the
  model output evaluated at the corresponding parameters.

  Users of heemod can create an object of this form, given an object
  produced by `run_psa` (`obj`, say), with
  [`import_heemod_inputs`](https://chjackson.github.io/voi/reference/import_heemod.md).

- pars:

  Either a character vector, or a list of character vectors.

  If a character vector is supplied, then a single, joint EVPPI
  calculation is done with for the parameters named in this vector.

  If a list of character vectors is supplied, then multiple EVPPI
  calculations are performed, one for each list component defined in the
  above vector form.

  `pars` must be specified if `inputs` is a matrix or data frame. This
  should then correspond to particular columns of `inputs`. If `inputs`
  is a vector, this is assumed to define the single parameter of
  interest, and then `pars` is not required.

- method:

  Character string indicating the calculation method. If one string is
  supplied, this is used for all calculations. A vector of different
  strings can be supplied if a different method is desired for different
  list components of `pars`.

  The default methods are based on nonparametric regression:

  `"gam"` for a generalized additive model implemented in the
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) function from the mgcv
  package. This is the default method for calculating the EVPPI of 4 or
  fewer parameters.

  `"gp"` for a Gaussian process regression, as described by Strong et
  al. (2014) and implemented in the SAVI package
  (<https://github.com/Sheffield-Accelerated-VoI/SAVI>). This is the
  default method for calculating the EVPPI of more than 4 parameters.

  `"inla"` for an INLA/SPDE Gaussian process regression method, from
  Heath et al. (2016).

  `"bart"` for Bayesian additive regression trees, using the dbarts
  package. Particularly suited for joint EVPPI of many parameters.

  `"earth"` for a multivariate adaptive regression spline with the earth
  package (Milborrow, 2019).

  `"so"` for the method of Strong and Oakley (2013). Only supported for
  single parameter EVPPI.

  `"sal"` for the method of Sadatsafavi et al. (2013). Only supported
  for single parameter EVPPI.

- se:

  If this is `TRUE`, calculate a standard error for the EVPPI if
  possible. Currently only supported for methods `"gam"`, `"earth"` and
  `method="bart"`. (In the latter method it is more correctly called a
  posterior standard deviation). These represent uncertainty about the
  parameters of the fitted regression model, and will naturally be lower
  when more simulations from the decision model are used to fit it. They
  do not represent uncertainty about the structure of the regression
  model, and are also typically small in the context of uncertainties
  arising from converting individual-level to population VoI.

- B:

  Number of parameter replicates for calculating the standard error.
  Only applicable to `method="gam"`. For `method="bart"` the analogous
  quantity is the number of MCMC samples, which is controlled by the
  `ndpost` argument to
  [`bart`](https://rdrr.io/pkg/dbarts/man/bart.html), which can be
  passed as an argument to `evppi`.

- nsim:

  Number of simulations from the decision model to use for calculating
  EVPPI. The first `nsim` rows of the objects in `inputs` and `outputs`
  are used.

- verbose:

  If `TRUE`, then messages are printed describing each step of the
  calculation, if the method supplies these. Can be useful to see the
  progress of slow calculations.

- check:

  If `TRUE`, then extra information about the estimation is saved inside
  the object that this function returns. This currently only applies to
  the regression-based methods `"gam"` and `"earth"` where the fitted
  regression model objects are saved. This allows use of the
  [`check_regression`](https://chjackson.github.io/voi/reference/check_regression.md)
  function, which produces some diagnostic checks of the regression
  models.

- ...:

  Other arguments to control specific methods.

  For `method="gam"`, the following arguments can be supplied:

  - `gam_formula`: a character string giving the right hand side of the
    formula supplied to the `gam()` function. By default, this is a
    tensor product of all the parameters of interest, e.g. if
    `pars = c("pi","rho")`, then `gam_formula` defaults to
    `t(pi, rho, bs="cr")`. The option `bs="cr"` indicates a cubic spline
    regression basis, which is more computationally efficient than the
    default "thin plate" basis. If there are four or more parameters of
    interest, then the additional argument `k=4` is supplied to `te()`,
    specifying a four-dimensional basis, which is currently the default
    in the SAVI package.

    If there are spaces in the variable names in `inputs`, then these
    should be converted to underscores before forming an explicit
    `gam_formula`.

  For `method="earth"`:

  - `degree`. Maximum degree of interactions to include in the
    regression. By default this is 1. A better-fitting regression
    approximation to EVPPI might be obtained by using `degree=2` or
    higher, at the expense of fitting speed.

  For `method="gp"`, the following arguments can be supplied:

  - `gp_hyper_n`: number of samples to use to estimate the
    hyperparameters in the Gaussian process regression method. By
    default, this is the minimum of the following three quantities: 30
    times the number of parameters of interest, 250, and the number of
    simulations being used for calculating EVPPI.

  - `maxSample`: Maximum sample size to employ for `method="gp"`. Only
    increase this from the default 5000 if your computer has sufficent
    memory to invert square matrices with this dimension.

  For `method="inla"`, the following arguments can be supplied, as
  described in detail in Baio, Berardi and Heath (2017):

  - `int.ord` (integer) maximum order of interaction terms to include in
    the regression predictor, e.g. if `int.ord=k` then all k-way
    interactions are used. Currently this applies to both effects and
    costs.

  - `cutoff` (default 0.3) controls the density of the points inside the
    mesh in the spatial part of the mode. Acceptable values are
    typically in the interval (0.1,0.5), with lower values implying more
    points (and thus better approximation and greatercomputational
    time).

  - `convex.inner` (default = -0.4) and `convex.outer` (default = -0.7)
    control the boundaries for the mesh. These should be negative values
    and can be decreased (say to -0.7 and -1, respectively) to increase
    the distance between the points and the outer boundary, which also
    increases precision and computational time.

  - `robust`. if `TRUE` then INLA will use a t prior distribution for
    the coefficients of the linear predictor, rather than the default
    normal distribution.

  - `h.value` (default=0.00005) controls the accuracy of the INLA
    grid-search for the estimation of the hyperparameters. Lower values
    imply a more refined search (and hence better accuracy), at the
    expense of computational speed.

  - `plot_inla_mesh` (default `FALSE`) Produce a plot of the mesh.

  - `max.edge` Largest allowed triangle edge length when constructing
    the mesh, passed to
    [`inla.mesh.2d`](https://rdrr.io/pkg/INLA/man/inla.mesh.2d.html).

  - `pfc_struc` Variance structure to pass to `pfc` in the ldr package
    for principal fitted components. The default `"AIC"` selects the one
    that fits best given two basis terms. Change this to, e.g. `"iso"`,
    `"aniso"` or `"unstr"` if an "Error in eigen..." is obtained.

  For any of the nonparametric regression methods:

  - `ref` The reference decision option used to define the incremental
    net benefit, cost or effects before performing nonparametric
    regression. Either an integer column number, or the name of the
    column from `outputs`.

  For `method="so"`:

  - `n.blocks` Number of blocks to split the sample into. Required.

  For `method="sal"`:

  - `n.seps` Number of separators (default 1).

## Value

A data frame with a column `pars`, indicating the parameter(s), and a
column `evppi`, giving the corresponding EVPPI.

If `outputs` is of "cost-effectiveness analysis" form, so that there is
one EVPPI per willingness-to-pay value, then a column `k` identifies the
willingness-to-pay.

If standard errors are requested, then the standard errors are returned
in the column `se`.

## References

Heath, A., Kunst, N., & Jackson, C. (eds.). (2024). Value of Information
for Healthcare Decision-Making. CRC Press.

Strong, M., Oakley, J. E., & Brennan, A. (2014). Estimating
multiparameter partial expected value of perfect information from a
probabilistic sensitivity analysis sample: a nonparametric regression
approach. Medical Decision Making, 34(3), 311-326.

Heath, A., Manolopoulou, I., & Baio, G. (2016). Estimating the expected
value of partial perfect information in health economic evaluations
using integrated nested Laplace approximation. Statistics in Medicine,
35(23), 4264-4280.

Baio, G., Berardi, A., & Heath, A. (2017). Bayesian cost-effectiveness
analysis with the R package BCEA. New York: Springer.

Milborrow, S. (2019) earth: Multivariate Adaptive Regression Splines. R
package version 5.1.2. Derived from mda:mars by Trevor Hastie and Rob
Tibshirani. Uses Alan Miller's Fortran utilities with Thomas Lumley's
leaps wrapper. https://CRAN.R-project.org/package=earth.

Strong, M., & Oakley, J. E. (2013). An efficient method for computing
single-parameter partial expected value of perfect information. Medical
Decision Making, 33(6), 755-766. Chicago

Sadatsafavi, M., Bansback, N., Zafari, Z., Najafzadeh, M., & Marra, C.
(2013). Need for speed: an efficient algorithm for calculation of
single-parameter expected value of partial perfect information. Value in
Health, 16(2), 438-448.
