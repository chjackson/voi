# Calculate the expected value of partial perfect information for an estimation problem

Calculate the expected value of partial perfect information for an
estimation problem. This computes the expected reduction in variance in
some quantity of interest with perfect information about a parameter or
parameters of interest.

## Usage

``` r
evppivar(
  outputs,
  inputs,
  pars = NULL,
  method = NULL,
  nsim = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- outputs:

  a vector of values for the quantity of interest, sampled from the
  uncertainty distribution of this quantity that is induced by the
  uncertainty about the parameters. This can also be a data frame with
  one column.

  Typically this will come from a Monte Carlo sample, where we first
  sample from the uncertainty distributions of the parameters, and then
  compute the quantity of interest as a function of the parameters. It
  might also be produced by a Markov Chain Monte Carlo sample from the
  joint distribution of parameters and outputs.

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

- nsim:

  Number of simulations from the decision model to use for calculating
  EVPPI. The first `nsim` rows of the objects in `inputs` and `outputs`
  are used.

- verbose:

  If `TRUE`, then messages are printed describing each step of the
  calculation, if the method supplies these. Can be useful to see the
  progress of slow calculations.

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

## References

Jackson, C., Presanis, A., Conti, S., & De Angelis, D. (2019). Value of
information: Sensitivity analysis and research design in Bayesian
evidence synthesis. Journal of the American Statistical Association,
114(528), 1436-1449.

Jackson, C., Johnson, R., de Nazelle, A., Goel, R., de Sa, T. H.,
Tainio, M., & Woodcock, J. (2021). A guide to value of information
methods for prioritising research in health impact modelling.
Epidemiologic Methods, 10(1).

Jackson, C. H., Baio, G., Heath, A., Strong, M., Welton, N. J., &
Wilson, E. C. (2022). Value of Information analysis in models to inform
health policy. Annual Review of Statistics and its Application, 9,
95-118.
