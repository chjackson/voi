# Traditional two-level Monte Carlo estimator of EVPPI.

Traditional two-level Monte Carlo estimator of the expected value of
partial perfect information from a decision-analytic model. Only useful
in the simplest of examples. For realistically complex examples, the
methods implemented in the
[`evppi`](https://chjackson.github.io/voi/reference/evppi.md) function,
based on regression, will usually be much more computationally
efficient.

## Usage

``` r
evppi_mc(
  model_fn,
  par_fn,
  pars,
  nouter,
  ninner,
  k = NULL,
  mfargs = NULL,
  verbose = FALSE
)
```

## Arguments

- model_fn:

  A function to evaluate a decision-analytic model at a given set of
  parameters. This should have one argument per parameter, and return
  either:

  (net benefit format) a vector giving the net benefit for each decision
  option, or

  (cost-effectiveness analysis format) a matrix or data frame with two
  rows, and one column for each decision option. If the rows have names
  `"e"` and `"c"` then these are assumed to be the effects and costs
  respectively.

  Otherwise, the first row is assumed to be the effects, and the second
  the costs.

- par_fn:

  A function to generate a random sample of values for the parameters of
  `model_fn`. This should return a matrix or a data frame with named
  columns matching the arguments of `model_fn`.

  If any required arguments to `model_fn` are not supplied in this
  return value, then `evppi_mc` looks for them in the list supplied as
  the `mfargs` argument.

  If any required arguments are not found in the results of `par_fn` or
  `mfargs`, and if `model_fn` defines default values for those
  arguments, then those default values are used.

  The first argument of `par_fn` should be an integer `n` denoting the
  number of random values to draw for each parameter. The object
  returned by `par_fn` should then have `n` rows, and one column for
  each parameter. If one value is drawn, then `par_fn` is also allowed
  to return a vector, but this should still be named.

  The parameters may be correlated. If we wish to compute the EVPPI for
  a parameter which is correlated with a different parameter q, then
  `par_fn` must have an argument with the name of that parameter. If
  that argument is set to a fixed value, then `par_fn` should return a
  sample drawn conditionally on that value. If that argument is not
  supplied, then `par_fn` must return a sample drawn from the marginal
  distribution. See the vignette for an example.

- pars:

  A character vector giving the parameters of interest, for which the
  EVPPI is required. This should correspond to an explicit argument to
  `model_fn`.

  The parameters of interest are assumed to have uncertainty
  distributions that are independent of those of the other parameters.

- nouter:

  Number of outer samples

- ninner:

  Number of inner samples

- k:

  Vector of willingness-to-pay values. Only used if `model_fn` is in
  cost-effectiveness analyis format.

- mfargs:

  Named list of additional arguments to supply to `model_fn`.

- verbose:

  Set to `TRUE` to print some additional messages to help with
  debugging.

## Value

A data frame with a column `pars`, indicating the parameter(s), and a
column `evppi`, giving the corresponding EVPPI.

If `outputs` is of "cost-effectiveness analysis" form, so that there is
one EVPPI per willingness-to-pay value, then a column `k` identifies the
willingness-to-pay.

## Details

See the [package overview / Get Started
vignette](https://chjackson.github.io/voi/articles/voi.html#evppimc) for
an example of using this function.
