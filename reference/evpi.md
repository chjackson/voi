# Calculate the expected value of perfect information from a decision model

Calculate the expected value of perfect information from a decision
model using standard Monte Carlo simulation

## Usage

``` r
evpi(outputs, nsim = NULL)
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

- nsim:

  Number of simulations from the decision model to use for calculating
  EVPPI. The first `nsim` rows of the objects in `inputs` and `outputs`
  are used.

## Value

The expected value of perfect information, either as a single value, or
a data frame indicating the value for each willingness-to-pay.
