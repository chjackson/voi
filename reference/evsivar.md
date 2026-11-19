# Calculate the expected value of sample information for an estimation problem

Calculate the expected value of sample information for an estimation
problem. This computes the expected reduction in variance in some
quantity of interest from a study of a certain design that informs the
parameters of interest.

## Usage

``` r
evsivar(
  outputs,
  inputs,
  study = NULL,
  datagen_fn = NULL,
  pars = NULL,
  n = 100,
  aux_pars = NULL,
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
  uncertainty about the parameters.

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

- study:

  Name of one of the built-in study types supported by this package for
  EVSI calculation. If this is supplied, then the columns of `inputs`
  that correspond to the parameters governing the study data should be
  identified in `pars`.

  Current built-in studies are

  `"binary"` A study with a binary outcome observed on one sample of
  individuals. Requires one parameter: the probability of the outcome.
  The sample size is specifed in the `n` argument to
  [`evsi()`](https://chjackson.github.io/voi/reference/evsi.md), and the
  binomially-distributed outcome is named `X1`.

  `"trial_binary"` Two-arm trial with a binary outcome. Requires two
  parameters: the probability of the outcome in arm 1 and 2
  respectively. The sample size is the same in each arm, specifed in the
  `n` argument to
  [`evsi()`](https://chjackson.github.io/voi/reference/evsi.md), and the
  binomial outcomes are named `X1` and `X2` respectively.

  `"normal_known"` A study of a normally-distributed outcome, with a
  known standard deviation, on one sample of individuals. Likewise the
  sample size is specified in the `n` argument to
  [`evsi()`](https://chjackson.github.io/voi/reference/evsi.md). The
  standard deviation defaults to 1, and can be changed by specifying
  `sd` as a component of the `aux_pars` argument, e.g.
  `evsi(..., aux_pars=list(sd=2))`.

  Either `study` or `datagen_fn` should be supplied to
  [`evsi()`](https://chjackson.github.io/voi/reference/evsi.md).

  For the EVSI calculation methods where explicit Bayesian analyses of
  the simulated data are performed, the prior parameters for these
  built-in studies are supplied in the `analysis_args` argument to
  [`evsi()`](https://chjackson.github.io/voi/reference/evsi.md). These
  assume Beta priors for probabilities, and Normal priors for the mean
  of a normal outcome.

- datagen_fn:

  If the proposed study is not one of the built-in types supported, it
  can be specified in this argument as an R function to sample predicted
  data from the study. This function should have the following
  specification:

  1.  the function's first argument should be a data frame of parameter
      simulations, with one row per simulation and one column per
      parameter. The parameters in this data frame must all be found in
      `inputs`, but need not necessarily be in the same order or include
      all of them.

  2.  the function should return a data frame.

  3.  the returned data frame should have number of rows equal to the
      number of parameter simulations in `inputs`.

  4.  if `inputs` is considered as a sample from the posterior, then
      `datagen_fn(inputs)` returns a corresponding sample from the
      posterior predictive distribution, which includes two sources of
      uncertainty: (a) uncertainty about the parameters and (b) sampling
      variation in observed data given fixed parameter values.

  5.  the function can optionally have more than one argument. If so,
      these additional arguments should be given default values in the
      definition of `datagen_fn`. If there is an argument called `n`,
      then it is interpreted as the sample size for the proposed study.

- pars:

  Character vector identifying which parameters are learned from the
  proposed study. This is required for the moment matching and
  importance sampling methods, and these should be columns of `inputs`.
  This is not required for the nonparametric regression methods.

- n:

  Sample size of future study, or vector of alternative sample sizes.
  This is understood by the built-in study designs. For studies
  specified by the user with `datagen_fn`, if `datagen_fn` has an
  argument `n`, then this is interpreted as the sample size. However if
  calling `evsi` for a user-specified design where `datagen_fn` does not
  have an `n` argument, then any `n` argument supplied to `evsi` will be
  ignored.

  Currently this shortcut is not supported if more than one quantity is
  required to describe the sample size, for example, trials with
  unbalanced arms. In that case, you will have to hard-code the required
  sample sizes into `datagen_fn`.

  For the nonparametric regression and importance sampling methods, the
  computation is simply repeated for each sample size supplied here.

  The moment matching method uses a regression model to estimate the
  dependency of the EVSI on the sample size, hence to enable EVSI to be
  calculated efficiently for any number of sample sizes (Heath et al.
  2019).

- aux_pars:

  A list of additional fixed arguments to supply to the function to
  generate the data, whether that is a built-in study design or
  user-defined function supplied in `datagen_fn`. For example,
  `evsi(..., aux_pars = list(sd=2))` defines the fixed standard
  deviation in the `"normal_known"` model.

- method:

  See [`evsi`](https://chjackson.github.io/voi/reference/evsi.md), only
  nonparametric regression methods are currently supported in `evsivar`.

- nsim:

  Number of simulations from the model to use for calculating EVPPI. The
  first `nsim` rows of the objects in `inputs` and `outputs` are used.

- verbose:

  If `TRUE`, then messages are printed describing each step of the
  calculation, if the method supplies these. Can be useful to see the
  progress of slow calculations.

- ...:

  Other arguments understood by specific methods, e.g. `gam_formula` and
  other controlling options (see
  [`evppi`](https://chjackson.github.io/voi/reference/evppi.md)) can be
  passed to the nonparametric regression used inside the moment matching
  method.

## Value

A data frame with a column `pars`, indicating the parameter(s), and a
column `evsi`, giving the corresponding EVSI. If there are EVSI
estimates for multiple sample sizes, the sample size is returned in the
column `n`.

## References

Jackson, C., Presanis, A., Conti, S., & De Angelis, D. (2019). Value of
information: Sensitivity analysis and research design in Bayesian
evidence synthesis. Journal of the American Statistical Association,
114(528), 1436-1449.
