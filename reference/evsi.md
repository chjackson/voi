# Calculate the expected value of sample information from a decision-analytic model

Calculate the expected value of sample information from a
decision-analytic model

## Usage

``` r
evsi(
  outputs,
  inputs,
  study = NULL,
  datagen_fn = NULL,
  pars = NULL,
  pars_datagen = NULL,
  n = 100,
  aux_pars = NULL,
  method = NULL,
  likelihood = NULL,
  analysis_fn = NULL,
  analysis_args = NULL,
  model_fn = NULL,
  par_fn = NULL,
  Q = 50,
  npreg_method = "gam",
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

- study:

  Name of one of the built-in study types supported by this package for
  EVSI calculation. If this is supplied, then the columns of `inputs`
  that correspond to the parameters governing the study data should be
  identified in `pars`.

  Current built-in studies are

  `"binary"` A study with a binary outcome observed on one sample of
  individuals. Requires one parameter: the probability of the outcome.
  The sample size is specifed in the `n` argument to `evsi()`, and the
  binomially-distributed outcome is named `X1`.

  `"trial_binary"` Two-arm trial with a binary outcome. Requires two
  parameters: the probability of the outcome in arm 1 and 2
  respectively. The sample size is the same in each arm, specifed in the
  `n` argument to `evsi()`, and the binomial outcomes are named `X1` and
  `X2` respectively.

  `"normal_known"` A study of a normally-distributed outcome, with a
  known standard deviation, on one sample of individuals. Likewise the
  sample size is specified in the `n` argument to `evsi()`. The standard
  deviation defaults to 1, and can be changed by specifying `sd` as a
  component of the `aux_pars` argument, e.g.
  `evsi(..., aux_pars=list(sd=2))`.

  Either `study` or `datagen_fn` should be supplied to `evsi()`.

  For the EVSI calculation methods where explicit Bayesian analyses of
  the simulated data are performed, the prior parameters for these
  built-in studies are supplied in the `analysis_args` argument to
  `evsi()`. These assume Beta priors for probabilities, and Normal
  priors for the mean of a normal outcome.

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

- pars_datagen:

  Character vector identifying which columns of `inputs` are the
  parameters required to generate data from the proposed study. These
  should be columns of `inputs`.

  If `pars_datagen` is not supplied, then it is assumed to be the same
  as `pars`. Note that these can be different. Even if the study data
  are generated by a particular parameter, when analysing the data we
  could choose to ignore the information that the data provides about
  that parameter.

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

  Character string indicating the calculation method. Defaults to
  `"gam"`.

  All the nonparametric regression methods supported for
  [`evppi`](https://chjackson.github.io/voi/reference/evppi.md), that is
  `"gam","gp","earth","inla"`, can also be used for EVSI calculation by
  regressing on a summary statistic of the predicted data (Strong et al
  2015).

  `"is"` for importance sampling (Menzies 2016)

  `"mm"` for moment matching (Heath et al 2018)

  Note that the `"is"` and `"mm"` methods are used in conjunction with
  nonparametric regression, and the `gam_formula` argument can be
  supplied to `evsi` to specify this regression - see
  [`evppi`](https://chjackson.github.io/voi/reference/evppi.md) for
  documentation of this argument.

- likelihood:

  Likelihood function, required (and only required) for the importance
  sampling method when a study design other than one of the built-in
  ones is used. This should have two arguments, named as follows:

  `Y`: a one-row data frame of predicted data. Columns are defined by
  different outcomes in the data, with names matching the names of the
  data frame returned by `datagen_fn`.

  `inputs`. a data frame of simulated parameter values. Columns should
  correspond to different variables in `inputs`. The column names should
  all be found in the names of `inputs`, though they do not have to be
  in the same order, or include everything in `inputs`. The number or
  rows should be the same as the number of rows in `inputs`.

  The function should return a vector whose length matches the number of
  rows of the parameters data frame given as the second argument. Each
  element of the vector gives the likelihood of the corresponding set of
  parameters, given the data in the first argument. An example is given
  in the vignette.

  The likelihood can optionally have a `n` argument, which is
  interpreted as the sample size of the study. If the `n` argument to
  `evsi` is used then this is passed to the likelihood function.
  Conversely any `n` argument to `evsi` will be ignored by a likelihood
  function that does not have its own `n` argument.

  Note the definition of the likelihood should agree with the definition
  of `datagen_fn` to define a consistent sampling distribution for the
  data. No automatic check is performed for this.

- analysis_fn:

  Function which fits a Bayesian model to the generated data. Required
  for `method="mm"` if a study design other than one of the built-in
  ones is used. This should be a function that takes the following
  arguments:

  `data`: A data frame with names matching the output of `datagen_fn`

  `args`: A list with constants required in the Bayesian analysis, e.g.
  prior parameters, or options for the analysis, e.g. number of MCMC
  simulations. The component of this list called `n` is assumed to
  contain the sample size of the study.

  `pars` Names of the parameters whose posterior is being sampled.

  The function should return a data frame with names matching `pars`,
  containing a sample from the posterior distribution of the parameters
  given data supplied through `data`.

  `analysis_fn` is required to have all three of these arguments, but
  you do not need to use any elements of `args` or `pars` in the body of
  `analysis_fn`. Instead, sample sizes, prior parameters, MCMC options
  and parameter names can alternatively be hard-coded inside
  `analysis_fn`. Passing these through the function arguments (via the
  `analysis_args` argument to `evsi`) is only necessary if we want to
  use the same `analysis_fn` to do EVSI calculations with different
  sample sizes or other settings.

- analysis_args:

  List of arguments required for the Bayesian analysis of the predicted
  data, e.g. definitions of the prior and options to control sampling.
  Only used in `method="mm"`. This is required if the study design is
  one of the built-in ones specified in `study`. If a custom design is
  specifed through `analysis_fn`, then any constants needed in
  `analysis_fn` can either be supplied in `analysis_args`, or hard-coded
  in `analysis_fn` itself.

  For the built-in designs, the lists should have the following named
  components. An optional component `niter` in each case defines the
  posterior sample size (default 1000).

  `study="binary"`: `a` and `b`: Beta shape parameters

  `study="trial_binary"`: `a1` and `b1`: Beta shape parameters for the
  prior for the first arm, `a2` and `b2`: Beta shape parameters for the
  prior for the second arm.

  `study="normal_known"`: `prior_mean`, `prior_sd` (mean and standard
  deviation deviation of the Normal prior) and `sampling_sd` (SD of an
  individual-level normal observation, so that the sampling SD of the
  mean outcome over the study is `sampling_sd/sqrt(n)`.

- model_fn:

  Function which evaluates the decision-analytic model, given parameter
  values. Required for `method="mm"`. See
  [`evppi_mc`](https://chjackson.github.io/voi/reference/evppi_mc.md)
  for full documentation of the required specification of this function.

- par_fn:

  Function to simulate values from the uncertainty distributions of
  parameters needed by the decision-analytic model. Should take one
  argument and return a data frame with one row for each simulated
  value, and one column for each parameter. See
  [`evppi_mc`](https://chjackson.github.io/voi/reference/evppi_mc.md)
  for full specification.

- Q:

  Number of quantiles to use in `method="mm"`.

- npreg_method:

  Method to use to calculate the EVPPI, for those methods that require
  it. This is passed to
  [`evppi`](https://chjackson.github.io/voi/reference/evppi.md) as the
  `method` argument.

- nsim:

  Number of simulations from the model to use for calculating EVPPI. The
  first `nsim` rows of the objects in `inputs` and `outputs` are used.

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

  Other arguments understood by specific methods, e.g. `gam_formula` and
  other controlling options (see
  [`evppi`](https://chjackson.github.io/voi/reference/evppi.md)) can be
  passed to the nonparametric regression used inside the moment matching
  method.

## Value

A data frame with a column `pars`, indicating the parameter(s), and a
column `evsi`, giving the corresponding EVPPI. If the EVSI for multiple
sample sizes was requested, then the sample size is returned in the
column `n`, and if `outputs` is of "cost-effectiveness analysis" form,
so that there is one EVPPI per willingness-to-pay value, then a column
`k` identifies the willingness-to-pay.

## Details

See the [package overview / Get Started
vignette](https://chjackson.github.io/voi/articles/voi.html#evsi) for
some examples of using this function.

## References

Heath, A., Kunst, N., & Jackson, C. (eds.). (2024). Value of Information
for Healthcare Decision-Making. CRC Press.

Strong, M., Oakley, J. E., Brennan, A., & Breeze, P. (2015). Estimating
the expected value of sample information using the probabilistic
sensitivity analysis sample: a fast, nonparametric regression-based
method. Medical Decision Making, 35(5), 570-583.

Menzies, N. A. (2016). An efficient estimator for the expected value of
sample information. Medical Decision Making, 36(3), 308-320.

Heath, A., Manolopoulou, I., & Baio, G. (2018). Efficient Monte Carlo
estimation of the expected value of sample information using moment
matching. Medical Decision Making, 38(2), 163-173.

Heath, A., Manolopoulou, I., & Baio, G. (2019). Estimating the expected
value of sample information across different sample sizes using moment
matching and nonlinear regression. Medical Decision Making, 39(4),
347-359.
