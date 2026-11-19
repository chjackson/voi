# Package index

## Package overview

- [`voi`](https://chjackson.github.io/voi/reference/voi-package.md)
  [`voi-package`](https://chjackson.github.io/voi/reference/voi-package.md)
  : Methods to calculate the Expected Value of Information

## Calculating the expected value of information

- [`evpi()`](https://chjackson.github.io/voi/reference/evpi.md) :
  Calculate the expected value of perfect information from a decision
  model
- [`evppi()`](https://chjackson.github.io/voi/reference/evppi.md) :
  Calculate the expected value of partial perfect information from a
  decision-analytic model
- [`evppi_mc()`](https://chjackson.github.io/voi/reference/evppi_mc.md)
  : Traditional two-level Monte Carlo estimator of EVPPI.
- [`evppivar()`](https://chjackson.github.io/voi/reference/evppivar.md)
  : Calculate the expected value of partial perfect information for an
  estimation problem
- [`plot(`*`<evppi>`*`)`](https://chjackson.github.io/voi/reference/plot.evppi.md)
  : Plot EVPPI estimates
- [`evsi()`](https://chjackson.github.io/voi/reference/evsi.md) :
  Calculate the expected value of sample information from a
  decision-analytic model
- [`evsivar()`](https://chjackson.github.io/voi/reference/evsivar.md) :
  Calculate the expected value of sample information for an estimation
  problem
- [`enbs()`](https://chjackson.github.io/voi/reference/enbs.md) :
  Expected net benefit of sampling
- [`enbs_opt()`](https://chjackson.github.io/voi/reference/enbs_opt.md)
  : Determine the optimum sample size in an analysis of the expected net
  benefit of sampling
- [`pop_voi()`](https://chjackson.github.io/voi/reference/pop_voi.md) :
  Population expected value of information

## Utilities

- [`import_heemod_outputs()`](https://chjackson.github.io/voi/reference/import_heemod.md)
  [`import_heemod_inputs()`](https://chjackson.github.io/voi/reference/import_heemod.md)
  : Import results of probabilistic analysis from heemod
- [`all_interactions()`](https://chjackson.github.io/voi/reference/all_interactions.md)
  : Generate a string with all interactions of a certain degree, to be
  used in a GAM formula
- [`check_regression()`](https://chjackson.github.io/voi/reference/check_regression.md)
  : Check the fit of a regression model used to estimate EVPPI or EVSI

## Example decision models

- [`chemo_cea`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_nb`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_pars`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_cea_501`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_constants`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_evsi_or`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_pars_fn()`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_model_nb()`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_model_cea()`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_model_lor_nb()`](https://chjackson.github.io/voi/reference/chemo_model.md)
  [`chemo_model_lor_cea()`](https://chjackson.github.io/voi/reference/chemo_model.md)
  : Chemotherapy cost-effectiveness model
