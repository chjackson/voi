# Determine the optimum sample size in an analysis of the expected net benefit of sampling

The optimum sample size for a given willingness to pay is determined
either by a simple search over the supplied ENBS estimates for different
sample sizes, or by a regression and interpolation method.

## Usage

``` r
enbs_opt(x, pcut = 0.05, smooth = FALSE, smooth_df = NULL, keep_preds = FALSE)
```

## Arguments

- x:

  Data frame containing a set of ENBS estimates for different sample
  sizes, which will be optimised over. Usually this is for a common
  willingness-to-pay. The required components are `enbs` and `n`.

- pcut:

  Cut-off probability which defines a "near-optimal" sample size. The
  minimum and maximum sample size for which the ENBS is within `pcut`
  (by default 5%) of its maximum value will be determined.

- smooth:

  If `TRUE`, then the maximum ENBS is determined after fitting a
  nonparametric regression to the data frame `x`, which estimates and
  smooths the ENBS for every integer sample size in the range of `x$n`.
  The regression is done using the default settings of
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) from the mgcv package.

  If this is `FALSE`, then no smoothing or interpolation is done, and
  the maximum is determined by searching over the values supplied in
  `x`.

- smooth_df:

  Basis dimension for the smooth regression. Passed as the `k` argument
  to the `s()` term in [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).
  Defaults to 6, or the number of unique sample sizes minus 1 if this is
  lower. Set to a higher number if you think the smoother does not
  capture the relation of ENBS to sample size accurately enough.

- keep_preds:

  If `TRUE` and `smooth=TRUE` then the data frame of predictions from
  the smooth regression model is stored in the `"preds"` attribute of
  the result.

## Value

A data frame with one row, and the following columns:

`ind`: An integer index identifying, e.g. the willingness to pay and
other common characteristics of the ENBS estimates (e.g. incident
population size, decision time horizon). This is copied from `x$ind`.

`enbsmax`: the maximum ENBS

`nmax`: the sample size at which this maximum is achieved

`nlower`: the lowest sample size for which the ENBS is within

`pcut` (default 5%) of its maximum value

`nupper`: the corresponding highest ENBS
