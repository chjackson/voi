# Expected net benefit of sampling

Calculates the expected net benefit of sampling for a typical study to
inform a health economic evaluation, given estimates of the per-person
expected value of sample information, decision population size and study
setup and per-participant costs. The optimal sample size for each
willingness-to-pay, population size and time horizon is also determined.

## Usage

``` r
enbs(
  evsi,
  costs_setup,
  costs_pp,
  pop,
  time,
  dis = 0.035,
  smooth = FALSE,
  smooth_df = NULL,
  pcut = 0.05
)
```

## Arguments

- evsi:

  Data frame giving estimates of the expected value of sample
  information, as returned by
  [`evsi`](https://chjackson.github.io/voi/reference/evsi.md). This may
  contain multiple estimates, one for each sample size and willingness
  to pay.

- costs_setup:

  Setup costs of the study. This can either be a constant, or a vector
  of two elements giving a 95% credible interval (with mean defined by
  the midpoint), or a vector of three elements assumed to define the
  mean and 95% credible interval.

- costs_pp:

  Per-participant costs of the study, supplied in the same format as
  `cost_setup`.

- pop:

  Size of the population who would be affected by the decision.

- time:

  Time horizon over which discounting will be applied.

- dis:

  Discount rate used when converting per-person to population EVSI.

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

- pcut:

  Cut-off probability which defines a "near-optimal" sample size. The
  minimum and maximum sample size for which the ENBS is within `pcut`
  (by default 5%) of its maximum value will be determined.

## Value

Data frame with components `enbs` giving the ENBS, and `sd` giving the
corresponding standard deviation. The rows of the data frame correspond
to the rows of `evsi`, and any `n` and `k` are inherited from `evsi`.
Additional columns include:

`pce`: the probability that the study is cost-effective, i.e. that the
ENBS is positive, obtained from a normal distribution defined by the
estimate and standard deviation.

`enbsmax`: The maximum ENBS for each willingness-to-pay `k`.

`nmax`: The sample size `n` at which this maximum is achieved.

A second data frame is returned as the `"enbsmax"` attribute. This has
one row per willingness-to-pay (`k`), giving the optimal ENBS
(`enbsmax`) the optimal sample size (`nmax`) and an interval estimate
for the optimal sample size (`nlower` to `nupper`).

If `pop`, `time` or `dis` were supplied as vectors of more than one
element, then additional columns will be returned in these data frames
to identify the population, time or discount rate for each ENBS
calculation. An index `ind` is also returned to identify the unique
combination that each row refers to.

## Details

`pop`,`time` and `dis` may be supplied as vectors of different lengths.
In that case, the ENBS is calculated for all possible combinations of
the values in these vectors.

## References

Value of Information for Healthcare Decision Making (CRC Press, eds.
Heath, Kunst and Jackson: forthcoming)
