# Population expected value of information

Convert per-person expected value of information to the population
expected value of information, given a discount rate over some time
horizon.

## Usage

``` r
pop_voi(voi, pop, time, dis = 0.035)
```

## Arguments

- voi:

  Vector of estimates of any per-person value of information measure,
  e.g. the `evsi` column of the data frame returned by
  [`evsi`](https://chjackson.github.io/voi/reference/evsi.md) or the
  correspondingly-named columns of the data frames returned by
  [`evppi`](https://chjackson.github.io/voi/reference/evppi.md) or
  [`evpi`](https://chjackson.github.io/voi/reference/evpi.md).

- pop:

  Size of the population who would be affected by the decision.

- time:

  Time horizon over which discounting will be applied.

- dis:

  Discount rate used when converting per-person to population EVSI.

## Value

A vector of population VoI estimates.

## Details

Calculated as `voi*pop/dis*(1 - exp(-dis*time))`, or `voi*pop` if the
discount rate is zero. This is a continuous-time variant of the typical
discrete-time discounting formula.

Any arguments may be supplied as vectors, in which case, all arguments
are replicated to the length of the longest argument.
