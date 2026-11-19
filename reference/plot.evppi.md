# Plot EVPPI estimates

Plot EVPPI estimates as simple dot or curve plots.

## Usage

``` r
# S3 method for class 'evppi'
plot(x, type = NULL, order = FALSE, top = NULL, ...)
```

## Arguments

- x:

  Object returned from
  [`evppi`](https://chjackson.github.io/voi/reference/evppi.md).

- type:

  `"dots"` for a dot plot of the EVPPI by parameter. If `x` includes
  multiple willingness-to-pay values for the same parameter, these are
  shown as multiple dots.

  `"curves"` for a plot of EVPPI against willingness-to-pay, with
  different parameters distinguished as different curves. This is only
  applicable if there are multiple willingness-to-pay values included in
  `x`.

- order:

  For dot plots, order the plot with highest EVPPI values at the top.

- top:

  A positive integer. If specified, for example as `top=5` then only
  five parameters are included in the plot, those with the top five
  maximum EVPPI values by parameter.

- ...:

  Other arguments (currently unused).

## Value

A `ggplot2` object.

## Details

These plotting functions are intended for quick interactive exploration
of EVPPI results, so they deliberately have limited options for
customising them.

For publication quality graphics, it is advised to use `ggplot2` by hand
on the data returned by `evppi`. Examine the code for `plot_evppi_dots`
and `plot_evppi_curves` to see how these plots might be constructed.
