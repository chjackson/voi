# Check that a decision-analytic model function is of the appropriate form. Detect if it returns net benefit or CEA format note check_outputs adds a class

Check that a decision-analytic model function is of the appropriate
form. Detect if it returns net benefit or CEA format note check_outputs
adds a class

## Usage

``` r
check_model_fn(
  model_fn,
  par_fn,
  mfargs = NULL,
  outputs = NULL,
  verbose = FALSE
)
```

## Value

A modified copy of `model_fn` with a `class` attribute indicating
whether it is in net benefit `"nb"` or cost-effectiveness `"cea"`
format, and an attribute `"nopt"` giving the number of decision options.
