# Import results of probabilistic analysis from heemod

[heemod](https://CRAN.R-project.org/package=heemod) is a package for
constructing common forms of health economic decision models. The
outputs from probabilistic analysis of these models can be imported
using these functions, to allow Value of Information measures to be
calculated for them using the voi package.

## Usage

``` r
import_heemod_outputs(obj, k = NULL)

import_heemod_inputs(obj)
```

## Arguments

- obj:

  Object returned by the
  [`run_psa`](https://aphp.github.io/heemod/reference/run_psa.html)
  function in heemod, containing samples from probabilistic analysis of
  a decision model.

- k:

  Vector of willingness-to-pay values. The default is inherited from the
  [`bcea`](https://n8thangreen.github.io/BCEA/reference/bcea.html)
  function from the BCEA package.

## Value

`import_heemod_outputs` produces a list of model outputs in
"cost-effectiveness analysis" format, that can be supplied as the
`outputs` argument to
[`evppi`](https://chjackson.github.io/voi/reference/evppi.md) and
similar functions in the voi package. Both the heemod and BCEA packages
need to be installed to use this.

`import_heemod_inputs` produces a data frame with samples of parameter
values under uncertainty, that can be supplied as the `inputs` argument
to [`evppi`](https://chjackson.github.io/voi/reference/evppi.md) and
similar functions in voi.
