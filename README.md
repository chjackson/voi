# voi: a generic package to calculate the expected value of information

* A common interface for several methods to calculate the

  - Expected Value of (Partial) Perfect Information and the 

  - Expected Value of Sample Information 

* A project of the [ConVOI Group: the Collaborative Network for Value of Information](https://www.convoi-group.org/)


## Comparison with other packages

`voi` is pure "command-based" R, with no web interface like [SAVI](https://github.com/Sheffield-Accelerated-VoI/SAVI).

* The R commands in `voi` are clean and consistent: they all have the same basic interface, so you can switch between computational methods easily.

* Outputs are all in "tidy" data frames for consistency, and to facilitate post-processing and plotting with modern tools such as [ggplot2](https://ggplot2.tidyverse.org/).

#### EVPPI computation

* `voi` includes all the [EVPPI computation](https://chjackson.github.io/voi/articles/voi.html#evppi) methods that are in [SAVI](https://github.com/Sheffield-Accelerated-VoI/SAVI) (GAM and Gaussian process regression), and includes the INLA method from [BCEA](https://cran.r-project.org/package=BCEA).

* Some other nonparametric regression methods ([`"earth"`](https://chjackson.github.io/voi/articles/voi.html#earth), [`"bart"`](https://chjackson.github.io/voi/articles/voi.html#bart)) are included in `voi`, which may improve efficiency for multiparameter EVPPI computation problems with large numbers of parameters.

#### EVSI computation

* `voi` is the first package to implement a range of [EVSI computation](https://chjackson.github.io/voi/articles/voi.html#evsi) methods: nonparametric regression, moment matching and importance sampling.   A simple model for the [expected net benefit of sampling](https://chjackson.github.io/voi/articles/plots.html) is also included.

#### In summary 

* `voi` will not benefit you if you want a web interface, or if you just need single-parameter EVPPI and are happy with SAVI/BCEA.

* `voi` will benefit you if you want to calculate EVSI, or multiparameter EVPPI with large numbers of parameters. 








## Installation

```
remotes::install_github("chjackson/voi")
 ```

## User guide 

[`voi` for Value of Information calculation: package overview](https://chjackson.github.io/voi/articles/voi.html)


## Source code

[Github repository](https://github.com/chjackson/voi)


<!-- badges: start -->
[![R-CMD-check](https://github.com/chjackson/voi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chjackson/voi/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/chjackson/voi/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/chjackson/voi/actions/workflows/test-coverage.yaml)
[![DOI](https://zenodo.org/badge/227181181.svg)](https://zenodo.org/badge/latestdoi/227181181)
<!-- badges: end -->
