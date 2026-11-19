# Methods to calculate the Expected Value of Information

[`evppi`](https://chjackson.github.io/voi/reference/evppi.md) calculates
the expected value of partial perfect information from a
decision-analytic model. The default, recommended computation methods
are based on nonparametric regression.
[`evpi`](https://chjackson.github.io/voi/reference/evpi.md) is also
provided for the expected value of perfect information.

[`evsi`](https://chjackson.github.io/voi/reference/evsi.md) calculates
the expected value of sample information. Currently this implements the
same set of nonparametric regression methods as in
[`evppi`](https://chjackson.github.io/voi/reference/evppi.md), and
methods based on moment matching and importance sampling.
[`enbs`](https://chjackson.github.io/voi/reference/enbs.md) can then be
used to calculate and optimise the expected net benefit of sampling for
a simple study with a fixed upfront cost and per-participant costs.

[`evppi`](https://chjackson.github.io/voi/reference/evppi.md) and
[`evsi`](https://chjackson.github.io/voi/reference/evsi.md) both require
a sample of inputs and outputs from a Monte Carlo probabilistic analysis
of a decision-analytic model.

Analogous functions
[`evppivar`](https://chjackson.github.io/voi/reference/evppivar.md) and
[`evsivar`](https://chjackson.github.io/voi/reference/evsivar.md)
calculate the EVPPI and EVSI for models used for estimation rather than
decision-making. The value of information is measured by expected
reductions in variance of an uncertain model output of interest.

A pure "brute-force" Monte Carlo method for EVPPI calculation is
provided in
[`evppi_mc`](https://chjackson.github.io/voi/reference/evppi_mc.md),
though this is usually computationally impractical.

The [package overview / Get Started
vignette](https://chjackson.github.io/voi/articles/voi.html) gives
worked examples of the use of all of these functions.

## References

Heath, A., Kunst, N., & Jackson, C. (eds.). (2024). Value of Information
for Healthcare Decision-Making. CRC Press.

Heath, A., Manolopoulou, I., & Baio, G. (2017). A review of methods for
analysis of the expected value of information. Medical Decision Making,
37(7), 747-758.

Heath, A., Kunst, N., Jackson, C., Strong, M., Alarid-Escudero, F.,
Goldhaber-Fiebert, J. D., Baio, G. Menzies, N.A, Jalal, H. (2020).
Calculating the Expected Value of Sample Information in Practice:
Considerations from 3 Case Studies. Medical Decision Making, 40(3),
314-326.

Kunst, N., Wilson, E. C., Glynn, D., Alarid-Escudero, F., Baio, G.,
Brennan, A., Fairley, M., Glynn, D., Goldhaber-Fiebert, J. D., Jackson,
C., Jalal, H., Menzies, N. A., Strong, M., Thom, H., Heath, A. (2020).
Computing the Expected Value of Sample Information Efficiently:
Practical Guidance and Recommendations for Four Model-Based Methods.
Value in Health, 3(6), 734-742.

## See also

Useful links:

- <https://chjackson.github.io/voi/>

- Report bugs at <https://github.com/chjackson/voi/issues>

## Author

**Maintainer**: Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

Authors:

- Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>

- Anna Heath <anna.heath@sickkids.ca>

Other contributors:

- Gianluca Baio <g.baio@ucl.ac.uk> (Author of code taken from the BCEA
  package) \[contributor\]

- Mark Strong <mark.strong@sheffield.ac.uk> (Author of code taken from
  the SAVI package) \[contributor\]

- Kofi Placid Adragni (Author of code taken from the ldr package)
  \[contributor\]

- Andrew Raim (Author of code taken from the ldr package)
  \[contributor\]
