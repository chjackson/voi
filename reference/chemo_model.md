# Chemotherapy cost-effectiveness model

An artificial health economic decision model with a typical Markov model
structure, used for illustrating Value of Information methods. Functions
are provided to generate model parameters and evaluate the model, and
samples from probabilistic analysis of the model are provided as
built-in datasets.

## Usage

``` r
chemo_cea

chemo_nb

chemo_pars

chemo_cea_501

chemo_constants

chemo_evsi_or

chemo_pars_fn(n)

chemo_model_nb(
  p_side_effects_t1,
  p_side_effects_t2,
  p_hospitalised_total,
  p_died,
  lambda_home,
  lambda_hosp,
  c_home_care,
  c_hospital,
  c_death,
  u_recovery,
  u_home_care,
  u_hospital,
  rate_longterm
)

chemo_model_cea(
  p_side_effects_t1,
  p_side_effects_t2,
  p_hospitalised_total,
  p_died,
  lambda_home,
  lambda_hosp,
  c_home_care,
  c_hospital,
  c_death,
  u_recovery,
  u_home_care,
  u_hospital,
  rate_longterm
)

chemo_model_lor_nb(
  p_side_effects_t1,
  logor_side_effects,
  p_hospitalised_total,
  p_died,
  lambda_home,
  lambda_hosp,
  c_home_care,
  c_hospital,
  c_death,
  u_recovery,
  u_home_care,
  u_hospital,
  rate_longterm
)

chemo_model_lor_cea(
  p_side_effects_t1,
  logor_side_effects,
  p_hospitalised_total,
  p_died,
  lambda_home,
  lambda_hosp,
  c_home_care,
  c_hospital,
  c_death,
  u_recovery,
  u_home_care,
  u_hospital,
  rate_longterm
)
```

## Format

An object of class `list` of length 33.

An object of class `evsi` (inherits from `data.frame`) with 15030 rows
and 3 columns.

Samples of 10000 from probabilistic analysis of this model are made
available in the package, in the following data objects:

`chemo_pars`: Sample from the distributions of the parameters, as a data
frame with names as documented above.

`chemo_cea`: List with components `e` (sampled effects), `c` (sampled
costs), and `k` (a set of five equally-spaced willingess-to-pay values
from 10000 to 50000 pounds). The effects and costs are data frames with
two columns, one for each decision option.

`chemo_nb`: Data frame with two columns, giving the net monetary benefit
for each decision option, at a willingness-to-pay of 20000 pounds.

`chemo_cea_501`: List with components `e` (sampled effects), `c`
(sampled costs), and `k` (a set of 501 willlingess-to-pay values from
10000 to 50000) This is provided to facilitate illustrations of plots of
VoI measures against willingness-to-pay.

The following additional data objects are supplied:

`chemo_constants` includes various constants required by the code.

`chemo_evsi_or` is the result of an EVSI analysis to estimate the
expected value of a two-arm trial, with a binary outcome, to estimate
the log odds ratio of side effects. This object is a data frame with
three columns, giving the sample size per arm (`n`), willingness-to-pay
(`k`) and the corresponding EVSI (`evsi`).

## Arguments

- n:

  Number of samples to generate from the uncertainty distribution of the
  parameters in `chemo_pars_fn`.

- p_side_effects_t1:

  Probability of side effects under treatment 1

- p_side_effects_t2:

  Probability of side effects under treatment 2

- p_hospitalised_total:

  Probability of hospitalisation in the year after receiving treatment

- p_died:

  Probability of death in the year after receiving treatment

- lambda_home:

  Recovery probability for someone treated at home

- lambda_hosp:

  Recovery probability for someone treated in hospital who does not die

- c_home_care:

  Cost of a yearly period under treatment at home

- c_hospital:

  Cost of hospital treatment

- c_death:

  Cost of death

- u_recovery:

  Utility of a period in the recovery state

- u_home_care:

  Utility of home care state

- u_hospital:

  Utility of hospital state

- rate_longterm:

  Long term mortality rate

- logor_side_effects:

  Log odds ratio of side effects for treatment 2 compared to 1

## Value

Two alternative functions are provided to evaluate the decision model
for given parameters.

`chemo_model_nb` returns a vector with elements giving the net monetary
benefit for standard of care and novel treatment, respectively, at a
willingness-to-pay of 20,000 pounds per QALY.

`chemo_model_cea` returns a matrix with:

- two rows, the first for expected costs and the second for expected
  effects (QALYs) over the fifty year time horizon, and

- two columns, the first for the "standard of care" decision option, and
  the second for the novel treatment.

`chemo_model_lor_nb` and `chemo_model_lor_cea` are the same model, but
parameterised in terms of the probability of side effects for the
standard of care `p_side_effects_t1` and the log odds ratio of side
effects between treatment groups `logor_side_effects`, rather than in
terms of `p_side_effects_t1` and `p_side_effects_t2`

`chemo_pars_fn` generates a sample from the uncertainty distribution of
the parameters in the chemotherapy model . This returns a data frame
with parameters matching the arguments of `chemo_model_nb`, and the
following additional derived parameters:

- `p_side_effects_t2`:

- `p_hospitalised_total`: probability of hospitalisation over the 50
  year time horizon

- `p_died`: probability of death over the time horizon, given
  hospitalisation

- `lambda_home`: conditional probability that a patient recovers given
  they are not hospitalised

- `lambda_hosp`: conditional probability that a patient in hospital
  recovers given they do not die

## Details

For more details, refer to Heath et al. (forthcoming book...)

## References

Value of Information for Healthcare Decision Making (CRC Press, eds.
Heath, Kunst and Jackson: forthcoming)
