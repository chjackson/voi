##' Chemotherapy cost-effectiveness model 
##'
##' An artificial health economic decision model with a typical Markov model structure, used for illustrating Value of Information methods.
##' Functions are provided to generate model parameters and evaluate the model, and samples from probabilistic analysis of the model are
##' provided as built-in datasets.
##'
##' TODO refer to VoI book for more details
##'
##' @param p_side_effects_t1 Probability of side effects under treatment 1
##' @param p_side_effects_t2 Probability of side effects under treatment 2. (this is unused - TODO remove it) 
##' @param logor_side_effects Log odds ratio of side effects for treatment 2 compared to 1
##' @param p_hospitalised_total Probability of hospitalisation in the year after receiving treatment
##' @param p_died Probability of death in the year after receiving treatment
##' @param lambda_home Recovery probability for someone treated at home
##' @param lambda_hosp Recovery probability for someone treated in hospital who does not die
##' @param c_home_care Cost of a yearly period under treatment at home 
##' @param c_hospital Cost of hospital treatment 
##' @param c_death Cost of death 
##' @param u_recovery Utility of a period in the recovery state 
##' @param u_home_care Utility of home care state 
##' @param u_hospital Utility of hospital state
##' @param rate_longterm Long term mortality rate
##'
##' @param n Number of samples to generate from the uncertainty distribution of the parameters in \code{chemo_pars_fn}. 
##' 
##' @return
##' Two alternative functions are provided to evaluate the decision model for given parameters.
##' 
##' \code{chemo_model_nb} returns a vector with elements giving the net monetary benefit for standard of care
##' and novel treatment, respectively, at a willingness-to-pay of 20,000 pounds per QALY.
##'
##' \code{chemo_model_cea} returns a matrix with:
##'
##' * two rows, the first for expected costs and the second for expected effects (QALYs) over the fifty year time horizon, and 
##' 
##' * two columns, the first for the "standard of care" decision option, and the second for the novel
##' treatment.
##'
##' \code{chemo_model_lor_nb} and \code{chemo_model_lor_cea} are the same model, but parameterised in terms of
##' the probability of side effects for the standard of care \code{p_side_effects_t1} and the log odds ratio
##' of side effects between treatment groups \code{logor_side_effects}, rather than in terms of
##' \code{p_side_effects_t1} and \code{p_side_effects_t2}
##' 
##' \code{chemo_pars_fn} generates a sample from the uncertainty distribution of the parameters in the chemotherapy model . This returns a data frame with parameters matching the arguments of
##' \code{\link{chemo_model_nb}}, and the following additional derived parameters:
##'
##' * `p_side_effects_t2`: 
##' 
##' * `p_hospitalised_total`: probability of hospitalisation over the 50 year time horizon
##' 
##' * `p_died`: probability of death over the time horizon, given hospitalisation
##' 
##' * `lambda_home`: conditional probability that a patient recovers given they are not hospitalised
##' 
##' * `lambda_hosp`: conditional probability that a patient in hospital recovers given they do not die
##'
##' \code{chemo_constants} includes various constants required by the code.
##'
##' @format Samples of 10000 from probabilistic analysis of this model are made available in the package, in the
##' following data objects:
##' 
##' \code{chemo_pars}: Sample from the distributions of the parameters, as a data frame with names as documented above.
##'
##' \code{chemo_cea}: List with components `e` (sampled effects), `c` (sampled costs), and `k` (a set of five
##' equally-spaced willingess-to-pay values from 10000 to 50000 pounds).   The effects and costs are data frames
##' with two columns, one for each decision option. 
##'
##' \code{chemo_nb}: Data frame with two columns, giving the net monetary benefit for each decision option,
##' at a willingness-to-pay of 20000 pounds. 
##' 
##' @name chemo_model
NULL

##' @rdname chemo_model
##' @export
chemo_pars_fn <- function(n){
    generate_psa_parameters(n)
}

##' @rdname chemo_model
##' @export
chemo_model_nb <- function(p_side_effects_t1, p_side_effects_t2,
                           p_hospitalised_total, p_died,
                           lambda_home, lambda_hosp,
                           c_home_care, c_hospital, c_death,
                           u_recovery, u_home_care, u_hospital,
                           rate_longterm)
{
    if (length(p_side_effects_t1) > 1)
        stop("This function is not vectorised, and parameters should be supplied as scalars")

    ce <- chemo_model_cea(p_side_effects_t1 = p_side_effects_t1,
                          p_side_effects_t2 = p_side_effects_t2,
                          p_hospitalised_total = p_hospitalised_total,
                          p_died = p_died,
                          lambda_home = lambda_home,
                          lambda_hosp = lambda_hosp,
                          c_home_care = c_home_care,
                          c_hospital = c_hospital,
                          c_death = c_death,
                          u_recovery = u_recovery,
                          u_home_care = u_home_care,
                          u_hospital = u_hospital,
                          rate_longterm = rate_longterm)
    ce[1,]*20000 - ce[2,]
}

##' @rdname chemo_model
##' @export
chemo_model_cea <- function(p_side_effects_t1, p_side_effects_t2,
                            p_hospitalised_total, p_died,
                            lambda_home, lambda_hosp,
                            c_home_care, c_hospital, c_death,
                            u_recovery, u_home_care, u_hospital,
                            rate_longterm)
{
  if (length(p_side_effects_t1) > 1)
    stop("This function is not vectorised, and parameters should be supplied as scalars")

    odds2 <- p_side_effects_t2 / (1 - p_side_effects_t2)
    odds1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
    logor_side_effects <- log(odds2 / odds1)
  
    ce <- calculate_costs_effects(p_side_effects_t1 = p_side_effects_t1,
                                  p_hospitalised_total = p_hospitalised_total,
                                  p_died = p_died,
                                  lambda_home = lambda_home,
                                  lambda_hosp = lambda_hosp,
                                  c_home_care = c_home_care,
                                  c_hospital = c_hospital,
                                  c_death = c_death,
                                  u_recovery = u_recovery,
                                  u_home_care = u_home_care,
                                  u_hospital = u_hospital,
                                  logor_side_effects = logor_side_effects,
                                  rate_longterm = rate_longterm)
  ce
}

##' @rdname chemo_model
##' @export
chemo_model_lor_nb <- function(p_side_effects_t1, logor_side_effects, 
                               p_hospitalised_total, p_died,
                               lambda_home, lambda_hosp,
                               c_home_care, c_hospital, c_death,
                               u_recovery, u_home_care, u_hospital,
                               rate_longterm)
{
  if (length(p_side_effects_t1) > 1)
    stop("This function is not vectorised, and parameters should be supplied as scalars")

  ce <- chemo_model_lor_cea(p_side_effects_t1 = p_side_effects_t1,
                            logor_side_effects = logor_side_effects,
                            p_hospitalised_total = p_hospitalised_total,
                            p_died = p_died,
                            lambda_home = lambda_home,
                            lambda_hosp = lambda_hosp,
                            c_home_care = c_home_care,
                            c_hospital = c_hospital,
                            c_death = c_death,
                            u_recovery = u_recovery,
                            u_home_care = u_home_care,
                            u_hospital = u_hospital,
                            rate_longterm = rate_longterm)
  ce[1,]*20000 - ce[2,]
}

##' @rdname chemo_model
##' @export
chemo_model_lor_cea <- function(p_side_effects_t1, logor_side_effects, 
                                p_hospitalised_total, p_died,
                                lambda_home, lambda_hosp,
                                c_home_care, c_hospital, c_death,
                                u_recovery, u_home_care, u_hospital,
                                rate_longterm)
{
  if (length(p_side_effects_t1) > 1)
    stop("This function is not vectorised, and parameters should be supplied as scalars")
  odds1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  odds2 <- odds1 * exp(logor_side_effects)
  p_side_effects_t2 <- odds2 / (1 + odds2)

  ce <- chemo_model_cea(p_side_effects_t1 = p_side_effects_t1,
                        p_side_effects_t2 = p_side_effects_t2,
                        p_hospitalised_total = p_hospitalised_total,
                        p_died = p_died,
                        lambda_home = lambda_home,
                        lambda_hosp = lambda_hosp,
                        c_home_care = c_home_care,
                        c_hospital = c_hospital,
                        c_death = c_death,
                        u_recovery = u_recovery,
                        u_home_care = u_home_care,
                        u_hospital = u_hospital,
                        rate_longterm = rate_longterm)

  ce
}

##' @rdname chemo_model
"chemo_constants"
