##' Chemotherapy cost-effectiveness model 
##'
##' An artificial health economic decision model with a typical Markov model structure, used for illustrating Value of Information methods.
##' Functions are provided to generate model parameters and evaluate the model, and samples from probabilistic analysis of the model are
##' provided as built-in datasets.
##'
##' @param p_side_effects_t1 Probability of side effects under treatment 1
##' @param p_side_effects_t2 Probability of side effects under treatment 2
##' @param logor_side_effects Log odds ratio for side effects between treatment 1 and treatment 2
##' @param p_home_home Annual transition probability from home to home
##' @param p_home_hospital Transition probability from home to hospital
##' @param p_home_recover Transition probability from home to recovery
##' @param p_hospital_hospital Transition probability from hospital to hospital
##' @param p_hospital_recover Transition probability from hospital to recovery
##' @param p_hospital_dead Transition probability from hospital to death
##' @param c_home_care Cost of a yearly period under treatment at home 
##' @param c_hospital Cost of hospital treatment 
##' @param c_death Cost of death 
##' @param u_recovery Utility of a period in the recovery state 
##' @param u_home_care Utility of home care state 
##' @param u_hospital Utility of hospital state
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
##' * `logor_side_effects`: log odds ratio of side effects for treatment 2 compared to 1
##' 
##' * `p_hospitalised_total`: probability of hospitalisation over the 50 year time horizon
##' 
##' * `p_died`: probability of death over the time horizon, given hospitalisation
##' 
##' * `lambda_home`: conditional probability that a patient recovers given they are not hospitalised
##' 
##' * `lambda_hosp`: conditional probability that a patient in hospital recovers given they do not die
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
chemo_model_nb <- function(p_side_effects_t1, 
                           p_side_effects_t2,
                           p_home_home, p_home_hospital, p_home_recover,
                           p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                           c_home_care, c_hospital, c_death,
                           u_recovery, u_home_care, u_hospital)
{
    if (length(p_side_effects_t1) > 1)
        stop("This function is not vectorised, and parameters should be supplied as scalars")
    ce <- calculate_costs_effects(p_side_effects_t1, p_side_effects_t2,
                                  p_home_home, p_home_hospital, p_home_recover,
                                  p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                                  c_home_care, c_hospital, c_death,
                                  u_recovery, u_home_care, u_hospital)
    calculate_net_benefit(ce, wtp=20000)
}

##' @rdname chemo_model
##' @export
chemo_model_cea <- function(p_side_effects_t1, 
                            p_side_effects_t2,
                            p_home_home, p_home_hospital, p_home_recover,
                            p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                            c_home_care, c_hospital, c_death,
                            u_recovery, u_home_care, u_hospital)
{
    if (length(p_side_effects_t1) > 1)
        stop("This function is not vectorised, and parameters should be supplied as scalars")
    ce <- calculate_costs_effects(p_side_effects_t1, p_side_effects_t2,
                                  p_home_home, p_home_hospital, p_home_recover,
                                  p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                                  c_home_care, c_hospital, c_death,
                                  u_recovery, u_home_care, u_hospital)
    t(ce[1,,])
}

##' @rdname chemo_model
##' @export
chemo_model_lor_nb <- function(p_side_effects_t1, logor_side_effects,
                               p_home_home, p_home_hospital, p_home_recover,
                               p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                               c_home_care, c_hospital, c_death,
                               u_recovery, u_home_care, u_hospital)
{
  if (length(p_side_effects_t1) > 1)
    stop("This function is not vectorised, and parameters should be supplied as scalars")
  odds1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  odds2 <- odds1 * exp(logor_side_effects)
  p_side_effects_t2 <- odds2 / (1 + odds2)
  ce <- calculate_costs_effects(p_side_effects_t1, p_side_effects_t2,
                                p_home_home, p_home_hospital, p_home_recover,
                                p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                                c_home_care, c_hospital, c_death,
                                u_recovery, u_home_care, u_hospital)
  calculate_net_benefit(ce, wtp=20000)
}

##' @rdname chemo_model
##' @export
chemo_model_lor_cea <- function(p_side_effects_t1, logor_side_effects,
                                p_home_home, p_home_hospital, p_home_recover,
                                p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                                c_home_care, c_hospital, c_death,
                                u_recovery, u_home_care, u_hospital)
{
  if (length(p_side_effects_t1) > 1)
    stop("This function is not vectorised, and parameters should be supplied as scalars")
  odds1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  odds2 <- odds1 * exp(logor_side_effects)
  p_side_effects_t2 <- odds2 / (1 + odds2)
  ce <- calculate_costs_effects(p_side_effects_t1, p_side_effects_t2,
                                p_home_home, p_home_hospital, p_home_recover,
                                p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                                c_home_care, c_hospital, c_death,
                                u_recovery, u_home_care, u_hospital)
  t(ce[1,,])
}
