##' Chemotherapy cost-effectiveness model 
##'
##' An artificial health economic decision model with a typical Markov model structure, used for illustrating Value of Information methods. 
##'
##' @param p_side_effects_t1 Probability of side effects under treatment 1
##' @param p_side_effects_t2 Probability of side effects under treatment 2
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
##' @return
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
##' \code{chemo_pars_fn} generates a sample from the uncertainty distribution of the parameters in the chemotherapy model . This returns a data frame with parameters matching the arguments of
##' \code{\link{chemo_model_nb}}.  Additional derived parameters are: 
##'
##' * `rr_side_effects`: relative risk of side effects for treatment 2 compared to 1
##' 
##' * `p_hospitalised_total`: probability of hospitalisation over the 50 year time horizon
##' 
##' * `p_died`: probability of death over the time horizon, given hospitalisation
##' 
##' * `lambda_home`: conditional probability that a patient recovers given they are not hospitalised
##' 
##' * `lambda_hosp`: conditional probability that a patient in hospital recovers given they do not die
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
    ce <- calculate_costs_effects(p_side_effects_t1, 
                                  p_side_effects_t2,
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
    ce <- calculate_costs_effects(p_side_effects_t1, 
                                  p_side_effects_t2,
                                  p_home_home, p_home_hospital, p_home_recover,
                                  p_hospital_hospital, p_hospital_recover, p_hospital_dead,
                                  c_home_care, c_hospital, c_death,
                                  u_recovery, u_home_care, u_hospital)
    t(ce[1,,])
}
