if (require("heemod",quietly=TRUE)){

  set.seed(1)
  
  ### START OF CODE FROM help(run_psa) in heemod   

  mod1 <- define_strategy(
    transition = define_transition(
      .5, .5,
      .1, .9
    ),
    define_state(
      cost = cost_init + age * 5,
      ly = 1
    ),
    define_state(
      cost = cost_init + age,
      ly = 0
    )
  )
  
  mod2 <- define_strategy(
    transition = define_transition(
      p_trans, C,
      .1, .9
    ),
    define_state(
      cost = 789 * age / 10,
      ly = 1
    ),
    define_state(
      cost = 456 * age / 10,
      ly = 0
    )
    
  )
  
  res2 <- run_model(
    mod1, mod2,
    parameters = define_parameters(
      age_init = 60,
      cost_init = 1000,
      age = age_init + markov_cycle,
      p_trans = .7
    ),
    init = 1:0,
    cycles = 10,
    cost = cost,
    effect = ly
  )
  
  rsp <- define_psa(
    age_init ~ normal(60, 10),
    cost_init ~ normal(1000, 100),
    p_trans ~ binomial(.7, 100),
    correlation = matrix(c(
      1,  .4, 0,
      .4, 1,  0,
      0,  0,  1
    ), byrow = TRUE, ncol = 3)
  )
  
  
  # with run_model result
  # (only 10 resample for speed)
  ndt1 <- run_psa(res2, psa = rsp, N = 10)

  ### END OF CODE FROM help(run_psa) in heemod   
  
  
  if (requireNamespace("BCEA",quietly=TRUE)){
    test_that("heemod interface", {
      outputs <- import_heemod_outputs(ndt1, k=seq(10000,50000,by=10000))
      inputs <- import_heemod_inputs(ndt1)
      
      expect_equal(evpi(outputs)$evpi[2], 4883, tol=1) 
      expect_equal(evppi(outputs, inputs, pars="p_trans")$evppi[2], 4895, tol=1)
    })
  } 
}
