test_that("moment matching method",{
  evm <- evsi(outputs=chemo_nb, inputs=chemo_pars, 
              pars="p_side_effects_t1",
              method = "mm",
              study =  "binary",
              n = c(100), Q = 5, 
              analysis_args = list(a=53, b=60, n=100),
              model_fn = chemo_model_nb,  par_fn = chemo_pars_fn)
  expect_true(is.data.frame(evm) && nrow(evm)==1)
  
  evm <- evsi(outputs=chemo_cea, inputs=chemo_pars, 
              pars="p_side_effects_t1",
              method = "mm",
              study =  "binary",
              n = c(100), Q = 5, 
              analysis_args = list(a=53, b=60, n=100),
              model_fn = chemo_model_cea,  par_fn = chemo_pars_fn)
  expect_true(is.data.frame(evm) && nrow(evm)==5)
})


test_that("moment matching method with multiple sample sizes",{
  ev <-  evsi(outputs=chemo_nb, inputs=chemo_pars, pars="p_side_effects_t2", 
              method="mm", study="binary", n=seq(10, 110, by=10), Q=10, 
              analysis_args = list(a=53, b=60),
              model_fn = chemo_model_nb, par_fn = chemo_pars_fn)
  expect_true(is.data.frame(ev))
  evc <- evsi(outputs=chemo_cea, inputs=chemo_pars, pars="p_side_effects_t2", 
              method="mm", study="binary", n=seq(10, 100, by=10), Q=10, 
              analysis_args = list(a=53, b=60),
              model_fn = chemo_model_cea, par_fn = chemo_pars_fn)
  expect_true(is.data.frame(evc))
})

test_that("errors in moment matching method",{
  expect_error(
    evsi(outputs=chemo_cea, inputs=chemo_pars, pars="p_side_effects_t1",
         method = "mm", study =  "binary", n = c(100), Q = 5, 
         analysis_fn = analysis_binary, 
         analysis_args = list(a=53, b=60, n=100),
         model_fn = chemo_model_nb,  par_fn = chemo_pars_fn),
    "output of model_fn should have two rows"
  )

  datagen_fn <- function(inputs, n=100){
    nsim <- nrow(inputs)
    data.frame(x1 = rbinom(nsim, size=n, prob=inputs[,"p_side_effects_t1"]))
  }
  
  expect_error(  
    evsi(outputs=chemo_nb, inputs=chemo_pars, 
         pars="p_side_effects_t1",
         method = "mm",
         datagen_fn = datagen_fn, 
         model_fn = chemo_model_nb,  par_fn = chemo_pars_fn),
    "`analysis_fn` should be supplied")
  
  expect_error(  
    evsi(outputs=chemo_nb, inputs=chemo_pars, 
         pars="p_side_effects_t1",
         method = "mm",
         datagen_fn = datagen_fn, 
         model_fn = chemo_model_nb,  par_fn = chemo_pars_fn),
    "`analysis_fn` should be supplied")

})


test_that("user-written analysis function",{

  if (requireNamespace("rjags", quietly = TRUE)) {
    
    datagen_fn <- function(inputs, n=100){
      p1 <- inputs[,"p_side_effects_t1"]
      logor <- inputs[,"logor_side_effects"]
      odds1 <- p1 / (1 - p1)
      odds2 <- odds1 * exp(logor)
      p2 <- odds2 / (1 + odds2)
      nsim <- nrow(inputs)
      data.frame(y1 = rbinom(nsim, n, p1),
                 y2 = rbinom(nsim, n, p2))
    }
    
    ## Constants hard coded in analysis_fn 
    analysis_fn <- function(data, args, pars){
      dat <- list(y=c(data[,"y1"], data[,"y2"]))
      design <- list(n = rep(args$n, 2))
      priors <- list(a1=1, b1=1, mu=0, sigma=2)
      jagsdat <- c(dat, design, priors)
      or_jagsmod <- "
      model {
        y[1] ~ dbinom(p[1], n[1])
        y[2] ~ dbinom(p[2], n[2])
        p[1] <- p1
        p[2] <- odds[2] / (1 + odds[2])
        p1 ~ dbeta(a1, b1)
        odds[1] <- p[1] / (1 - p[1])
        odds[2] <- odds[1] * exp(logor)
        logor ~ dnorm(mu, 1/sigma^2)
      }
      "
      or.jag <- rjags::jags.model(textConnection(or_jagsmod), 
                                  data=jagsdat, inits=list(logor=0, p1=0.5), quiet = TRUE)
      update(or.jag, 100, progress.bar="none")
      sam <- rjags::coda.samples(or.jag, c("logor"), 500, progress.bar="none")
      data.frame(logor_side_effects = as.numeric(sam[[1]][,"logor"]))
    }

    ev <- evsi(chemo_nb, chemo_pars, 
               pars="logor_side_effects", 
               pars_datagen = c("p_side_effects_t1", "logor_side_effects"), 
               method="mm",
               datagen_fn = datagen_fn, analysis_fn = analysis_fn, 
               n = 100, Q = 10, 
               model_fn = chemo_model_lor_nb, par_fn = chemo_pars_fn)
    expect_true(is.numeric(ev$evsi))
    
    ## Constants specified through analysis_args
    analysis_fn_dynamic <- function(data, args, pars){
      dat <- list(y=c(data[,"y1"], data[,"y2"]))
      design <- list(n = rep(args$n, 2))
      priors <- list(a1=args$a1, b1=args$b1, mu=0, sigma=2)
      jagsdat <- c(dat, design, priors)
      or_jagsmod <- "
      model {
        y[1] ~ dbinom(p[1], n[1])
        y[2] ~ dbinom(p[2], n[2])
        p[1] <- p1
        p[2] <- odds[2] / (1 + odds[2])
        p1 ~ dbeta(a1, b1)
        odds[1] <- p[1] / (1 - p[1])
        odds[2] <- odds[1] * exp(logor)
        logor ~ dnorm(mu, 1/sigma^2)
      }
      "
      or.jag <- rjags::jags.model(textConnection(or_jagsmod), 
                                  data=jagsdat, inits=list(logor=0, p1=0.5), quiet = TRUE)
      update(or.jag, 100, progress.bar="none")
      sam <- rjags::coda.samples(or.jag, c("logor"), 500, progress.bar="none")
      res <- data.frame(logor_side_effects = as.numeric(sam[[1]][,"logor"]))
      names(res) <- pars
      res
    }
    
    ev <- evsi(chemo_nb, chemo_pars, pars="logor_side_effects", 
               pars_datagen=c("p_side_effects_t1","logor_side_effects"), 
               method="mm",
               datagen_fn = datagen_fn, analysis_fn = analysis_fn_dynamic, 
               n = c(100,300), Q = 10, analysis_args = list(a1=2, b1=2),
               model_fn = chemo_model_lor_nb, par_fn = chemo_pars_fn)
    expect_lt(ev$evsi[1], ev$evsi[2])
    
  }
  else {
    message("Not running test of user-written analysis function, since rjags is not installed")
  }
})

test_that("errors in analysis function",{

  datagen_fn <- function(inputs, n=100){
    p1 <- inputs[,"p_side_effects_t1"]
    logor <- inputs[,"logor_side_effects"]
    odds1 <- p1 / (1 - p1)
    odds2 <- odds1 * exp(logor)
    p2 <- odds2 / (1 + odds2)
    data.frame(y1 = rbinom(10000, n, p1),
               y2 = rbinom(10000, n, p2))
  }

  expect_error(
    evsi(chemo_nb, chemo_pars, pars="p_side_effects_t1", method="mm",
         datagen_fn = datagen_fn, analysis_fn = function(){}, 
         n = 10, 
         model_fn = chemo_model_nb, par_fn = chemo_pars_fn),
    "`analysis_fn` should have arguments `data`,`args`,`pars` in that order")

  expect_error(
    evsi(chemo_nb, chemo_pars, pars="logor_side_effects", method="mm",
         datagen_fn = datagen_fn, analysis_fn = function(data,args,pars){}, 
         n = 10, 
         model_fn = chemo_model_nb, par_fn = chemo_pars_fn),
    "Parameters logor_side_effects are not arguments to model_fn")

  analysis_fn_wrongpars <- function(data,args,pars){
    data.frame(wrongpars = rep(0, 100))
  }

  expect_error(
    evsi(chemo_nb, chemo_pars, pars="logor_side_effects", method="mm",
         datagen_fn = datagen_fn, analysis_fn = analysis_fn_wrongpars, 
         n = 10, 
         model_fn = chemo_model_lor_nb, par_fn = chemo_pars_fn),
    "Parameters logor_side_effects not found in data frame returned by `analysis_fn`")
})
