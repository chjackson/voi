
test_that("enbs",{
  nbs <- enbs(chemo_evsi_or, costs_setup = c(5e6, 1e7), costs_pp = c(28000,42000), 
              pop = 46000, time = 10)
  expect_equal(nbs$enbs[1], -9250000)
  nbs <- enbs(chemo_evsi_or, costs_setup = c(7500000), costs_pp = c(28000,42000), 
              pop = 46000, time = 10)
  expect_equal(nbs$enbs[1], -9250000)
  cm <- 7500000
  nbs <- enbs(chemo_evsi_or, costs_setup = c(cm/2, cm, cm*2), costs_pp = c(28000,42000), 
              pop = 46000, time = 10)
  expect_equal(nbs$enbs[1], -9250000)
})

test_that("enbs maxima are computed correctly",{
  evsi_permuted <- chemo_evsi_or[sample(1:nrow(chemo_evsi_or), replace=FALSE),]
  cm <- 7500000
  nbs <- enbs(evsi_permuted, costs_setup = c(cm/2, cm, cm*2), 
              costs_pp = c(28000,42000), 
              pop = 46000, time = 10)
  nbs20 <- nbs[nbs$k==20000,] 
  nmax <- nbs20$enbsmax[1]
  expect_equal(max(nbs20$enbs), nmax)
  expect_equal(nbs20$n[which.max(nbs20$enbs)], nbs20$nmax[1])
})
  

test_that("enbs_opt",{
  nbs <- enbs(chemo_evsi_or, costs_setup = c(5e6, 1e7),
              costs_pp = c(28000,42000), 
              pop = 46000, time = 10) %>%
    filter(k==20000)
  nbopt <- enbs_opt(nbs)
  expect_equal(nbopt$nmax, 450)
  nbopts <- enbs_opt(nbs, smooth=TRUE)
  expect_equal(nbopts$nmax, 427)
  nbopts <- enbs_opt(nbs, smooth=TRUE, smooth_df=4)
  expect_equal(nbopts$nmax, 578)
  nbopt2 <- enbs_opt(nbs, pcut=0.01)
  expect_lt(nbopt2$nupper, nbopt$nupper)
})

test_that("pop_voi",{
  evs <- chemo_evsi_or %>% filter(k==20000) %>% pull(evsi) %>% head(3)
  pv <- pop_voi(evs, pop=100, time=10, dis=0.035)
  pv0 <- pop_voi(evs, pop=100, time=10, dis=0)
  pvm <- pop_voi(evs, pop=100, time=10, dis=c(0.035,0,0))
  expect_equal(pvm[1], pv[1])  
  expect_equal(pvm[2], pv0[2])  
  expect_equal(pvm[3], pv0[3])  
})
