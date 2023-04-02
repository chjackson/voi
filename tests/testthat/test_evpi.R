test_that("EVPI",{
  expect_equal(evpi(chemo_nb), 368.6051, tolerance=0.1)
  expect_equal(evpi(chemo_cea)$evpi, 
               c(8.76486716445652, 368.605096225918, 206.468704777362, 150.472904023947, 
                   126.983450042782), tolerance=0.01)
})
