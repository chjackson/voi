test_that("EVPI",{
  expect_equal(evpi(chemo_nb), 1535.81127362873, tolerance=0.1)
  expect_equal(evpi(chemo_cea)$evpi, 
               c(550.146022178757, 1535.81127362873, 968.76194762066, 759.004292830592, 
                 667.890283294953), tolerance=0.01)
})
