test_that("EVPI",{
    expect_equal(evpi(chemo_nb), 42.178, tolerance=0.01)
    expect_equal(evpi(chemo_cea)$evpi, c(1.02877281101973, 12.3534189284037, 42.1780795210652, 
                                      58.9895661745759, 40.7463368296667), tolerance=0.01)
})
