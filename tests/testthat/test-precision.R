test_that("precision requirement matches the paper and tau^2 = A for the unweighted fit", {
  expect_equal(ss_recal_precision(0.81, 0.0916), 2986)
  k <- criteria_constants(0.81, 0.0916, "unweighted")
  expect_equal(k$tau^2, k$A, tolerance = 1e-8)
  kw <- criteria_constants(0.75, 0.2, "weights")
  expect_equal(kw$tau^2, kw$A, tolerance = 1e-8)
  expect_gt(ss_recal_precision(0.81, 0.0916, epsilon = 0.05), ss_recal_precision(0.81, 0.0916))
})
