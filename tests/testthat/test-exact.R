test_that("exact criteria match the reference implementation (trial)", {
  skip_on_cran()
  rr <- ref[sapply(ref, function(x) x$world == "IST")]
  a <- ss_c1(0.81, 0.0916, 71); b <- ss_c2(0.81, 0.0916, 71)
  for (r in rr) {
    expect_equal(a$n[a$method == r$method], ceiling(r$c1))
    expect_equal(b$n[b$method == r$method], ceiling(r$c2))
  }
  # agreement with samplesizedev (GitHub commit 2a62e37), which returned 7318 by simulation
  expect_lt(abs(a$n[a$method == "unweighted"] - 7318), 0.01 * 7318)
})

test_that("probability versions match (trial, unweighted)", {
  skip_on_cran()
  r <- ref[[which(sapply(ref, function(x) x$world == "IST" && x$method == "unweighted"))]]
  expect_equal(ss_c1(0.81, 0.0916, 71, "unweighted", assurance = 0.9)$n, ceiling(r$c1p))
  expect_equal(ss_c2(0.81, 0.0916, 71, "unweighted", assurance = 0.9)$n, ceiling(r$c2p))
  expect_equal(ss_c1(0.81, 0.0916, 71, "unweighted", assurance = 0.9, form = "closed")$n, ceiling(r$c1pcf))
})

test_that("predicted slope at the C1 size is the target", {
  skip_on_cran()
  n <- ss_c1(0.81, 0.0916, 71, "unweighted")$n
  pp <- predict_performance(0.81, 0.0916, 71, n = n, method = "unweighted")
  expect_gte(pp$slope, 0.9); expect_lt(pp$slope - 0.9, 1e-3)
  expect_gt(pp$slope_plattK, 1)
})

test_that("exact C1 matches the reference in 36 designs (prevalence x AUC x method)", {
  skip_on_cran()
  p <- 100
  for (i in seq_len(nrow(ref_epv))) {
    r <- ref_epv[i, ]
    a <- ss_c1(r$auc, r$pi, p, method = r$method)
    n_ref <- r$epv * p / r$pi
    expect_lt(abs(a$n - ceiling(n_ref)), 2, label = sprintf("pi %.3f, AUC %.2f, %s", r$pi, r$auc, r$method))
  }
})
