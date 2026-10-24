test_that("closed-form constants match the reference implementation", {
  for (r in ref) {
    k <- criteria_constants(r$c, r$pi, r$method)
    expect_equal(k$A, r$A, tolerance = 1e-6)
    expect_equal(k$nu0, r$nu0, tolerance = 1e-6)
    expect_equal(k$tau, r$tau, tolerance = 1e-6)
  }
})

test_that("separability threshold matches", {
  for (r in ref[sapply(ref, function(x) x$method == "unweighted")]) expect_equal(kappa_sep(r$c, r$pi), r$ksep, tolerance = 1e-5)
})

test_that("closed forms match", {
  for (w in c("IST", "NoShow")) {
    rr <- ref[sapply(ref, function(x) x$world == w)]
    a <- ss_c1(rr[[1]]$c, rr[[1]]$pi, rr[[1]]$p, form = "closed")
    b <- ss_c2(rr[[1]]$c, rr[[1]]$pi, rr[[1]]$p, form = "closed")
    for (r in rr) {
      expect_equal(a$n[a$method == r$method], ceiling(r$c1cf))
      expect_equal(b$n[b$method == r$method], ceiling(r$c2cf))
    }
  }
})

test_that("Gaussian quadrature integrates correctly", {
  Q <- normal_nodes()
  expect_equal(sum(Q$w), 1)
  expect_equal(sum(Q$w * Q$x^2), 1, tolerance = 1e-10)
  expect_equal(sum(Q$w * Q$x^4), 3, tolerance = 1e-10)
})
