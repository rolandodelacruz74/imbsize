# Expected validation performance at a given development size.

logloss_mix <- function(slope, offset, m, q, b, pi) {
  Q <- normal_nodes(); ll <- 0
  for (cls in 1:2) {
    y <- if (cls == 1) 1 else -1; py <- if (cls == 1) pi else 1 - pi
    z <- slope * (y * m + sqrt(q) * Q$x) + offset
    ll <- ll + py * sum(Q$w * (if (y == 1) log1p(exp(-z)) else log1p(exp(z))))
  }
  ll
}

#' Expected performance of a fit at a given development size
#'
#' Limiting validation calibration slope without recalibration, AUC (unchanged by recalibration), calibration slope
#' after Platt scaling on \eqn{K}-fold cross-fitted scores, and log-loss after Platt scaling on leave-one-out
#' scores, for a development sample of size \code{n}.
#'
#' @inheritParams ss_c1
#' @param n Development sample size.
#' @param folds Number of folds \eqn{K} for cross-fitted Platt scaling.
#' @return A data frame with one row per method.
#' @examples
#' \donttest{predict_performance(0.81, 0.0916, 71, n = 3876)}
#' @export
predict_performance <- function(cstatistic, prevalence, parameters, n,
                                method = c("unweighted", "weights", "undersampling"), folds = 5) {
  check_inputs(cstatistic, prevalence, parameters)
  s <- stats::qnorm(cstatistic) / sqrt(2)
  rows <- lapply(method, function(m) {
    sp <- problem_spec(m, prevalence); k <- parameters / n * sp$kmult
    ks <- kappa_sep(cstatistic, sp$pi)
    if (k >= ks) {
      warning(sprintf("%s: kappa = %.4f >= kappa_sep = %.4f; the unpenalized fit does not exist in the limit", m, k, ks))
      return(data.frame(method = m, kappa = k, kappa_sep = ks, slope = NA, auc = NA, slope_plattK = NA, logloss_loo = NA))
    }
    fp <- fixed_point(s, sp$pi, k, sp$c)
    cs <- 2 * fp$m / fp$q
    auc <- stats::pnorm(sqrt(2) * fp$m / sqrt(fp$q))
    kf <- k * folds / (folds - 1)
    pk <- if (kf < ks) { ff <- fixed_point(s, sp$pi, kf, sp$c); cs / (2 * ff$m / ff$q) } else NA
    ll <- logloss_mix(cs, logit(prevalence), fp$m, fp$q, fp$b, prevalence)
    data.frame(method = m, kappa = k, kappa_sep = ks, slope = cs, auc = auc, slope_plattK = pk, logloss_loo = ll)
  })
  out <- do.call(rbind, rows); attr(out, "n") <- n; attr(out, "folds") <- folds
  out
}
