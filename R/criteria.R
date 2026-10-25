# Criteria C1 and C2, exact and closed form, and their probability versions.

check_inputs <- function(cstatistic, prevalence, parameters) {
  if (!is.numeric(cstatistic) || cstatistic <= 0.5 || cstatistic >= 1) stop("cstatistic must lie in (0.5, 1)")
  if (!is.numeric(prevalence) || prevalence <= 0 || prevalence >= 0.5)
    stop("prevalence must lie in (0, 0.5); for a common outcome, code the rarer class as the event")
  if (!is.numeric(parameters) || parameters < 2) stop("parameters must be a number of predictor columns of at least 2")
  invisible(TRUE)
}

# Limiting slope and direction quality of method m at sample size n.
limits_at <- function(s, prevalence, parameters, n, method) {
  sp <- problem_spec(method, prevalence)
  k <- parameters / n * sp$kmult
  fp <- fixed_point(s, sp$pi, k, sp$c)
  list(cs = 2 * fp$m / fp$q, rho = fp$alpha / sqrt(fp$q), kappa_fit = k, res = fp$res)
}

# Smallest n with g(n) >= 0, g increasing in n; search on log n below the separability threshold.
solve_n <- function(g, n_lo, n_hi = 1e7) {
  if (g(n_hi) < 0) return(NA_real_)
  r <- stats::uniroot(function(ln) g(exp(ln)), c(log(n_lo), log(n_hi)), tol = 1e-6)
  exp(r$root)
}

# Exact criterion as one system: the fixed-point equations plus "criterion = target", with log kappa as a fifth
# unknown. crit(a, sigma, kappa) returns the criterion value minus its target. Starts from the fixed point at the
# closed-form kappa; falls back to a root search in n if the joint Newton fails.
solve_exact <- function(s, sp, parameters, crit, kappa0, fallback) {
  fp <- fixed_point(s, sp$pi, kappa0, sp$c)
  x0 <- c(fp$alpha, log(fp$sigma), NA, fp$b, log(kappa0))
  base <- paste(signif(c(s, sp$pi, sp$c), 12), collapse = "_")
  pp <- .imbsize_env$path[[base]]
  x0[3] <- pp$x[[which.min(abs(pp$k - kappa0))]][3]
  fn <- function(par) {
    k <- exp(par[5])
    c(fp_residuals(par[1:4], s, sp$pi, k, sp$c), crit(par[1], exp(par[2]), k))
  }
  sol <- newton_solve(x0, fn, tol = 1e-11)
  if (is.finite(sol$res) && sol$res < 1e-9) return(parameters * sp$kmult / exp(sol$x[5]))
  fallback()
}

min_n <- function(prevalence, parameters, method, ksep) {
  sp <- problem_spec(method, prevalence)
  parameters * sp$kmult / (0.95 * ksep) # n at which the fitted problem reaches 95% of the threshold
}

#' Criterion C1: calibration slope
#'
#' Development sample size giving an expected calibration slope of at least \code{slope} for a model used without
#' recalibration, fitted without rebalancing, with class weights or with random undersampling.
#'
#' @inheritParams kappa_sep
#' @param parameters Number of predictor columns after coding.
#' @param method One or more of \code{"unweighted"}, \code{"weights"} and \code{"undersampling"}.
#' @param slope Target expected calibration slope \eqn{S}.
#' @param assurance If not \code{NULL}, the probability with which the slope must reach \code{slope}
#'   (probability version: proportional-regime mean with standard deviation \eqn{cs(\kappa)\tau/\sqrt{n}}).
#' @param form \code{"exact"} (fixed point, recommended) or \code{"closed"} (hyperbolic closed form, which
#'   slightly underestimates the exact size).
#' @return A data frame with one row per method: total sample size \code{n}, events, EPV and the multiplier relative
#'   to the unweighted fit.
#' @examples
#' ss_c1(0.81, 0.0916, 71, method = "unweighted", form = "closed")
#' @export
ss_c1 <- function(cstatistic, prevalence, parameters, method = c("unweighted", "weights", "undersampling"),
                  slope = 0.9, assurance = NULL, form = c("exact", "closed")) {
  check_inputs(cstatistic, prevalence, parameters); form <- match.arg(form)
  s <- stats::qnorm(cstatistic) / sqrt(2)
  ms <- union("unweighted", method)
  d <- parameters - 1 # kappa = (p - 1) / n: the leading terms are driven by the p - 1 orthogonal coordinates
  n <- sapply(ms, function(m) {
    cst <- criteria_constants(cstatistic, prevalence, m); sp <- problem_spec(m, prevalence)
    z <- if (is.null(assurance)) 0 else stats::qnorm(assurance)
    # probability version: sd of the slope = cs(kappa) * tau / sqrt(n)
    if (form == "closed") {
      Ak <- cst$A * d * sp$kmult
      if (z == 0) return(Ak * slope / (1 - slope))
      return(solve_n(function(n) (1 - z * cst$tau / sqrt(n)) / (1 + Ak / n) - slope, 10))
    }
    ks <- kappa_sep(cstatistic, if (m == "undersampling") 0.5 else prevalence)
    nlo <- min_n(prevalence, d, m, ks)
    g <- function(n) limits_at(s, prevalence, d, n, m)$cs * (1 - z * cst$tau / sqrt(n)) - slope
    Ak <- cst$A * d * sp$kmult
    n0 <- if (z == 0) Ak * slope / (1 - slope) else solve_n(function(n) (1 - z * cst$tau / sqrt(n)) / (1 + Ak / n) - slope, 10)
    k0 <- min(d * sp$kmult / n0, 0.9 * ks)
    crit <- function(a, sg, k) 2 * a * s / (a^2 + sg^2) * (1 - z * cst$tau / sqrt(d * sp$kmult / k)) - slope
    solve_exact(s, sp, d, crit, k0, function() solve_n(g, nlo))
  })
  out_table(n[method], n["unweighted"], prevalence, parameters, "C1")
}

#' Criterion C2: loss of discrimination
#'
#' Development sample size for which the expected AUC of a model that will be recalibrated (by Platt scaling on
#' leave-one-out scores) is within \code{auc_loss} of the anticipated AUC.
#'
#' @inheritParams ss_c1
#' @param auc_loss Tolerated loss of AUC \eqn{\delta}.
#' @return As for \code{\link{ss_c1}}.
#' @examples
#' ss_c2(0.81, 0.0916, 71, method = "unweighted", form = "closed")
#' @export
ss_c2 <- function(cstatistic, prevalence, parameters, method = c("unweighted", "weights", "undersampling"),
                  auc_loss = 0.01, assurance = NULL, form = c("exact", "closed")) {
  check_inputs(cstatistic, prevalence, parameters); form <- match.arg(form)
  s <- stats::qnorm(cstatistic) / sqrt(2)
  delta <- if (is.null(assurance)) auc_loss else auc_loss * (parameters - 1) / stats::qchisq(assurance, parameters - 1)
  target <- cstatistic - delta
  ms <- union("unweighted", method)
  d <- parameters - 1
  n <- sapply(ms, function(m) {
    sp <- problem_spec(m, prevalence)
    if (form == "closed") {
      cst <- criteria_constants(cstatistic, prevalence, m)
      rs <- stats::qnorm(target) / (sqrt(2) * s)
      return(d * sp$kmult * cst$nu0 / (s^2 * (rs^-2 - 1)))
    }
    ks <- kappa_sep(cstatistic, if (m == "undersampling") 0.5 else prevalence)
    nlo <- min_n(prevalence, d, m, ks)
    g <- function(n) stats::pnorm(sqrt(2) * s * limits_at(s, prevalence, d, n, m)$rho) - target
    cst <- criteria_constants(cstatistic, prevalence, m)
    rs <- stats::qnorm(target) / (sqrt(2) * s)
    n0 <- d * sp$kmult * cst$nu0 / (s^2 * (rs^-2 - 1))
    k0 <- min(d * sp$kmult / n0, 0.9 * ks)
    crit <- function(a, sg, k) stats::pnorm(sqrt(2) * s * a / sqrt(a^2 + sg^2)) - target
    solve_exact(s, sp, d, crit, k0, function() solve_n(g, nlo))
  })
  out_table(n[method], n["unweighted"], prevalence, parameters, "C2")
}

out_table <- function(n, n_unw, prevalence, parameters, crit) {
  data.frame(criterion = crit, method = names(n), n = ceiling(n), events = ceiling(n * prevalence),
             EPV = round(n * prevalence / parameters, 1), multiplier = round(n / n_unw, 2), row.names = NULL)
}
