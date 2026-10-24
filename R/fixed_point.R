# Proportional-regime limits of the (class-weighted) unpenalized logistic fit under the two-class Gaussian model
# with a common covariance (equations (E1)-(E4) of the companion paper; Proposition A.2 of the paper).

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

# t solving t + lc * expit(-t) = v (proximal point of the logistic loss), elementwise: bisection then Newton.
prox_logistic <- function(v, lc) {
  # g(t) = t - v - lc * expit(-t) is increasing with g' >= 1 and its root lies in [v, v + lc]:
  # safeguarded Newton from one fixed-point step, then bisection for any element that has not converged.
  lo <- v; hi <- v + lc
  t <- v + lc * expit(-v)
  for (k in 1:12) {
    e <- expit(-t)
    t <- t - (t - v - lc * e) / (1 + lc * e * (1 - e))
    t <- pmin(pmax(t, lo), hi)
  }
  bad <- which(!is.finite(t) | abs(t - v - lc * expit(-t)) > 1e-10 * (1 + abs(t)))
  if (length(bad)) {
    l <- lo[bad]; h <- hi[bad]; vv <- v[bad]; lcb <- if (length(lc) > 1) lc[bad] else lc
    for (k in 1:60) {
      mid <- 0.5 * (l + h); f <- mid - vv - lcb * expit(-mid); neg <- f < 0; neg[is.na(neg)] <- FALSE
      l[neg] <- mid[neg]; h[!neg] <- mid[!neg]
    }
    t[bad] <- 0.5 * (l + h)
  }
  t
}

fp_residuals <- function(par, s, pi, kappa, c) {
  Q <- normal_nodes()
  a <- par[1]; sig <- exp(par[2]); lam <- exp(par[3]); b <- par[4]
  tau <- sqrt(a^2 + sig^2)
  e <- numeric(4)
  for (cls in 1:2) {
    y <- if (cls == 1) 1 else -1
    py <- if (cls == 1) pi else 1 - pi
    cy <- if (cls == 1) c else 1
    v <- a * s + y * b + tau * Q$x
    t <- prox_logistic(v, lam * cy)
    psi <- cy * expit(-t)
    l2 <- cy * expit(t) * expit(-t)
    dpsi <- -l2 / (1 + lam * l2)
    e[1] <- e[1] + py * sum(Q$w * y * psi)
    e[2] <- e[2] + py * sum(Q$w * (s * psi + a * dpsi))
    e[3] <- e[3] + py * sum(Q$w * psi^2)
    e[4] <- e[4] + py * sum(Q$w / (1 + lam * l2))
  }
  c(e[1], e[2], lam^2 * e[3] - kappa * sig^2, e[4] - (1 - kappa))
}

# Damped Newton with a forward-difference Jacobian.
newton_solve <- function(x0, fn, tol = 1e-12, maxit = 60) {
  x <- x0; f <- fn(x); nf <- max(abs(f))
  for (it in seq_len(maxit)) {
    if (!all(is.finite(f)) || nf < tol) break
    d <- length(x); J <- matrix(0, d, d)
    for (j in seq_len(d)) {
      h <- 1e-7 * max(1, abs(x[j])); xh <- x; xh[j] <- xh[j] + h
      J[, j] <- (fn(xh) - f) / h
    }
    step <- tryCatch(solve(J, -f), error = function(e) NULL)
    if (is.null(step)) break
    tstep <- 1
    repeat {
      xn <- x + tstep * step; fnew <- fn(xn); nn <- max(abs(fnew))
      if (all(is.finite(fnew)) && nn < nf) break
      tstep <- tstep / 2
      if (tstep < 1e-6) break
    }
    if (tstep < 1e-6) break
    x <- xn; f <- fnew; nf <- nn
  }
  list(x = x, res = nf)
}

#' Proportional-regime limit of the logistic fit
#'
#' Solves the fixed-point equations of the unpenalized, possibly class-weighted, logistic fit under a two-class
#' Gaussian model with a common covariance, in the regime \eqn{p/n \to \kappa}.
#'
#' @param s Separation, \eqn{s = \Phi^{-1}(\mathrm{AUC})/\sqrt2}.
#' @param prevalence Prevalence of the problem the model is fitted on.
#' @param kappa Ratio of predictor columns to patients of the fitted problem.
#' @param c Weight of the events (1 for the unweighted fit).
#' @return A list with \code{alpha}, \code{sigma}, \code{b}, \code{m} (\eqn{= \alpha s}), \code{q}
#'   (\eqn{= \alpha^2 + \sigma^2}) and the maximal residual \code{res}.
#' @keywords internal
fixed_point <- function(s, prevalence, kappa, c = 1) {
  key <- paste(signif(c(s, prevalence, kappa, c), 12), collapse = "_")
  if (!is.null(.imbsize_env$fp[[key]])) return(.imbsize_env$fp[[key]])
  pi <- prevalence
  solve_at <- function(x, k, lc) newton_solve(x, function(par) fp_residuals(par, s, pi, k, exp(lc)))
  advance <- function(x, k0, l0, k1, l1, depth = 0) {
    sol <- solve_at(x, k1, l1)
    if (!all(is.finite(sol$x))) sol <- list(x = x, res = Inf)
    if (sol$res <= 1e-9 || depth >= 4) return(sol)
    km <- 0.5 * (k0 + k1); lm <- 0.5 * (l0 + l1)
    mid <- advance(x, k0, l0, km, lm, depth + 1)
    advance(mid$x, km, lm, k1, l1, depth + 1)
  }
  base <- paste(signif(c(s, prevalence, c), 12), collapse = "_")
  path <- .imbsize_env$path[[base]]
  if (!is.null(path) && length(path$k)) {
    # warm start from the nearest kappa already solved for this (s, prevalence, c)
    j <- which.min(abs(log(path$k) - log(kappa)))
    x <- path$x[[j]]; kp <- path$k[j]
    for (k in exp(seq(log(kp), log(kappa), length.out = 6))[-1]) { sol <- advance(x, kp, log(c), k, log(c)); x <- sol$x; kp <- k }
    if (sol$res > 1e-9) path <- NULL
  }
  if (is.null(path) || !length(path$k) || sol$res > 1e-9) {
    k0 <- min(0.002, 0.5 * kappa)
    x <- c(2 * s, log(sqrt(20 * k0)), log(k0 / (pi * (1 - pi))), logit(pi))
    sol <- solve_at(x, k0, 0); x <- sol$x
    if (kappa > k0) {
      ks <- exp(seq(log(k0), log(kappa), length.out = 15)); kp <- k0
      for (k in ks[-1]) { sol <- advance(x, kp, 0, k, 0); x <- sol$x; kp <- k }
    }
  if (c != 1) {
    ls <- seq(0, log(c), length.out = max(3, ceiling(abs(log(c)) * 5)))
    for (j in seq_along(ls)[-1]) { sol <- advance(x, kappa, ls[j - 1], kappa, ls[j]); x <- sol$x }
  }
  }
  if (is.null(.imbsize_env$path)) .imbsize_env$path <- list()
  pp <- .imbsize_env$path[[base]]; if (is.null(pp)) pp <- list(k = numeric(0), x = list())
  if (sol$res <= 1e-9) { pp$k <- c(pp$k, kappa); pp$x[[length(pp$x) + 1]] <- x; .imbsize_env$path[[base]] <- pp }
  out <- list(alpha = x[1], sigma = exp(x[2]), b = x[4], m = x[1] * s, q = x[1]^2 + exp(2 * x[2]), res = sol$res)
  if (is.null(.imbsize_env$fp)) .imbsize_env$fp <- list()
  .imbsize_env$fp[[key]] <- out
  out
}

# E[(c - Z)_+^2] for Z ~ N(0, 1).
epos2 <- function(c) (1 + c^2) * stats::pnorm(c) + c * stats::dnorm(c)

#' Separability threshold
#'
#' Ratio \eqn{\kappa = p/n} above which the unpenalized logistic fit does not exist in the proportional limit
#' (the reciprocal of the threshold of Mignacco et al. and Kammoun and Alouini).
#'
#' @param cstatistic Anticipated AUC (C-statistic).
#' @param prevalence Prevalence.
#' @return The threshold \eqn{\kappa_{\mathrm{sep}}}.
#' @examples
#' kappa_sep(0.81, 0.0916)
#' @export
kappa_sep <- function(cstatistic, prevalence) {
  s <- stats::qnorm(cstatistic) / sqrt(2); pi <- prevalence
  f <- function(x) {
    r <- tanh(x[1])
    (pi * epos2(-r * s - x[2]) + (1 - pi) * epos2(-r * s + x[2])) / (1 - r^2)
  }
  best <- Inf
  for (a in c(0, 0.5, 1, 2)) for (b in c(-3, -1, 0, 1)) {
    o <- stats::optim(c(a, b), f, method = "Nelder-Mead", control = list(reltol = 1e-14, maxit = 5000))
    best <- min(best, o$value)
  }
  best
}
