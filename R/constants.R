# One-dimensional integrals of the closed forms: A_c (first-order constant of the calibration slope), nu_0(c)
# (design effect of the score weights at kappa = 0) and tau_c (fixed-p standard deviation constant of the slope).

problem_spec <- function(method, prevalence) {
  switch(method,
         unweighted = list(pi = prevalence, c = 1, kmult = 1),
         weights = list(pi = prevalence, c = (1 - prevalence) / prevalence, kmult = 1),
         undersampling = list(pi = 0.5, c = 1, kmult = 1 / (2 * prevalence)),
         stop("method must be 'unweighted', 'weights' or 'undersampling'"))
}

const_A <- function(s, pi, c) {
  nd <- panel_nodes(s); x <- nd$x
  phi <- stats::dnorm(x + s) * nd$w
  f <- 2 * s * x + logit(pi) + log(c)
  p <- expit(f); q1 <- expit(-f)
  Wq <- phi * p; Z <- sum(Wq)
  m <- function(h) sum(Wq * h) / Z
  mf <- m(f); varf <- m((f - mf)^2)
  qq <- p * q1; covqf <- m((qq - m(qq)) * (f - mf))
  Epc <- m(p + c * q1)
  (Epc - (c - 1) * covqf) / ((1 - pi) * Z * varf)
}

const_nu0 <- function(s, pi, c) {
  nd <- panel_nodes(s); x <- nd$x
  Ep <- 0; Ep2 <- 0
  for (cls in 1:2) {
    mu <- if (cls == 1) s else -s; py <- if (cls == 1) pi else 1 - pi
    wx <- stats::dnorm(x - mu) * nd$w * py
    eta <- expit(2 * s * x + logit(pi)); f <- 2 * s * x + logit(pi) + log(c)
    Ep <- Ep + sum(wx * (eta * c * expit(-f) + (1 - eta) * expit(f)))
    Ep2 <- Ep2 + sum(wx * (eta * c^2 * expit(-f)^2 + (1 - eta) * expit(f)^2))
  }
  Ep2 / Ep^2
}

const_tau0 <- function(s, pi, c) {
  nd <- panel_nodes(s); x <- nd$x
  A <- matrix(0, 2, 2); B <- matrix(0, 2, 2)
  for (cls in 1:2) {
    mu <- if (cls == 1) s else -s; py <- if (cls == 1) pi else 1 - pi; yv <- if (cls == 1) 1 else 0
    wx <- stats::dnorm(x - mu) * nd$w * py
    f <- 2 * s * x + logit(pi) + log(c); p <- expit(f); cy <- if (yv == 1) c else 1
    X <- rbind(1, x)
    A <- A + (X * rep(wx * cy * p * (1 - p), each = 2)) %*% t(X)
    B <- B + (X * rep(wx * cy^2 * (yv - p)^2, each = 2)) %*% t(X)
  }
  Ai <- solve(A); V <- Ai %*% B %*% Ai
  sqrt(V[2, 2]) / (2 * s)
}

#' Constants of the closed-form criteria
#'
#' The first-order constant \eqn{A} of the calibration slope, the design effect \eqn{\nu_0} of the score weights and
#' the fixed-\eqn{p} standard deviation constant \eqn{\tau} of the slope, for one fitting method. For undersampling,
#' \eqn{A} and \eqn{\nu_0} are those of the balanced problem and \eqn{\tau} is on the scale of the original
#' sample size (it includes the factor \eqn{(2\pi)^{-1/2}}).
#'
#' @inheritParams kappa_sep
#' @param method One of \code{"unweighted"}, \code{"weights"} or \code{"undersampling"}.
#' @return A list with \code{A}, \code{nu0}, \code{tau} and \code{s}.
#' @examples
#' criteria_constants(0.81, 0.0916, "unweighted")
#' @export
criteria_constants <- function(cstatistic, prevalence, method = "unweighted") {
  s <- stats::qnorm(cstatistic) / sqrt(2); sp <- problem_spec(method, prevalence)
  tau <- const_tau0(s, sp$pi, sp$c)
  if (method == "undersampling") tau <- tau / sqrt(2 * prevalence)
  list(A = const_A(s, sp$pi, sp$c), nu0 = const_nu0(s, sp$pi, sp$c), tau = tau, s = s)
}
