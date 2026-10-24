# Gauss-Legendre nodes by the Golub-Welsch algorithm, and expectations over N(0, 1).

gauss_legendre <- function(n) {
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
}

.imbsize_env <- new.env(parent = emptyenv())

# Nodes and normalized weights for E[f(Z)], Z ~ N(0, 1): Gauss-Hermite rule (probabilists' Hermite polynomials,
# Golub-Welsch), with getOption("imbsize.nodes", 120) nodes.
normal_nodes <- function() {
  N <- getOption("imbsize.nodes", 120)
  if (is.null(.imbsize_env$normal) || length(.imbsize_env$normal$x) != N) {
    i <- seq_len(N - 1)
    J <- matrix(0, N, N); J[cbind(i, i + 1)] <- sqrt(i); J[cbind(i + 1, i)] <- sqrt(i)
    e <- eigen(J, symmetric = TRUE); o <- order(e$values)
    w <- e$vectors[1, o]^2
    .imbsize_env$normal <- list(x = e$values[o], w = w / sum(w))
    .imbsize_env$fp <- NULL; .imbsize_env$path <- NULL
  }
  .imbsize_env$normal
}

# Composite Gauss-Legendre nodes on [-s - 12, s + 12] (20 nodes per panel), for the one-dimensional integrals.
panel_nodes <- function(s) {
  if (is.null(.imbsize_env$gl20)) .imbsize_env$gl20 <- gauss_legendre(20)
  g <- .imbsize_env$gl20
  lo <- -s - 12; hi <- s + 12
  h <- min(0.5, 1 / (4 * s)); m <- ceiling((hi - lo) / h)
  e <- seq(lo, hi, length.out = m + 1)
  mid <- (e[-1] + e[-(m + 1)]) / 2; half <- (e[-1] - e[-(m + 1)]) / 2
  list(x = as.vector(outer(g$x, half) + rep(mid, each = 20)),
       w = as.vector(outer(g$w, half)))
}
