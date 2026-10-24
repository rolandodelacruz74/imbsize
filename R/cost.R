#' Cost of rebalancing
#'
#' Factor by which class weights and random undersampling multiply the development sample required by C1 or C2,
#' relative to the unweighted fit, over a grid of prevalences and anticipated AUCs. The ratios do not depend on the
#' number of predictor columns.
#'
#' @param cstatistic Vector of anticipated AUCs.
#' @param prevalence Vector of prevalences.
#' @param criterion \code{"C1"} or \code{"C2"}.
#' @param slope,auc_loss Targets of C1 and C2.
#' @param form \code{"closed"} (fast) or \code{"exact"}.
#' @return A data frame of class \code{imbsize_cost}, with a \code{plot} method.
#' @examples
#' rc <- rebalancing_cost(c(0.7, 0.8, 0.9), c(0.05, 0.1, 0.2))
#' plot(rc)
#' @export
rebalancing_cost <- function(cstatistic, prevalence, criterion = c("C1", "C2"), slope = 0.9, auc_loss = 0.01,
                             form = c("closed", "exact")) {
  criterion <- match.arg(criterion); form <- match.arg(form)
  g <- expand.grid(cstatistic = cstatistic, prevalence = prevalence)
  res <- lapply(seq_len(nrow(g)), function(i) {
    f <- if (criterion == "C1") ss_c1 else ss_c2
    a <- if (criterion == "C1") f(g$cstatistic[i], g$prevalence[i], 100, slope = slope, form = form)
         else f(g$cstatistic[i], g$prevalence[i], 100, auc_loss = auc_loss, form = form)
    data.frame(cstatistic = g$cstatistic[i], prevalence = g$prevalence[i],
               weights = a$multiplier[a$method == "weights"], undersampling = a$multiplier[a$method == "undersampling"],
               epv_unweighted = a$EPV[a$method == "unweighted"])
  })
  out <- do.call(rbind, res); attr(out, "criterion") <- criterion
  class(out) <- c("imbsize_cost", "data.frame"); out
}

#' @rdname rebalancing_cost
#' @param x An \code{imbsize_cost} object.
#' @param ... Further graphical arguments passed to \code{plot}.
#' @export
plot.imbsize_cost <- function(x, ...) {
  aucs <- sort(unique(x$cstatistic)); cols <- grDevices::colorRampPalette(c("#6baed6", "#08306b"))(length(aucs))
  yl <- range(1, x$weights, x$undersampling)
  graphics::plot(range(x$prevalence), yl, type = "n", log = "x", xlab = "Prevalence",
                 ylab = "Multiplier of the required size", main = attr(x, "criterion"), ...)
  graphics::abline(h = 1, col = "grey60")
  for (j in seq_along(aucs)) {
    d <- x[x$cstatistic == aucs[j], ]; d <- d[order(d$prevalence), ]
    graphics::lines(d$prevalence, d$weights, col = cols[j], lwd = 2)
    graphics::lines(d$prevalence, d$undersampling, col = cols[j], lwd = 2, lty = 2)
  }
  graphics::legend("topright", c(paste("AUC", aucs), "Class weights", "Undersampling"),
                   col = c(cols, "grey30", "grey30"), lty = c(rep(1, length(aucs)), 1, 2), lwd = 2, bty = "n", cex = 0.8)
  invisible(x)
}
