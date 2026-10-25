#' Sample size for a logistic risk model fitted with or without rebalancing
#'
#' Computes criteria C1 (calibration slope of a model used without recalibration) and C2 (loss of AUC of a model that
#' will be recalibrated) for the unweighted fit, class weights and random undersampling.
#'
#' The criteria assume Gaussian predictors given the outcome, a correctly specified logistic model and independent
#' patients. With mostly binary predictors and rare categories they can recommend too few patients (5--17\% too few
#' for C1 in two clinical data sets); simulation on pilot covariates is advisable when they are available. The criteria of Riley et al. (2019)
#' on the optimism of \eqn{R^2} and on the precision of the overall risk still apply.
#'
#' @inheritParams ss_c1
#' @inheritParams ss_c2
#' @return An object of class \code{imbsize}: a data frame with one row per criterion and method.
#' @examples
#' imbsize(cstatistic = 0.81, prevalence = 0.0916, parameters = 71, form = "closed")
#' @export
imbsize <- function(cstatistic, prevalence, parameters, method = c("unweighted", "weights", "undersampling"),
                    slope = 0.9, auc_loss = 0.01, assurance = NULL, form = c("exact", "closed")) {
  form <- match.arg(form)
  out <- rbind(ss_c1(cstatistic, prevalence, parameters, method, slope, assurance, form),
               ss_c2(cstatistic, prevalence, parameters, method, auc_loss, assurance, form))
  attr(out, "inputs") <- list(cstatistic = cstatistic, prevalence = prevalence, parameters = parameters, slope = slope,
                              auc_loss = auc_loss, assurance = assurance, form = form)
  class(out) <- c("imbsize", "data.frame"); out
}

#' @export
print.imbsize <- function(x, ...) {
  a <- attr(x, "inputs")
  cat(sprintf("Sample size for a logistic risk model (%s criteria)\n", a$form))
  cat(sprintf("Inputs: C-statistic %.3f, prevalence %.4f, %d predictor columns%s\n", a$cstatistic, a$prevalence,
              as.integer(a$parameters), if (is.null(a$assurance)) "" else sprintf(", assurance %.2f", a$assurance)))
  cat(sprintf("C1: expected calibration slope >= %.2f (model used without recalibration)\n", a$slope))
  cat(sprintf("C2: loss of AUC <= %.3f (model that will be recalibrated)\n\n", a$auc_loss))
  print(as.data.frame(unclass(x))[, c("criterion", "method", "n", "events", "EPV", "multiplier")], row.names = FALSE)
  cat("\nWith mostly binary or rare predictors the criteria may recommend too few patients; simulate on pilot covariates if available.\n")
  invisible(x)
}
