# imbsize

<!-- badges: start -->
[![R-CMD-check](https://github.com/rolandodelacruz74/imbsize/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rolandodelacruz74/imbsize/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Sample size for developing logistic risk prediction models fitted without rebalancing, with class weights or with
random undersampling.

- **C1** requires an expected calibration slope of at least *S* (models used without recalibration).
- **C2** requires that the AUC lost by estimating a model that will be recalibrated be at most *δ*.

Both come in exact form (fixed point of the proportional-asymptotic limits under a Gaussian model) and closed form,
with probability ("assurance") versions, and express the cost of rebalancing as a multiplier of the required size.

```r
# install.packages("remotes"); remotes::install_github("rolandodelacruz74/imbsize")
library(imbsize)
imbsize(cstatistic = 0.81, prevalence = 0.0916, parameters = 71)
```

The criteria assume Gaussian predictors given the outcome, a correctly specified logistic model and independent
patients. With mostly binary predictors and rare categories they can recommend too few patients; simulate on pilot
covariates when they are available.

Reference: de la Cruz R. Sample-size criteria for logistic risk models that will be recalibrated or rebalanced
(manuscript).
