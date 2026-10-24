# imbsize

Sample size for developing logistic risk prediction models fitted without rebalancing, with class weights or with
random undersampling.

- **C1** requires an expected calibration slope of at least *S* (models used without recalibration).
- **C2** requires that the AUC lost by estimating a model that will be recalibrated be at most *δ*.

Both come in exact form (fixed point of the proportional-asymptotic limits under a Gaussian model) and closed form,
with probability ("assurance") versions, and express the cost of rebalancing as a multiplier of the required size.

```r
# install.packages("remotes"); remotes::install_github("<user>/imbsize")
library(imbsize)
imbsize(cstatistic = 0.81, prevalence = 0.0916, parameters = 71)
```

The criteria assume Gaussian predictors given the outcome, a correctly specified logistic model and independent
patients. With mostly binary predictors and rare categories, add about 20% or simulate on pilot covariates.

Reference: de la Cruz R. Sample-size criteria for the calibration and discrimination of logistic risk models fitted
with or without rebalancing (manuscript).
