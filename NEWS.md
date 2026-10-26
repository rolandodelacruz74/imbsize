# imbsize 0.3.0

* New `ss_recal_precision()`: size for the slope of leave-one-out Platt scaling to lie in [1 - epsilon, 1 + epsilon] with a given probability, n = A_unw (z/epsilon)^2, using tau^2 = A (exact for the unweighted fit, numerical for balanced weights). It binds when there are few columns; take the larger of it and the C2 size.
* Documentation: "exact" criteria are now called fixed-point criteria.

# imbsize 0.2.0

* The limits are evaluated at kappa = (p - 1)/n: the leading terms are fixed-dimension quantities driven by the p - 1 coordinates orthogonal to the best direction. This removes most of the conservativeness with few columns; sizes change by the factor (p - 1)/p.
* Probability version of C1: the standard deviation of the slope is cs(kappa) * tau / sqrt(n), which accounts for the inflation of the fit; evaluated directly in simulation, it attains the requested probability in Gaussian designs.
* The exact criteria are recommended; the closed forms bracket them (resummation below, linear form above).

# imbsize 0.1.0

* First version: criteria C1 and C2 (exact and closed form, probability versions), expected performance at a given size, cost of rebalancing, separability threshold.
* Exact criteria solved as one system (fixed-point equations plus the target) by Newton, with safeguarded Newton for the proximal point and a Gauss-Hermite rule (option `imbsize.nodes`, default 120): about 0.3 s per method.
