# imbsize 0.1.0

* First version: criteria C1 and C2 (exact and closed form, probability versions), expected performance at a given size, cost of rebalancing, separability threshold.
* Exact criteria solved as one system (fixed-point equations plus the target) by Newton, with safeguarded Newton for the proximal point and a Gauss-Hermite rule (option `imbsize.nodes`, default 120): about 0.3 s per method.
