# Propagate a per-draw metric across the posterior (B2), optionally with the Bayesian bootstrap (B3).

Propagate a per-draw metric across the posterior (B2), optionally with
the Bayesian bootstrap (B3).

## Usage

``` r
propagate_metric(
  ps,
  metric_fn,
  weights = c("none", "bootstrap"),
  draws = NULL,
  seed = 1L,
  progress = FALSE
)
```

## Arguments

- ps:

  a `psurv` object.

- metric_fn:

  function(draw, w) -\> numeric scalar OR named numeric vector, where
  `draw` is the list from psurv_draw() and `w` is a length-n weight
  vector (sums to 1). Metric functions may ignore `w` if they cannot use
  case weights.

- weights:

  "none" (parameter uncertainty only), "bootstrap" (add a Dirichlet(1)
  reweight per draw), or a user matrix `[S x n]`.

- draws:

  optional integer subset / count of draws to use (thinning).

- seed:

  RNG seed for the bootstrap weights (reproducibility).

- progress:

  print a dot every 10% of draws.

## Value

a `posterior_metric` object: matrix `[S_used x K]` of draws plus the
metric names; summarise with summary(). The summary reports non-finite
draw counts and flags metrics with more than 5% non-finite draws by
default.
