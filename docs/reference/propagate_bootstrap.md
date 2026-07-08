# Propagate a metric with the Bayesian bootstrap (B2 + B3 in one call). Thin alias for propagate_metric(..., weights = "bootstrap").

Propagate a metric with the Bayesian bootstrap (B2 + B3 in one call).
Thin alias for propagate_metric(..., weights = "bootstrap").

## Usage

``` r
propagate_bootstrap(ps, metric_fn, draws = NULL, seed = 1L, progress = FALSE)
```

## Arguments

- ps:

  a `psurv` object.

- metric_fn:

  function(draw, w) -\> numeric scalar OR named numeric vector.

- draws:

  optional integer subset / count of draws to use (thinning).

- seed:

  RNG seed for the bootstrap weights (reproducibility).

- progress:

  print a dot every 10% of draws.
