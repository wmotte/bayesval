# Convenience: propagate several named metric functions in one pass over the draws (cheaper than one pass each — the survival object is only sliced once).

Convenience: propagate several named metric functions in one pass over
the draws (cheaper than one pass each — the survival object is only
sliced once).

## Usage

``` r
propagate_many(
  ps,
  metric_fns,
  weights = c("none", "bootstrap"),
  draws = NULL,
  seed = 1L,
  progress = FALSE
)
```

## Arguments

- ps:

  a `psurv`.

- metric_fns:

  named list of function(draw, w) -\> scalar/named-vector.

- weights:

  "none" (parameter uncertainty only), "bootstrap" (add a Dirichlet(1)
  reweight per draw), or a user matrix `[S x n]`.

- draws:

  optional integer subset / count of draws to use (thinning).

- seed:

  RNG seed for the bootstrap weights (reproducibility).

- progress:

  print a dot every 10% of draws.
