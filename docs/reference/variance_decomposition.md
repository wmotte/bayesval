# Decompose posterior variance into parameter vs evaluation-sample sources.

Runs the same metric three ways and reports interval widths so the
coverage study (Part C) can decide whether to down-weight one source:

- parameter-only (weights = "none")

- bootstrap-only (fix the draw index, vary only Dirichlet weights)

- combined (weights = "bootstrap")

## Usage

``` r
variance_decomposition(ps, metric_fn, draws = NULL, seed = 1L)
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

## Value

data.frame of source, sd, and 95% interval width per metric column.
