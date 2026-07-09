# Bayesian-bootstrap optimism correction for a scalar performance estimator.

Bayesian-bootstrap optimism correction for a scalar performance
estimator.

## Usage

``` r
bb_optimism(
  fit_psurv,
  data,
  metric_fn,
  B = 50L,
  summarise = stats::median,
  seed = 1L
)
```

## Arguments

- fit_psurv:

  function(data, weights) -\> a `psurv` fitted on `data` (weights are
  Dirichlet case weights; pass to the fitter or ignore if it resamples
  internally).

- data:

  the evaluation data.frame.

- metric_fn:

  a per-draw metric function metric_fn(draw, w) -\> scalar
  (single-value; if it returns a vector the first element is used).

- B:

  number of Bayesian-bootstrap replicates.

- summarise:

  function to collapse a metric's posterior draws to a point (default
  median).

- seed:

  RNG seed.

## Value

list with \$apparent (on full data), \$optimism (mean and per-rep), and
\$corrected = apparent - mean(optimism).
