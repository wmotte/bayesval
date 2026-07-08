# Fit a psurv with transparent disk caching.

On a cache hit the stored psurv is returned immediately; on a miss the
model is fitted via fit_bayes_survival() and written to
`<cache_dir>/psurv-<key>.rds`. The key is derived from the inputs (or
supply your own via `key`).

## Usage

``` r
fit_bayes_survival_cached(
  formula,
  data,
  horizons,
  backend = "coxph",
  S = 500L,
  ...,
  cache_dir = ".psurv_cache",
  key = NULL,
  refresh = FALSE,
  verbose = TRUE
)
```

## Arguments

- formula:

  survival formula, Surv(time, status) ~ covariates.

- data:

  training/evaluation data.frame.

- horizons:

  `numeric[H]` evaluation horizons.

- backend:

  "stan_surv" (rstanarm M-spline PH; primary), "pwexp"
  (piecewise-exponential Poisson-trick PH; robustness comparator,
  survival-only), "bart" (nonparametric BART survival), or "coxph"
  (Cox-Laplace approximation; fallback / testing).

- S:

  number of posterior draws (coxph backend) / thinning target.

- ...:

  passed to fit_bayes_survival() (and folded into the key).

- cache_dir:

  directory for cache files (created if needed).

- key:

  explicit cache key (default: hash of the inputs).

- refresh:

  TRUE forces a re-fit even on a cache hit.

- verbose:

  message on hit/miss.

## Value

a `psurv` object.
