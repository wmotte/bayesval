# Piecewise-exponential (Poisson-trick) Bayesian PH model -\> psurv.

A genuine Bayesian proportional-hazards analogue with a
piecewise-constant baseline hazard, fit WITHOUT Stan: split the
follow-up at cut points into person-time pieces, fit a Poisson GLM with
piece-specific intercepts (the log baseline hazards) and covariate
log-HRs, and draw the posterior from the GLM's Laplace approximation
N(coef, vcov). This is the plan's robustness comparator (B1) and shows
the framework is model-agnostic — it produces the identical `psurv`
object as stan_surv or the Cox-Laplace fallback, so B2-B6 are unchanged.
Needs only the `survival` package.

## Usage

``` r
psurv_from_pwexp(
  formula,
  data,
  horizons,
  cuts = NULL,
  n_pieces = 8L,
  S = 500L,
  seed = 1L
)
```

## Arguments

- formula:

  Surv(time, status) ~ covariates.

- data:

  the data.frame.

- horizons:

  `numeric[H]` evaluation horizons.

- cuts:

  interior cut points; NULL =\> event-time quantiles giving `n_pieces`
  pieces.

- n_pieces:

  number of baseline-hazard pieces when `cuts` is NULL.

- S, seed:

  posterior draws and RNG seed.
