# BART survival (nonparametric Bayesian) model -\> psurv.

The plan's nonparametric comparator (B1): demonstrates the framework is
model-agnostic by feeding a Bayesian machine-learning survival model
through the SAME `psurv` interface. Requires the `BART` package
(Sparapani et al.); guarded like the stan_surv backend so the framework
still loads without it.

## Usage

``` r
psurv_from_bart(formula, data, horizons, time, status, S = 500L, seed = 1L)
```

## Arguments

- formula:

  Surv(time, status) ~ covariates.

- data:

  the data.frame.

- horizons:

  `numeric[H]` evaluation horizons.

- time, status:

  response column names.

- S:

  posterior draws to thin to.

- seed:

  RNG seed.
