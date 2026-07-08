# External-validation psurv: evaluate a coxph's Laplace posterior on NEW data.

Fits nothing new — takes the training coxph, draws coefficients from its
Laplace posterior N(beta_hat, V) and uses its Breslow baseline, then
predicts S_i(t)^(s) and eta_i^(s) on `newdata`. This is how B3's
external-set validation design is realised (train on cohort A, validate
on cohort B).

## Usage

``` r
psurv_external(fit, newdata, time, status, horizons, S = 500L, seed = 1L)
```

## Arguments

- fit:

  a survival::coxph fit (from the training cohort).

- newdata:

  external validation data.frame (with the model covariates).

- time, status:

  column names of the response in newdata.

- horizons:

  `numeric[H]` evaluation horizons.

- S, seed:

  posterior draws and RNG seed.
