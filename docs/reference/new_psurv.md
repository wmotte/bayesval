# Construct a posterior survival object.

Construct a posterior survival object.

## Usage

``` r
new_psurv(surv, eta, time, status, horizons, x = NULL, label = "model")
```

## Arguments

- surv:

  array `[S x n x H]` of predicted survival S_i(t_h)^(s) in (0,1\].

- eta:

  matrix `[S x n]` of linear predictors eta_i^(s).

- time:

  `numeric[n]` observed follow-up time Ytilde_i = min(T_i, C_i).

- status:

  `integer[n]` event indicator Delta_i (1 = event, 0 = censored).

- horizons:

  `numeric[H]` the evaluation horizons t_h (same order as dim 3).

- x:

  optional data.frame`[n x p]` of covariates (needed by some
  calibration/PPC steps that refit on eta).

- label:

  character tag identifying the model that produced the draws.

## Value

an object of class "psurv".
