# External-validation psurv from a fitted stan_pwexp model on new data.

Uses the training posterior draws of beta and the baseline log-hazards;
only the covariate design changes (train on cohort A, validate on cohort
B, B3).

## Usage

``` r
stan_pwexp_predict(obj, newdata, time, status, horizons = NULL)
```

## Arguments

- obj:

  a fitted `stan_pwexp` object from fit_stan_pwexp().

- newdata:

  the external validation data.frame.

- time, status:

  names of the time / status columns in `newdata`.

- horizons:

  `numeric[H]` evaluation horizons (default: the training horizons).
