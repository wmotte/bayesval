# Build a psurv from an rstanarm::stan_surv fit (the primary Bayesian PH model).

Uses posterior_survfit() to obtain S_i(t)^(s) and the posterior linear
predictor. Requires rstanarm; kept thin so the model layer is swappable.

## Usage

``` r
psurv_from_stan_surv(fit, newdata, time, status, horizons, draws = NULL)
```

## Arguments

- fit:

  a stan_surv fit.

- newdata:

  evaluation data (defaults to the fit's model frame).

- time, status:

  names of the time / status columns in newdata.

- horizons:

  `numeric[H]` evaluation horizons.

- draws:

  optional number of posterior draws to thin to.
