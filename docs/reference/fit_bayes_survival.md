# Fit a Bayesian survival model and return a `psurv`.

Fit a Bayesian survival model and return a `psurv`.

## Usage

``` r
fit_bayes_survival(
  formula,
  data,
  horizons,
  backend = c("stan_surv", "stan_pwexp", "coxph", "pwexp", "bart"),
  time = NULL,
  status = NULL,
  S = 500L,
  ...
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

  "stan_pwexp" (primary rstan piecewise-exponential PH), "stan_surv"
  (rstanarm M-spline PH if available), "pwexp" (piecewise-exponential
  Poisson-trick PH; survival-only), "bart" (nonparametric BART
  survival), or "coxph" (Cox-Laplace approximation; fallback / testing).

- time, status:

  column names (needed for the stan_surv prediction step).

- S:

  number of posterior draws (coxph backend) / thinning target.

- ...:

  passed to the backend fitter (priors, chains, cuts, etc.).
