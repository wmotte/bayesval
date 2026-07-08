# Prior-sensitivity analysis (checklist 1): refit the primary full-Bayes PH model (fit_stan_pwexp) under several log-HR prior SDs and compare the posterior of each validation metric. A robust conclusion should be stable across weakly-informative priors; a large shift flags prior-dependence.

Uses `beta_scale` (the Normal(0, beta_scale) SD on the
log-hazard-ratios) as the prior knob — the working analogue of
stan_surv's coefficient prior, since CRAN rstanarm ships no stan_surv.

## Usage

``` r
prior_sensitivity(
  formula,
  data,
  horizons,
  metric_fns,
  scales = c(tight = 1, default = 2.5, loose = 5),
  flag = 0.02,
  ...
)
```

## Arguments

- formula, data, horizons:

  as in fit_stan_pwexp.

- metric_fns:

  named list of per-draw metric functions to compare.

- scales:

  named numeric of prior SDs; default contrasts the plan default (2.5)
  with a tight (1.0) and a loose (5.0) log-HR prior.

- flag:

  threshold on the posterior-median range across priors that triggers
  the prior-sensitive message (default 0.02).

- ...:

  extra args passed to fit_stan_pwexp (chains, iter, n_pieces, ...).

## Value

data.frame: prior x metric x (median, lower, upper), plus a printed note
if any metric's median moves by more than `flag` across priors.
