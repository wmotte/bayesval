# Fit a full-Bayes piecewise-exponential PH model and return a rich object.

Fit a full-Bayes piecewise-exponential PH model and return a rich
object.

## Usage

``` r
fit_stan_pwexp(
  formula,
  data,
  horizons,
  n_pieces = 8L,
  beta_scale = 2.5,
  chains = 4L,
  iter = 1000L,
  warmup = floor(iter/2),
  cores = chains,
  seed = 1L,
  S = NULL,
  refresh = 0L,
  adapt_delta = 0.95
)
```

## Arguments

- formula:

  Surv(time, status) ~ covariates.

- data:

  training data.frame.

- horizons:

  `numeric[H]` evaluation horizons.

- n_pieces:

  baseline-hazard pieces (cuts at event-time quantiles).

- beta_scale:

  prior SD on log-HRs (Normal(0, beta_scale)); plan default 2.5.

- chains, iter, warmup, cores, seed:

  rstan sampling controls.

- S:

  optional thinning to S posterior draws for the psurv array.

- refresh:

  rstan::sampling() progress-print interval (0 = silent).

- adapt_delta:

  target HMC acceptance rate (higher = smaller steps, fewer divergences,
  slower sampling); passed to rstan's `control` list.

## Value

list with \$psurv (internal), \$fit (stanfit), \$beta_draws,
\$log_h0_draws, \$cuts, \$edges, \$terms/\$xlevels/\$contrasts (for
external prediction), and \$diagnostics (max Rhat, min ESS,
divergences).
