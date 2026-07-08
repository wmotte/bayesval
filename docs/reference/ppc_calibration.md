# Posterior predictive calibration check at a horizon.

Posterior predictive calibration check at a horizon.

## Usage

``` r
ppc_calibration(
  ps,
  horizon = NULL,
  method = c("auto", "hare", "group"),
  draws = NULL,
  grid = seq(0.02, 0.98, by = 0.02),
  seed = 1L
)
```

## Arguments

- ps:

  a `psurv`.

- horizon:

  index into ps\$horizons (default: last).

- method:

  calibration-curve method passed to .calibration_curve.

- draws:

  optional draw subset/count (PPC is heavier; thin freely).

- grid:

  predicted-risk grid on which observed curves are evaluated for the
  credible band.

- seed:

  RNG seed.

## Value

list with \$ici_obs, \$ici_rep, \$ppp (posterior predictive p-value),
\$band (data.frame grid + observed-curve median/lower/upper), and
\$curve_obs (observed calibration curve on real data, per grid point).
