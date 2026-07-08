# The framework's default settings.

The framework's default settings.

## Usage

``` r
.BC_DEFAULTS
```

## Fields

- `seed`:

  global RNG seed for reproducible draws/bootstraps.

- `horizons`:

  default evaluation horizons t (NULL =\> choose from the data, see
  bc_horizons()).

- `horizons_q`:

  if horizons is NULL, the event-time quantiles used to pick data-driven
  horizons.

- `tau`:

  truncation time for time-dependent measures (NULL =\> the largest
  horizon).

- `S`:

  number of posterior draws (coxph-Laplace backend / thinning target for
  stan_surv).

- `thin`:

  keep every `thin`-th draw when summarising heavy metrics.

- `bins_dcal`:

  number of cells for D-calibration.

- `thresholds`:

  default decision-curve threshold probabilities.

- `ci_probs`:

  credible-interval probabilities (median + 95% CrI).
