# Posterior net-benefit curves over draws (B6).

Posterior net-benefit curves over draws (B6).

## Usage

``` r
posterior_net_benefit(
  ps,
  thresholds = seq(0.05, 0.6, by = 0.05),
  horizon = NULL,
  weights = c("none", "bootstrap"),
  draws = NULL,
  seed = 1L
)
```

## Arguments

- ps:

  a `psurv`.

- thresholds:

  risk thresholds p_t.

- horizon:

  index into ps\$horizons (default last).

- weights:

  "none" or "bootstrap" (adds evaluation-sample uncertainty).

- draws:

  optional draw subset/count.

- seed:

  RNG seed.

## Value

list with \$thresholds, \$model (band df), \$all (band df), and
\$p_useful = P(NB_model \> max(NB_all, 0)) per threshold.
