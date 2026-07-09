# ICI / E50 / E90 at a horizon (per draw): mean / 50th / 90th percentile of \|pi_i(t) - pi_hat^cal_i(t)\| (Austin, Harrell & van Klaveren 2020).

ICI / E50 / E90 at a horizon (per draw): mean / 50th / 90th percentile
of \|pi_i(t) - pi_hat^cal_i(t)\| (Austin, Harrell & van Klaveren 2020).

## Usage

``` r
metric_ici(draw, w = NULL, horizon = NULL, method = c("auto", "hare", "group"))
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

- horizon:

  index into draw\$horizons (default: last).

- method:

  "auto" (hare if installed else group), "hare", or "group". The grouped
  method is fully case-weighted. HARE smoothing is unweighted because
  [`polspline::hare()`](https://rdrr.io/pkg/polspline/man/hare.html) has
  no case-weight argument; the final ICI/E50/E90 aggregation still uses
  `w`.

## Value

c(ICI = , E50 = , E90 = )
