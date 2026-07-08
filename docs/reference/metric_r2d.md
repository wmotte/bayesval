# Royston-Sauerbrei R^2_D (per draw) — explained variation from the ordered linear predictor. D is the Cox coefficient on scaled normal rankit scores. R^2_D = (D^2/kappa^2) / (pi^2/6 + D^2/kappa^2), kappa = sqrt(8/pi).

Royston-Sauerbrei R^2_D (per draw) — explained variation from the
ordered linear predictor. D is the Cox coefficient on scaled normal
rankit scores. R^2_D = (D^2/kappa^2) / (pi^2/6 + D^2/kappa^2), kappa =
sqrt(8/pi).

## Usage

``` r
metric_r2d(draw, w = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

## Value

c(R2_D = value)
