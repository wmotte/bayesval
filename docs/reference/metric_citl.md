# Calibration-in-the-large at each horizon (per draw): mean predicted risk vs observed (weighted Kaplan-Meier) risk.

Calibration-in-the-large at each horizon (per draw): mean predicted risk
vs observed (weighted Kaplan-Meier) risk.

## Usage

``` r
metric_citl(draw, w = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

## Value

named vector CITL_t = predicted_mean_risk - observed_risk.
