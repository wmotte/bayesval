# Calibration slope (per draw): coefficient of a Cox model with eta as the sole covariate. slope \< 1 signals overfitting (van Houwelingen shrinkage).

Calibration slope (per draw): coefficient of a Cox model with eta as the
sole covariate. slope \< 1 signals overfitting (van Houwelingen
shrinkage).

## Usage

``` r
metric_cal_slope(draw, w = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

## Value

c(cal_slope = value)
