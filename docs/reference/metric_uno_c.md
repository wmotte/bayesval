# Uno's C (IPCW, truncated at tau) per draw. survival::concordance with timewt = "n/G2" is the Uno estimator; ymax = tau truncates.

Uno's C (IPCW, truncated at tau) per draw. survival::concordance with
timewt = "n/G2" is the Uno estimator; ymax = tau truncates.

## Usage

``` r
metric_uno_c(draw, w = NULL, tau = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

- tau:

  truncation time; defaults to the largest horizon.

## Value

c(C_uno = value)
