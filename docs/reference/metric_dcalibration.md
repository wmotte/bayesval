# D-calibration (per draw) — A2.6, Haider et al. (JMLR 2020). Distributional calibration of the WHOLE predicted survival distribution, not one horizon: if the model is calibrated, u_i = S_i(T_i) is Uniform(0,1). Bin `[0,1]` into `bins` equal cells; an event contributes 1 to its cell, a censored subject (known only that its true u \< S_i(c_i)) spreads its unit mass uniformly over `[0, S_i(c_i)]` across cells. The Pearson statistic vs the expected count n/bins tests uniformity (chi-square, bins-1 df).

Survival at each subject's own event/censoring time is interpolated from
the horizon-knot cumulative hazard, so D-calibration honours the full
curve even though the psurv array is stored at a few horizons.

## Usage

``` r
metric_dcalibration(draw, w = NULL, bins = 10L)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

- bins:

  number of equal-probability cells (default 10).

## Value

c(Dcal = chi-square statistic, Dcal_p = its p-value).
