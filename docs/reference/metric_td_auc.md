# Time-dependent (cumulative/dynamic) AUC at each horizon, IPCW, per draw. Cases: T_i \<= t & event; controls: T_j \> t. Uses the predicted risk 1 - S_i(t) as the marker and IPCW weights G-hat from the censoring KM.

Time-dependent (cumulative/dynamic) AUC at each horizon, IPCW, per draw.
Cases: T_i \<= t & event; controls: T_j \> t. Uses the predicted risk
1 - S_i(t) as the marker and IPCW weights G-hat from the censoring KM.

## Usage

``` r
metric_td_auc(draw, w = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

## Value

named vector AUC_t for each horizon plus iAUC (trapezoid over t).
