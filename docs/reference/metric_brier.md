# IPCW Brier score at each horizon (per draw). BS(t) = sum_i w_i \* Ŵ_i(t) \* (I(Ytilde_i \> t) - S_i(t))^2 Ŵ_i(t): 1/G(T_i-) for events before t, 1/G(t) for survivors past t, 0 for subjects censored before t.

IPCW Brier score at each horizon (per draw). BS(t) = sum_i w_i \* Ŵ_i(t)
\* (I(Ytilde_i \> t) - S_i(t))^2 Ŵ_i(t): 1/G(T_i-) for events before t,
1/G(t) for survivors past t, 0 for subjects censored before t.

## Usage

``` r
metric_brier(draw, w = NULL, Ghat = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

- Ghat:

  optional precomputed censoring KM (from .km_censoring); supplying it
  avoids recomputing the fixed censoring distribution on every draw.

## Value

named vector Brier_t plus IBS (integrated over horizons).
