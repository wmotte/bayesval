# Harrell's C (per draw). Wraps survival::concordance with case weights.

survival::concordance(Surv ~ marker) uses the convention that a *larger*
marker predicts *longer* survival. eta is a risk score (larger = shorter
survival), so we feed -eta; the returned concordance is then Harrell's C
for the risk score (equivalently 1 - concordance(Surv ~ eta)).

## Usage

``` r
metric_harrell_c(draw, w = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

## Value

c(C_harrell = value)
