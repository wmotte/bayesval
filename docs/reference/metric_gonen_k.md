# Gonen-Heller K (per draw) — model-based concordance from eta only, so it is independent of the study censoring distribution. Closed form (A1): `K = sum_{i!=j} w_i w_j I(eta_i > eta_j) / (1 + exp(eta_j - eta_i))` normalised by the pair weight mass. O(n^2); intended for moderate n.

Gonen-Heller K (per draw) — model-based concordance from eta only, so it
is independent of the study censoring distribution. Closed form (A1):
`K = sum_{i!=j} w_i w_j I(eta_i > eta_j) / (1 + exp(eta_j - eta_i))`
normalised by the pair weight mass. O(n^2); intended for moderate n.

## Usage

``` r
metric_gonen_k(draw, w = NULL, max_n = 3000L)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

- max_n:

  guard: above this n the pair loop is subsampled (with a note in the
  returned attribute) to keep the posterior loop tractable.

## Value

c(K_gonen = value)
