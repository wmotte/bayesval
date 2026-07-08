# Index of Prediction Accuracy (scaled Brier) at each horizon (per draw): IPA(t) = 1 - BS_model(t) / BS_null(t), null = marginal Kaplan-Meier.

Index of Prediction Accuracy (scaled Brier) at each horizon (per draw):
IPA(t) = 1 - BS_model(t) / BS_null(t), null = marginal Kaplan-Meier.

## Usage

``` r
metric_ipa(draw, w = NULL, Ghat = NULL)
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

named vector IPA_t.
