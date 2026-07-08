# Net benefit across thresholds for a single draw.

Net benefit across thresholds for a single draw.

## Usage

``` r
net_benefit_draw(draw, thresholds, horizon = NULL, w = NULL)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- thresholds:

  numeric vector of risk thresholds p_t.

- horizon:

  index into `draw$horizons` (default: last).

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

## Value

list(model = `numeric[T]`, all = `numeric[T]`) net benefit per
threshold.
