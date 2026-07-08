# Bayesian-bootstrap replicates of a metric on a SINGLE fixed prediction set.

Isolates evaluation-set (B3) uncertainty when the model is treated as
fixed (e.g. the true model in the coverage study): the metric is
recomputed under B fresh Dirichlet(1) weightings of the n subjects.

## Usage

``` r
bb_replicates(draw, metric_fn, B = 500L, seed = 1L)
```

## Arguments

- draw:

  a psurv_draw()-style list (fixed predictions + data).

- metric_fn:

  function(draw, w) -\> scalar or named vector.

- B:

  number of bootstrap replicates.

- seed:

  RNG seed.

## Value

matrix `[B x K]` of metric values.
