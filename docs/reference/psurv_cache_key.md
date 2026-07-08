# Deterministic cache key for a psurv fit from its defining inputs.

Deterministic cache key for a psurv fit from its defining inputs.

## Usage

``` r
psurv_cache_key(formula, data, horizons, backend, S, extra = NULL)
```

## Arguments

- formula, data, horizons, backend, S:

  the fit specification.

- extra:

  any additional backend arguments that change the result (priors,
  n_pieces, seed, ...).

## Value

a short hex string.
