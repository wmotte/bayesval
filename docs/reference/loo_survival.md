# PSIS-LOO elpd from a pointwise log-likelihood matrix `[S x n]`.

PSIS-LOO elpd from a pointwise log-likelihood matrix `[S x n]`.

## Usage

``` r
loo_survival(log_lik)
```

## Arguments

- log_lik:

  matrix draws x subjects, e.g. rstanarm::log_lik(fit).

## Value

the loo object (contains elpd_loo, p_loo, pareto-k diagnostics).
