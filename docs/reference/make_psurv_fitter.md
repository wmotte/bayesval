# Build a fit_psurv closure for bb_optimism(): function(data, weights) -\> psurv.

The Laplace/Cox backend accepts Dirichlet case weights directly, so the
Bayesian-bootstrap optimism loop reweights the fit coherently.

## Usage

``` r
make_psurv_fitter(formula, horizons, S = 300L, backend = c("coxph"))
```

## Arguments

- formula, horizons, S:

  as in fit_bayes_survival.

- backend:

  "coxph" (weighted refit) is the supported optimism backend.
