# MCMC convergence diagnostics for a fitted Bayesian survival model.

Extracts split-Rhat, bulk/tail ESS and the divergent-transition count.
Works on a stanreg (rstanarm) fit or anything the `posterior` package
can coerce to draws. For the Laplace backends (label starts with "coxph"
/ "pwexp") there is no MCMC and the function returns a row with
applicable = FALSE.

## Usage

``` r
mcmc_diagnostics(fit)
```

## Arguments

- fit:

  a stan_surv/stanreg fit, OR a `psurv` carrying attr("fit").

## Value

list(applicable, max_rhat, min_ess_bulk, min_ess_tail, n_divergent,
table = per-parameter data.frame).
