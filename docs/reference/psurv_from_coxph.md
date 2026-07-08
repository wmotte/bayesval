# Fallback constructor: posterior draws from a frequentist Cox fit.

This is the framework's *testing harness* and a legitimate large-sample
(Laplace / Bayesian-central-limit) approximation to the posterior of a
Cox model under flat priors: coefficients are drawn from N(beta_hat,
vcov_hat) and combined with the Breslow baseline cumulative hazard. It
lets the entire Part B pipeline (B2-B6) run and be unit-tested with only
the `survival` package installed — no Stan toolchain required. The
primary manuscript path remains psurv_from_stan_surv(); this documents
the model-agnostic claim.

## Usage

``` r
psurv_from_coxph(
  formula,
  data,
  horizons,
  S = 500L,
  baseline = c("fixed", "bootstrap"),
  seed = 1L
)
```

## Arguments

- formula:

  a survival::coxph formula, e.g. Surv(time, status) ~ x1 + x2.

- data:

  the evaluation data.frame.

- horizons:

  `numeric[H]` evaluation horizons.

- S:

  number of posterior draws.

- baseline:

  "fixed" (Breslow at the MLE) or "bootstrap" (Bayesian bootstrap of the
  baseline hazard for extra coherence).

- seed:

  RNG seed for reproducibility.
