# bayesval

Bayesian probabilistic validation of Cox proportional-hazards prediction
models. Every discrimination, calibration, and clinical-utility metric you
would normally compute once — Harrell/Uno/Gonen-Heller/Antolini concordance,
time-dependent AUC, calibration-in-the-large, calibration slope, ICI/E50/E90,
D-calibration, IPCW Brier/IBS/IPA, R²_D, decision-curve net benefit — is
instead reported as a **posterior distribution** with a 95% credible interval.

The single idea behind the package: every validation metric `M` is a
functional `M(S, data)` of a posterior survival object and the data. Evaluate
the existing frequentist estimator at each posterior draw; the collection of
per-draw values **is** the metric's posterior. No new point estimators —
`propagate_metric()` is the only wrapper. Two uncertainty sources are
propagated coherently: model-parameter uncertainty (posterior draws) and
evaluation-sample uncertainty (a Bayesian bootstrap over the validation set).

## Installation

```r
# install.packages("remotes")
remotes::install_github("wmotte/bayesval")
```

The modelling core needs only `survival`; everything below works without a
Stan toolchain via the `coxph` backend. The full-Bayes path (`rstan`), BART
comparator (`BART`), and a few optional calibration/decision-curve estimators
are Suggested, not required — install them only if you use that backend.

## Quick start

```r
library(bayesval)

set.seed(1); n <- 500
dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4))
dat$time   <- rexp(n, 0.15 * exp(0.8 * dat$x1 - 0.5 * dat$x2))
dat$status <- rbinom(n, 1, 0.7)
hz <- quantile(dat$time[dat$status == 1], c(.25, .5, .75))

ps <- fit_bayes_survival(survival::Surv(time, status) ~ x1 + x2, dat, hz,
                         backend = "coxph", S = 300)           # survival-only

disc <- propagate_many(ps, list(harrell  = metric_harrell_c,
                                uno      = metric_uno_c,
                                antolini = metric_antolini_c))
posterior_table(disc)                                          # median [95% CrI]

propagate_metric(ps, metric_harrell_c, weights = "bootstrap")  # + eval-set uncertainty
ppc_calibration(ps, horizon = 2)$ppp                            # calibration PPC
posterior_net_benefit(ps)$p_useful                              # decision-curve layer
```

## Model backends

Every backend returns the same model-agnostic `psurv` object (an `[S x n x
H]` posterior survival array plus the linear predictor), so every metric
above is backend-agnostic:

| backend | needs | role |
|---|---|---|
| `stan_pwexp` | `rstan` | **primary** full-Bayes piecewise-exponential PH (MCMC) |
| `pwexp` | `survival` | piecewise-exponential Poisson-trick PH (GLM + Laplace) |
| `coxph` | `survival` | Cox-Laplace approximation (fast fallback / no Stan) |
| `bart` | `BART` | nonparametric Bayesian survival comparator |
| `stan_surv` | `rstanarm` w/ `stan_surv` | M-spline PH (aspirational; not yet in CRAN rstanarm) |

MCMC fits are gated by `check_mcmc_gate()` (Rhat < 1.01, ESS, no divergences)
before any metric is reported; `prior_sensitivity()` refits under several
log-hazard-ratio prior scales to check robustness.

## Documentation

`?propagate_metric` is the entry point for writing a new metric; every
`metric_*()` function follows the same `function(draw, w)` signature and can
be passed to it directly. See the function reference (`help(package =
"bayesval")`) for the full API.

## Citation and reproducibility

`bayesval` is the software companion to:

> Otte WM, Kampman SL, [...] *A Bayesian framework for validating survival
> prediction models: posterior credible intervals for discrimination,
> calibration, and clinical utility.* in prep.

## License

MIT © Willem M. Otte
