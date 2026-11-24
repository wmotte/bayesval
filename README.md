# bayesval

`bayesval` propagates a range of standard performance measures over posterior draws from the parameters of a fitted Bayesian model, for assessment on a finite sample. Its main focus is the **external validation** of a survival model. Given a Cox proportional-hazards model fitted to other patients, one wants to know its discrimination, calibration and decision utility in a further sample of patients, often quite small and so with few events. Each standard finite-sample statistic is computed in the usual way over each posterior draw, resulting in posteriors and credible intervals (CIs) for all quantities of interest, in place of single values. Hence assessment against a pre-specified deployment gate for a given clinical setting reduces to one posterior probability of acceptability.

The same propagation extends naturally to apparent performance and internal and internal-external validation, all supported. However, **deciding on the suitability of a model generated from another source and patient population for local application, that is, external validation, is the primary motivation**, and all examples shown here externally validate a model developed elsewhere.

All of the following can be assessed: Harrell C-index, Uno, Gonen-Heller, Antolini statistics; time-dependent area under ROC curve; calibration-in-the-large; calibration slope; ICI; E50, E90; D-calibration; IPCW Brier score, integrated Brier score and IPA; R²_D; decision-curve net benefit. As above, each will be represented by a posterior distribution with credible interval. In general, any such statistic `M` of interest can be written as a functional `M(S, data)` of a posterior survival object and the validation data. So simply evaluating the standard estimator over each posterior draw in turn gives a set of posterior draws of the desired statistic.

There will therefore be posterior uncertainty due to uncertainty over the underlying model parameters, as well as due to sampling variability of the validation set, implemented as a Bayesian bootstrap. These two sources of variation are maintained separately. There is no proposal of any novel point estimators; standard statistics are merely evaluated over the posterior. See `propagate_metric()`, the propagation wrapper.

## Status

`bayesval` is currently a research code package, expected to evolve further, and accompanies the survival-validation manuscript listed below. All fast examples can be run using just the `coxph` backend, with no dependence on Stan or other MCMC sampling. Longer examples requiring full MCMC sampling use the `stan_pwexp` backend and require a working `rstan` toolchain; run times may be several minutes.

## Installation

Install the package from GitHub:

```r
# install.packages("remotes")
remotes::install_github("wmotte/bayesval")
```

Only imports are required for the modelling core; in addition to `survival`, including the tidyverse set `dplyr`, `purrr`, `tibble` and `ggplot2`. Optional backends to propagate from, or comparators, require further dependencies, as listed in Suggested:

```r
install.packages(c("posterior", "loo", "BART", "riskRegression", "pec",
                   "timeROC", "survC1", "CPE", "dcurves", "ggdist"))

# Optional, platform-sensitive: needed only for the full Bayesian Stan backend.
install.packages("rstan")
```

It would normally be expected that only relevant ones would be installed, according to the desired backends and metric families. Setting up dependencies to allow Stan models will be the longest task, although of course this can be skipped entirely, as for the quick start below, which uses only the survival-based `coxph` path.

## Quick start

```r
library(bayesval)

set.seed(1); n <- 500
dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4))
dat$time   <- rexp(n, 0.15 * exp(0.8 * dat$x1 - 0.5 * dat$x2))
dat$status <- rbinom(n, 1, 0.7)
hz <- as.numeric(quantile(dat$time[dat$status == 1], c(.25, .5, .75)))

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

`posterior_table()` returns one line for each statistic assessed, giving mean and median over the posterior, credible interval, draw counts and whether reliable. If generated using `weights = "bootstrap"`, uncertainty will be due to variation both over underlying parameters and over the validation sample, not parameters alone. `ppc_calibration()` gives a posterior predictive check p-value for calibration at the requested horizon.

## Model backends

A common, model-agnostic object `psurv` is returned from any backend: posterior draws of survival probability, dimension `[S x n x H]`, and corresponding linear predictor. All metrics are then calculated from this common representation, independent of backend.

| backend | needs | use when |
|---|---|---|
| `coxph` | `survival` | you want fast examples, teaching code, or a no-Stan fallback |
| `pwexp` | `survival` | you want a piecewise-exponential Poisson-trick PH model with a Laplace approximation |
| `stan_pwexp` | `rstan` | you want the primary full-likelihood Bayesian piecewise-exponential PH workflow |
| `bart` | `BART` | you want a nonparametric Bayesian survival comparator |
| `stan_surv` | `rstanarm` with `stan_surv` | aspirational/conditional; not available in the CRAN `rstanarm` release |

MCMC fits must pass convergence checks (`check_mcmc_gate()`; split-Rhat, ESS, number of divergent transitions) before any metric is reported. The effect of using more or less informative priors can be investigated by `prior_sensitivity()`, refitting with a range of prior scales for log-HRs.

## Documentation

The pkgdown site is the main reference:

- package site: <https://wmotte.github.io/bayesval/>
- function index: <https://wmotte.github.io/bayesval/reference/>
- local entry point for writing a new metric: `?propagate_metric`
- full package help: `help(package = "bayesval")`

The `summary()` method for posterior metric objects gives median and CI over all finite values, number finite out of total number of draws, and flags the result as unreliable if greater than 5% non-finite.

Worked examples, mostly external validations on public survival cohorts, can be found in `vignettes/`:

- `worked_example.R` — Rotterdam → GBSG breast cancer (the headline example)
- `showcase_pbc.R` — Mayo PBC, randomized → non-randomized patients
- `temporal_validation.R` — one Rotterdam-trained model on a later Rotterdam era against GBSG
- `worked_example_lung.R` — NCCTG lung, split by recruiting centre
- `worked_example_colon.R` — colon adjuvant trial, split by treatment arm

`showcase_small_cohort.R` applies criteria for acceptable deployment to an artificially generated small sample, knowing the true underlying model. Heavier Stan articles (`worked_example_stan.R`, `small_sample_stan_vs_laplace.R`) involve more extensive Stan sampling and are therefore statically pre-generated, rather than regenerated each time the pkgdown site is built (see `articles/README.md`).

## Development checks

For a lightweight local check without rebuilding vignettes:

```sh
R CMD build --no-build-vignettes --no-manual .
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --ignore-vignettes \
  bayesval_0.1.0.tar.gz
```

If the heavier Suggested dependencies are not installed, it is useful to set
`_R_CHECK_FORCE_SUGGESTS_=false`. All Suggested packages should be installed prior to a full CRAN-style check.

## Citation and reproducibility

`bayesval` is the software companion to:

> Otte WM, Kampman SL, [...] *A Bayesian framework for validating survival
> prediction models: posterior credible intervals for discrimination,
> calibration, and clinical utility.* in prep.

## License

GPL-3. The full license text is in [LICENSE.md](LICENSE.md).
