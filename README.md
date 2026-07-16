# bayesval

**Package website & documentation:** <https://wmotte.github.io/bayesval/index.html>

`bayesval` is an R package for the **external validation** of survival
prediction models: taking a Cox proportional-hazards model developed elsewhere
and asking how well it discriminates, calibrates, and supports clinical
decisions in a new — often small and event-poor — validation cohort. It keeps
the familiar finite-sample validation estimators that clinicians and
methodologists already use, but evaluates each one across the posterior survival
draws of a Bayesian model, so every performance measure is reported as a
posterior distribution with a credible interval instead of a single point
estimate. A pre-specified deployment gate then reads off as one posterior
probability that the model is fit for local use.

The same propagation works for apparent, internal, and internal-external
designs, but **external validation — deciding whether a model developed
elsewhere can be trusted in a new setting — is the problem the package is built
around**, and the one every worked example below illustrates.

This covers the full validation toolkit — Harrell/Uno/Gonen-Heller/Antolini
concordance, time-dependent AUC, calibration-in-the-large, calibration slope,
ICI/E50/E90, D-calibration, IPCW Brier/IBS/IPA, R²_D, and decision-curve net
benefit — each reported as a posterior distribution with a credible interval.
The core idea: every performance estimator `M` is treated as a functional
`M(S, data)` of a posterior survival object and the validation data. Evaluate
the existing estimator at each posterior draw; the collection of per-draw values
is the posterior distribution for that finite-sample estimator.

Two uncertainty sources are kept explicit: model-parameter uncertainty from the
posterior survival draws, and evaluation-sample uncertainty from a Bayesian
bootstrap over the validation set. No new point estimators are introduced;
`propagate_metric()` is the propagation wrapper.

## Status

`bayesval` is a research package under active development and accompanies the
survival-validation manuscript listed below. The fast examples use the
`coxph` backend and run without Stan. Full Bayesian MCMC examples use the
`stan_pwexp` backend and require a working `rstan` toolchain; they may take
several minutes.

## Installation

Install the package from GitHub:

```r
# install.packages("remotes")
remotes::install_github("wmotte/bayesval")
```

The modelling core needs only the imported package dependencies, including
`survival`, `dplyr`, `purrr`, `tibble`, and `ggplot2`. Optional backends and
comparators need additional Suggested packages:

```r
install.packages(c("posterior", "loo", "BART", "riskRegression", "pec",
                   "timeROC", "survC1", "CPE", "dcurves", "ggdist"))

# Optional, platform-sensitive: needed only for the full Bayesian Stan backend.
install.packages("rstan")
```

Install optional packages only for the backend or metric family you intend to
use. Stan setup is the slowest part of the toolchain; the quick start below uses
only the survival-based `coxph` path.

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

`posterior_table()` returns one row per metric with the posterior mean, median,
credible interval, draw counts, and a reliability flag. Adding
`weights = "bootstrap"` widens the target from model-parameter uncertainty alone
to model plus validation-sample uncertainty. `ppc_calibration()` returns a
posterior predictive p-value for calibration at the requested horizon.

## Model backends

Every backend returns the same model-agnostic `psurv` object: an `[S x n x H]`
posterior survival array plus the linear predictor. Downstream metrics are
therefore backend-agnostic.

| backend | needs | use when |
|---|---|---|
| `coxph` | `survival` | you want fast examples, teaching code, or a no-Stan fallback |
| `pwexp` | `survival` | you want a piecewise-exponential Poisson-trick PH model with a Laplace approximation |
| `stan_pwexp` | `rstan` | you want the primary full-likelihood Bayesian piecewise-exponential PH workflow |
| `bart` | `BART` | you want a nonparametric Bayesian survival comparator |
| `stan_surv` | `rstanarm` with `stan_surv` | aspirational/conditional; not available in the CRAN `rstanarm` release |

MCMC fits are gated by `check_mcmc_gate()` before any metric is reported. The
gate checks split-Rhat, effective sample size, and divergent transitions.
`prior_sensitivity()` refits under several log-hazard-ratio prior scales to
check robustness.

## Documentation

The pkgdown site is the main reference:

- package site: <https://wmotte.github.io/bayesval/>
- function index: <https://wmotte.github.io/bayesval/reference/>
- local entry point for writing a new metric: `?propagate_metric`
- full package help: `help(package = "bayesval")`

Posterior metric objects have `summary()` methods that report the median,
credible interval, finite-draw counts, and an unreliability flag when more than
5% of draws are non-finite.

Worked examples live in `vignettes/`, and most are external validations on
public survival cohorts:

- `worked_example.R` — Rotterdam → GBSG breast cancer (the headline example)
- `showcase_pbc.R` — Mayo PBC, randomized → non-randomized patients
- `temporal_validation.R` — one Rotterdam-trained model on a later Rotterdam era against GBSG
- `worked_example_lung.R` — NCCTG lung, split by recruiting centre
- `worked_example_colon.R` — colon adjuvant trial, split by treatment arm

`showcase_small_cohort.R` runs the deployment-gate showcase on a small simulated
cohort with a known oracle. Heavier Stan articles (`worked_example_stan.R`,
`small_sample_stan_vs_laplace.R`) are pre-rendered outside the normal pkgdown
build; see `articles/README.md` for that workflow.

## Development checks

For a lightweight local check without rebuilding vignettes:

```sh
R CMD build --no-build-vignettes --no-manual .
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --ignore-vignettes \
  bayesval_0.1.0.tar.gz
```

`_R_CHECK_FORCE_SUGGESTS_=false` is useful on machines without the heavier
Suggested packages. A complete CRAN-style check should install all Suggested
packages first.

## Citation and reproducibility

`bayesval` is the software companion to:

> Otte WM, Kampman SL, [...] *A Bayesian framework for validating survival
> prediction models: posterior credible intervals for discrimination,
> calibration, and clinical utility.* in prep.

## License

GPL-3. The full license text is in [LICENSE.md](LICENSE.md).
