# small_sample_stan_vs_laplace.R
# ---------------------------------------------------------------------------
# Reviewer-requested robustness check (Epidemiology R1 #11, R2 #15): the
# Rotterdam-to-GBSG agreement between the full-Bayes MCMC model and the
# Cox-Laplace approximation (worked_example_stan.R) is a large, well-identified
# example. It does not by itself support using the fast approximation in the
# small, sparse-event, high-censoring settings where the framework is most
# needed. This script fits BOTH backends on one small, high-censoring training
# cohort (n = 150, ~60% censored, moderate signal) and validates both on an
# independent small cohort of the same size, so any divergence between the
# genuine posterior and its large-sample approximation shows up directly in
# the propagated metrics, and MCMC pathologies (divergences, poor mixing) that
# a large well-identified fit would not reveal have a chance to appear.
#
#   Rscript vignettes/small_sample_stan_vs_laplace.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ---------------------------------------------------------------------------
suppressMessages({library(survival); library(rstan)})
source(file.path("..", "R", "load_all.R"), chdir = TRUE)
source(file.path("..", "simulation", "dgm.R"), chdir = TRUE)
rstan::rstan_options(auto_write = TRUE)
figdir <- "figures"; dir.create(figdir, showWarnings = FALSE)
set.seed(2027)

## ---- one small, sparse-event, high-censoring cohort pair ------------------
N_TRAIN <- 150L; N_VAL <- 150L; CENS <- 0.60
BETA <- c(b1 = 0.70, b2 = -0.50); SHAPE <- 1.3; RATE <- 0.14
horizons <- c(2, 5); h5 <- 2L; tau <- 5

train <- dgm_weibull(N_TRAIN, beta = BETA, shape = SHAPE, rate = RATE,
                     cens_rate = CENS)$data
val   <- dgm_weibull(N_VAL,   beta = BETA, shape = SHAPE, rate = RATE,
                     cens_rate = CENS)$data
cat(sprintf("Training cohort: n=%d, %d events (%.0f%% censored)\n",
           nrow(train), sum(train$status), 100 * mean(train$status == 0)))
cat(sprintf("Validation cohort: n=%d, %d events (%.0f%% censored)\n",
           nrow(val), sum(val$status), 100 * mean(val$status == 0)))

form <- Surv(time, status) ~ x1 + x2

## ---- Cox-Laplace backend ---------------------------------------------------
fit_cox <- coxph(form, data = train, x = TRUE, y = TRUE)
ps_cox  <- psurv_external(fit_cox, val, time = "time", status = "status",
                          horizons = horizons, S = 1000L, seed = 1L)

## ---- Genuine full-Bayes MCMC (fewer pieces: sparse events) -----------------
cat("Fitting full-Bayes piecewise-exponential PH (rstan) on the small cohort...\n")
t0 <- Sys.time()
obj <- fit_stan_pwexp(form, train, horizons = horizons, n_pieces = 4L,
                      beta_scale = 2.5, chains = 4L, iter = 4000L,
                      seed = 2027, adapt_delta = 0.995)
cat(sprintf("fitted in %.0f s\n", difftime(Sys.time(), t0, units = "secs")))
print(obj)
ps_stan <- stan_pwexp_predict(obj, val, time = "time", status = "status")

## ---- propagate the headline suite on both backends -------------------------
suite <- list(
  harrell = metric_harrell_c,
  uno     = function(d, w) metric_uno_c(d, w, tau = tau),
  slope   = metric_cal_slope,
  brier   = function(d, w) metric_brier(d, w)[c(paste0("Brier_t", h5), "IBS")])

post_cox  <- propagate_many(ps_cox,  suite, weights = "bootstrap", seed = 1)
post_stan <- propagate_many(ps_stan, suite, weights = "bootstrap", seed = 1)
tab_cox  <- posterior_table(post_cox)  |> dplyr::mutate(metric = sub("^[^.]+\\.", "", metric))
tab_stan <- posterior_table(post_stan) |> dplyr::mutate(metric = sub("^[^.]+\\.", "", metric))

fmt <- function(x) sprintf("%.3f [%.3f, %.3f]", x$median, x$lower, x$upper)
pick <- function(tab, nm) dplyr::filter(tab, metric == nm)
cmp <- tibble::tibble(
  metric     = c("C_harrell", "C_uno", "cal_slope", "IBS"),
  full_bayes = c(fmt(pick(tab_stan, "C_harrell")), fmt(pick(tab_stan, "C_uno")),
                fmt(pick(tab_stan, "cal_slope")), fmt(pick(tab_stan, "IBS"))),
  cox_laplace = c(fmt(pick(tab_cox, "C_harrell")), fmt(pick(tab_cox, "C_uno")),
                 fmt(pick(tab_cox, "cal_slope")), fmt(pick(tab_cox, "IBS"))))
cat("\n== Small-sample (n=150, 60% censored) full-Bayes MCMC vs Cox-Laplace ==\n")
print(as.data.frame(cmp), row.names = FALSE)
cat(sprintf("\nMCMC diagnostics on the small cohort: max Rhat %.3f, min ESS %.0f, divergences %d\n",
           obj$diagnostics$max_rhat, obj$diagnostics$min_ess, obj$diagnostics$divergences))

saveRDS(list(cmp = cmp, diagnostics = obj$diagnostics,
            train_n = nrow(train), train_events = sum(train$status),
            val_n = nrow(val), val_events = sum(val$status)),
       file.path(figdir, "small_sample_stan_vs_laplace_results.rds"))
cat("Done. Results in vignettes/figures/small_sample_stan_vs_laplace_results.rds\n")
