# Edge-case robustness (checklist 8): heavy censoring, tied event times, a
# subgroup with zero events, and a truncation tau beyond the last event. Metrics
# should return finite values (or a clean NA) without erroring.
library(testthat)
suppressMessages(library(survival))

make_ps <- function(dat, horizons, S = 80L, seed = 1L)
  fit_bayes_survival(Surv(time, status) ~ x1 + x2, data = dat,
                     horizons = horizons, backend = "coxph", S = S, seed = seed)

test_that("heavy censoring (85%) still yields finite discrimination + IBS", {
  set.seed(31); n <- 400
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.15 * exp(eta)); C <- rexp(n, 1.2)   # heavy censoring
  time <- pmin(T, C); status <- as.integer(T <= C)
  dat <- data.frame(time, status, x1, x2)
  expect_gt(mean(status == 0), 0.7)
  hz <- as.numeric(quantile(time[status == 1], c(0.25, 0.5, 0.75)))
  ps <- make_ps(dat, hz)
  s_uno <- summary(propagate_metric(ps, function(d, w)
    metric_uno_c(d, w, tau = max(hz))))
  s_ibs <- summary(propagate_metric(ps, metric_brier))
  expect_true(is.finite(s_uno$median))
  expect_true(all(is.finite(s_ibs$median)))
  expect_true(s_uno$median > 0.4 && s_uno$median < 1)
})

test_that("many tied event times do not break concordance or calibration", {
  set.seed(32); n <- 300
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.2 * exp(eta))
  C <- rexp(n, 0.1)
  time <- ceiling(pmin(T, C))                        # integer times => heavy ties
  status <- as.integer(T <= C)
  dat <- data.frame(time = pmax(time, 1), status, x1, x2)
  expect_true(max(table(dat$time)) > 10)             # genuine ties
  hz <- as.numeric(quantile(dat$time[dat$status == 1], c(0.3, 0.6, 0.9)))
  ps <- make_ps(dat, hz)
  har <- summary(propagate_metric(ps, metric_harrell_c))
  ant <- summary(propagate_metric(ps, metric_antolini_c))
  ici <- summary(propagate_metric(ps, function(d, w)
    metric_ici(d, w, method = "group")))
  expect_true(is.finite(har$median) && is.finite(ant$median))
  expect_true(all(is.finite(ici$median)))
})

test_that("a subgroup with zero events is handled by calibration + net benefit", {
  set.seed(33); n <- 300
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.15 * exp(eta)); C <- rexp(n, 0.1)
  time <- pmin(T, C); status <- as.integer(T <= C)
  # Force the lowest-risk decile to be event-free (censor them early).
  lo <- order(eta)[1:30]; status[lo] <- 0L; time[lo] <- pmin(time[lo], 0.05)
  dat <- data.frame(time, status, x1, x2)
  hz <- as.numeric(quantile(time[status == 1], c(0.25, 0.5, 0.75)))
  ps <- make_ps(dat, hz)
  ici <- summary(propagate_metric(ps, function(d, w)
    metric_ici(d, w, method = "group")))
  pnb <- posterior_net_benefit(ps, thresholds = seq(0.1, 0.5, 0.1),
                               draws = 50L)
  expect_true(all(is.finite(ici$median)))
  expect_true(all(is.finite(pnb$model$median)))
})

test_that("tau beyond the last event time does not error", {
  set.seed(34); n <- 250
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.2 * exp(eta)); C <- rexp(n, 0.15)
  time <- pmin(T, C); status <- as.integer(T <= C)
  dat <- data.frame(time, status, x1, x2)
  tau_big <- max(time) * 1.5                          # beyond last observation
  hz <- c(as.numeric(quantile(time[status == 1], c(0.5))), tau_big)
  ps <- make_ps(dat, hz)
  uno <- summary(propagate_metric(ps, function(d, w)
    metric_uno_c(d, w, tau = tau_big)))
  dcal <- summary(propagate_metric(ps, metric_dcalibration))
  expect_true(is.finite(uno$median))
  expect_true(all(is.finite(dcal$median[dcal$metric == "Dcal"])))
})

test_that("Gonen-Heller subsampling guard keeps K in (0,1) for large n", {
  set.seed(35); n <- 200
  eta <- rnorm(n)
  d <- list(eta = eta, n = n)
  k <- metric_gonen_k(d, max_n = 100L)[["K_gonen"]]   # force subsample
  expect_true(k > 0 && k < 1)
})
