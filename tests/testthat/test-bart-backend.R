# Nonparametric BART-survival backend (checklist 1): the framework is
# model-agnostic, so a Bayesian machine-learning survival model must flow through
# the same `psurv` interface and yield sensible posteriors. Skipped when BART is
# unavailable.
library(testthat)
suppressMessages(library(survival))

test_that("BART backend yields a valid psurv and coherent metrics", {
  skip_if_not_installed("BART")
  set.seed(6); n <- 250
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.15 * exp(eta)); C <- rexp(n, 0.10)
  time <- pmin(T, C); status <- as.integer(T <= C)
  dat <- data.frame(time, status, x1, x2)
  hz <- as.numeric(quantile(time[status == 1], c(0.25, 0.5, 0.75)))

  ps <- fit_bayes_survival(Surv(time, status) ~ x1 + x2, dat, hz,
                           backend = "bart", time = "time", status = "status",
                           S = 150L, seed = 1L)
  expect_s3_class(ps, "psurv")
  expect_equal(dim(ps$surv), c(nrow(ps$eta), n, length(hz)))
  expect_true(all(ps$surv >= 0 & ps$surv <= 1, na.rm = TRUE))

  har <- summary(propagate_metric(ps, metric_harrell_c))
  ibs <- summary(propagate_metric(ps, metric_brier))
  expect_true(is.finite(har$median) && har$median > 0.6 && har$median < 0.9)
  expect_true(all(is.finite(ibs$median)))

  # nonparametric BART should land near the Cox fit on PH-generated data.
  psc <- fit_bayes_survival(Surv(time, status) ~ x1 + x2, dat, hz,
                            backend = "coxph", S = 150L, seed = 1L)
  cox_har <- summary(propagate_metric(psc, metric_harrell_c))$median
  expect_lt(abs(har$median - cox_har), 0.05)
})
