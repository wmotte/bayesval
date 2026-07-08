# Reproducibility (checklist 8): with fixed seeds the whole pipeline must be
# deterministic. These are regression snapshots — if a change moves a posterior
# median, the test fails and the change must be justified (and the snapshot
# updated deliberately).
library(testthat)
suppressMessages(library(survival))

fixed_ps <- function() {
  set.seed(2024); n <- 400
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.15 * exp(eta)); C <- rexp(n, 0.10)
  time <- pmin(T, C); status <- as.integer(T <= C)
  dat <- data.frame(time, status, x1, x2)
  hz <- as.numeric(quantile(time[status == 1], c(0.25, 0.5, 0.75)))
  list(ps = fit_bayes_survival(Surv(time, status) ~ x1 + x2, data = dat,
                               horizons = hz, backend = "coxph", S = 200L,
                               seed = 99L), hz = hz)
}

med <- function(pm, name = NULL) {
  s <- summary(pm)
  if (is.null(name)) s$median else s$median[s$metric == name]
}

test_that("posterior medians match stored snapshots (fixed seed)", {
  f <- fixed_ps(); ps <- f$ps; hz <- f$hz
  tol <- 1e-5
  expect_equal(med(propagate_metric(ps, metric_harrell_c, seed = 5)),
               0.708853, tolerance = tol)
  expect_equal(med(propagate_metric(ps, function(d, w)
    metric_uno_c(d, w, tau = max(hz)), seed = 5)), 0.701676, tolerance = tol)
  expect_equal(med(propagate_metric(ps, metric_gonen_k, seed = 5)),
               0.694396, tolerance = tol)
  expect_equal(med(propagate_metric(ps, metric_antolini_c, seed = 5)),
               0.708853, tolerance = tol)
  expect_equal(med(propagate_metric(ps, metric_brier, seed = 5), "IBS"),
               0.172147, tolerance = tol)
  expect_equal(med(propagate_metric(ps, metric_r2d, seed = 5)),
               0.281126, tolerance = tol)
  expect_equal(med(propagate_metric(ps, metric_dcalibration, seed = 5), "Dcal"),
               2.277165, tolerance = 1e-4)
})

test_that("Bayesian bootstrap is reproducible and wider than parameter-only", {
  f <- fixed_ps(); ps <- f$ps
  expect_equal(med(propagate_metric(ps, metric_harrell_c, weights = "bootstrap",
                                    seed = 5)), 0.70751, tolerance = 1e-4)
  wn <- diff(unlist(summary(propagate_metric(ps, metric_harrell_c,
    weights = "none", seed = 5))[c("lower", "upper")]))
  wb <- diff(unlist(summary(propagate_metric(ps, metric_harrell_c,
    weights = "bootstrap", seed = 5))[c("lower", "upper")]))
  expect_gt(wb, wn)
})

test_that("re-fitting with the same seed gives identical draws", {
  a <- fixed_ps()$ps; b <- fixed_ps()$ps
  expect_identical(a$surv, b$surv)
  expect_identical(a$eta, b$eta)
})
