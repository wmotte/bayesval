# Full-Bayes backend validation (checklist 1): the psurv array from the genuine
# rstan piecewise-exponential PH model must match an independent reconstruction
# from the raw posterior draws (no indexing bug), must agree with a frequentist
# Cox prediction on the same data, and must pass through the MCMC-diagnostics
# gate. Skipped when rstan is unavailable (the survival-only backends carry the
# rest of the suite).
library(testthat)
suppressMessages(library(survival))

test_that("stan_pwexp psurv array matches a raw-draw reconstruction", {
  skip_if_not_installed("rstan")
  wd <- setwd(normalizePath(file.path("..", ".."))); on.exit(setwd(wd), add = TRUE)
  set.seed(7); n <- 250
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.15 * exp(eta)); C <- rexp(n, 0.10)
  time <- pmin(T, C); status <- as.integer(T <= C)
  dat <- data.frame(time, status, x1, x2)
  hz <- as.numeric(quantile(time[status == 1], c(0.25, 0.5, 0.75)))

  obj <- fit_stan_pwexp(Surv(time, status) ~ x1 + x2, dat, hz, n_pieces = 6L,
                        chains = 2L, iter = 700L, seed = 1L)
  ps <- obj$psurv
  X <- model.matrix(~ x1 + x2, dat)[, -1, drop = FALSE]
  edges <- obj$edges

  # Independent reconstruction for a few (draw, subject, horizon) cells.
  set.seed(123)
  for (rep in 1:15) {
    s <- sample.int(nrow(obj$beta_draws), 1)
    i <- sample.int(n, 1); h <- sample.int(length(hz), 1)
    width <- pmax(0, pmin(hz[h], edges[-1]) - edges[-length(edges)])
    H0 <- sum(width * exp(obj$log_h0_draws[s, ]))
    Sref <- exp(-H0 * exp(sum(X[i, ] * obj$beta_draws[s, ])))
    expect_equal(ps$surv[s, i, h], Sref, tolerance = 1e-9)
  }
})

test_that("stan_pwexp survival agrees with a Cox prediction (both PH)", {
  skip_if_not_installed("rstan")
  wd <- setwd(normalizePath(file.path("..", ".."))); on.exit(setwd(wd), add = TRUE)
  set.seed(8); n <- 400
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.15 * exp(eta)); C <- rexp(n, 0.10)
  time <- pmin(T, C); status <- as.integer(T <= C)
  dat <- data.frame(time, status, x1, x2)
  hz <- as.numeric(quantile(time[status == 1], c(0.5)))

  ps <- fit_bayes_survival(Surv(time, status) ~ x1 + x2, dat, hz,
                           backend = "stan_pwexp", n_pieces = 8L,
                           chains = 2L, iter = 900L, seed = 2L)
  bayes_S <- colMeans(ps$surv[, , 1])                 # posterior-mean S_i(t)

  cx <- coxph(Surv(time, status) ~ x1 + x2, dat, x = TRUE)
  bh <- basehaz(cx, centered = FALSE)
  H0 <- approx(bh$time, bh$hazard, xout = hz, method = "constant",
               rule = 2, f = 0)$y
  cox_S <- exp(-H0 * exp(as.numeric(cx$x %*% cx$coefficients)))

  # Two PH fits (full-Bayes vs partial-likelihood) should predict closely.
  expect_lt(mean(abs(bayes_S - cox_S)), 0.03)
  expect_gt(cor(bayes_S, cox_S), 0.99)
})

test_that("MCMC diagnostics gate runs on a real stanfit", {
  skip_if_not_installed("rstan")
  wd <- setwd(normalizePath(file.path("..", ".."))); on.exit(setwd(wd), add = TRUE)
  set.seed(9); n <- 250
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  T <- rexp(n, 0.15 * exp(0.8 * x1 - 0.5 * x2)); C <- rexp(n, 0.10)
  time <- pmin(T, C); status <- as.integer(T <= C)
  dat <- data.frame(time, status, x1, x2)
  hz <- as.numeric(quantile(time[status == 1], c(0.25, 0.5, 0.75)))

  ps <- fit_bayes_survival(Surv(time, status) ~ x1 + x2, dat, hz,
                           backend = "stan_pwexp", n_pieces = 6L,
                           chains = 2L, iter = 1200L, seed = 3L)
  d <- mcmc_diagnostics(ps)
  expect_true(d$applicable)
  expect_true(is.finite(d$max_rhat))
  expect_true(is.finite(d$min_ess_bulk))
  # divergence count now resolves for a raw/attached stanfit (not NA).
  expect_false(is.na(d$n_divergent))
})
