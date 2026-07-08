# Verification (plan section "Verification" + checklist 8): each per-draw metric
# must match its reference frequentist estimator on fixed data at a fixed
# linear predictor, to numerical tolerance.
library(testthat)
suppressMessages(library(survival))

make_draw <- function(seed = 1L, n = 400L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  eta <- 0.8 * x1 - 0.5 * x2
  T <- rexp(n, 0.15 * exp(eta)); C <- rexp(n, 0.10)
  time <- pmin(T, C); status <- as.integer(T <= C)
  hz <- as.numeric(quantile(time[status == 1], c(0.25, 0.5, 0.75)))
  fit <- coxph(Surv(time, status) ~ x1 + x2)
  bh <- basehaz(fit, centered = FALSE)
  H0 <- approx(bh$time, bh$hazard, xout = hz, method = "constant",
               rule = 2, f = 0)$y
  surv2d <- exp(-outer(exp(eta), H0))
  list(time = time, status = status, eta = eta, n = n, horizons = hz,
       surv2d = surv2d)
}

test_that("KM helpers return survival 1 before the first observed time", {
  expect_equal(.km_at(c(5, 6, 7), c(1, 0, 0), t = 1), 1)

  Ghat <- .km_censoring(c(5, 6, 7), c(0, 1, 1))
  expect_equal(.step_eval(Ghat, 1), 1)
})

test_that("early-horizon IPCW Brier uses unit censoring survival before first censoring", {
  time <- c(5, 6, 7, 8)
  status <- c(0, 1, 1, 1)
  S_t <- c(0.8, 0.9, 0.7, 0.6)
  w <- rep(1 / length(time), length(time))
  Ghat <- .km_censoring(time, status)

  got <- .brier_ipcw(time, status, S_t, t = 1, Ghat = Ghat, w = w)
  expect_equal(got, mean((1 - S_t)^2), tolerance = 1e-12)
})

test_that("Harrell C matches survival::concordance (risk orientation)", {
  d <- make_draw(1L)
  ref <- 1 - concordance(Surv(d$time, d$status) ~ d$eta)$concordance
  expect_equal(unname(metric_harrell_c(d)[["C_harrell"]]), ref, tolerance = 1e-8)
})

test_that("Uno C matches survival::concordance timewt = n/G2", {
  d <- make_draw(2L)
  tau <- max(d$horizons)
  ref <- 1 - concordance(Surv(d$time, d$status) ~ d$eta,
                         timewt = "n/G2", ymax = tau)$concordance
  expect_equal(unname(metric_uno_c(d, tau = tau)[["C_uno"]]), ref,
               tolerance = 1e-8)
})

test_that("calibration slope matches a Cox refit on eta", {
  d <- make_draw(3L)
  ref <- coxph(Surv(d$time, d$status) ~ d$eta)$coefficients[[1]]
  expect_equal(unname(metric_cal_slope(d)[["cal_slope"]]), ref, tolerance = 1e-6)
})

test_that("Gonen-Heller K reproduces its closed form and lies in (0,1)", {
  d <- make_draw(4L, n = 200L)
  eta <- d$eta; n <- length(eta)
  dmat <- outer(eta, eta, `-`)
  W <- matrix(1, n, n); diag(W) <- 0
  ref <- 2 * sum(W * ((dmat > 0) * plogis(dmat))) / sum(W)
  got <- metric_gonen_k(d)[["K_gonen"]]
  expect_equal(unname(got), ref, tolerance = 1e-8)
  expect_gt(got, 0.5)  # eta is genuinely predictive here
})

test_that("Brier at t reproduces a hand IPCW computation", {
  d <- make_draw(5L)
  h <- 2L; t <- d$horizons[h]; S_t <- d$surv2d[, h]
  Gfit <- survfit(Surv(d$time, 1 - d$status) ~ 1)
  Gt   <- approx(Gfit$time, Gfit$surv, xout = t, method = "constant",
                 rule = 2, f = 0)$y
  Gti  <- approx(Gfit$time, Gfit$surv, xout = d$time, method = "constant",
                 rule = 2, f = 0)$y
  ev <- (d$time <= t) & (d$status == 1); sp <- d$time > t
  w  <- numeric(d$n)
  w[ev] <- 1 / pmax(Gti[ev], 1e-8); w[sp] <- 1 / max(Gt, 1e-8)
  r2 <- numeric(d$n); r2[ev] <- (0 - S_t[ev])^2; r2[sp] <- (1 - S_t[sp])^2
  ref <- mean(w * r2)
  got <- metric_brier(d)[[paste0("Brier_t", h)]]
  expect_equal(unname(got), ref, tolerance = 1e-8)
})

test_that("Antolini C equals Harrell C under proportional hazards", {
  d <- make_draw(8L, n = 300L)
  har <- metric_harrell_c(d)[["C_harrell"]]
  ant <- metric_antolini_c(d)[["C_antolini"]]
  # Under PH the risk ordering is time-invariant, so the time-dependent
  # concordance must reduce to Harrell's (up to knot-interpolation ties).
  expect_equal(unname(ant), unname(har), tolerance = 0.01)
})

test_that("D-calibration reproduces its binning statistic by hand", {
  set.seed(21); n <- 400L
  # One horizon, all times == that horizon, all events => u_i = S_i(t) = surv2d.
  ui <- runif(n)
  d <- list(time = rep(2, n), status = rep(1L, n), n = n, horizons = 2,
            surv2d = matrix(ui, ncol = 1L),
            eta = rnorm(n))
  bins <- 10L
  edges <- seq(0, 1, length.out = bins + 1L)
  b <- pmin(pmax(findInterval(ui, edges, rightmost.closed = TRUE), 1L), bins)
  O <- as.numeric(table(factor(b, levels = seq_len(bins))))
  E <- n / bins
  ref_stat <- sum((O - E)^2 / E)
  got <- metric_dcalibration(d, bins = bins)
  expect_equal(unname(got[["Dcal"]]), ref_stat, tolerance = 1e-8)
  expect_equal(unname(got[["Dcal_p"]]),
               pchisq(ref_stat, df = bins - 1L, lower.tail = FALSE),
               tolerance = 1e-10)
})

test_that("D-calibration rejects a deliberately miscalibrated model", {
  d <- make_draw(9L, n = 500L)
  good <- metric_dcalibration(d)[["Dcal_p"]]
  d_bad <- d; d_bad$surv2d <- d$surv2d^3            # risks pushed too high
  bad  <- metric_dcalibration(d_bad)[["Dcal_p"]]
  expect_gt(good, 0.05)     # correct predictions: not rejected
  expect_lt(bad, 0.05)      # miscalibrated: rejected
})

test_that("net benefit at a threshold matches a hand computation", {
  d <- make_draw(10L, n = 400L)
  h <- 3L; t <- d$horizons[h]; pt <- 0.3
  d2 <- c(d, list(surv2d = d$surv2d))               # net_benefit_draw needs $H
  d2$H <- length(d$horizons)
  nb <- net_benefit_draw(d2, thresholds = pt, horizon = h)$model
  # hand: f = P(flagged), r = 1 - KM among flagged at t.
  risk <- 1 - d$surv2d[, h]; pos <- risk >= pt
  f <- mean(pos)
  km <- survfit(Surv(d$time[pos], d$status[pos]) ~ 1)
  r <- 1 - approx(km$time, km$surv, xout = t, method = "constant",
                  rule = 2, f = 0)$y
  ref <- f * r - f * (1 - r) * pt / (1 - pt)
  expect_equal(unname(nb), ref, tolerance = 1e-8)
})

test_that("IPA of the null model is 0 and R2_D matches its formula", {
  d <- make_draw(11L, n = 300L)
  # A model whose predictions ARE the marginal KM has IPA == 0 at every horizon.
  km_t <- vapply(d$horizons, function(t)
    .km_at(d$time, d$status, t), numeric(1))
  d_null <- d
  d_null$surv2d <- matrix(km_t, nrow = d$n, ncol = length(d$horizons),
                          byrow = TRUE)
  ipa <- metric_ipa(d_null)
  expect_true(all(abs(ipa) < 1e-6))

  z <- qnorm((rank(d$eta) - 0.375) / (d$n + 0.25)) / sqrt(8 / pi)
  D <- coxph(Surv(d$time, d$status) ~ z)$coefficients[[1]]
  ref <- (D^2 / (8 / pi)) / (pi^2 / 6 + D^2 / (8 / pi))
  expect_equal(unname(metric_r2d(d)[["R2_D"]]), ref, tolerance = 1e-6)
})

test_that("propagate_metric returns one row per used draw and summarises", {
  d <- make_draw(6L, n = 150L)
  S <- 40L
  surv <- array(rep(d$surv2d, each = S), dim = c(S, d$n, length(d$horizons)))
  eta  <- matrix(d$eta, S, d$n, byrow = TRUE)
  ps <- new_psurv(surv, eta, d$time, d$status, d$horizons)
  pm <- propagate_metric(ps, metric_harrell_c, weights = "none")
  expect_equal(nrow(unclass(pm)), S)
  s <- summary(pm)
  expect_true(all(c("median", "lower", "upper") %in% names(s)))
})

test_that("Bayesian bootstrap widens the interval vs parameter-only (B3)", {
  d <- make_draw(7L, n = 300L)
  set.seed(11)
  S <- 200L
  # Parameter spread: jitter eta across draws; survival held at its draw-1 value
  # so the ONLY extra width under bootstrap comes from the Dirichlet weights.
  etaS <- sweep(matrix(rnorm(S * d$n, 0, 0.05), S, d$n), 2, d$eta, `+`)
  surv <- array(rep(d$surv2d, each = S), dim = c(S, d$n, length(d$horizons)))
  ps <- new_psurv(surv, etaS, d$time, d$status, d$horizons)
  p_none <- summary(propagate_metric(ps, metric_harrell_c, weights = "none"))
  p_boot <- summary(propagate_metric(ps, metric_harrell_c, weights = "bootstrap"))
  expect_gt(p_boot$upper - p_boot$lower, p_none$upper - p_none$lower)
})
