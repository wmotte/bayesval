# Part C machinery checks: DGM correctness, oracle sanity, and that the
# coverage evaluator computes coverage/width/bias correctly on known input.
library(testthat)
suppressMessages(library(survival))
source(file.path("..", "..", "R", "load_all.R"), chdir = TRUE)
source(file.path("..", "..", "simulation", "dgm.R"), chdir = TRUE)
source(file.path("..", "..", "simulation", "evaluate_coverage.R"), chdir = TRUE)

test_that("censoring rate is calibrated to the target", {
  set.seed(1)
  for (cr in c(0.2, 0.5, 0.7)) {
    d <- dgm_weibull(20000, cens_rate = cr)
    expect_equal(mean(d$data$status == 0), cr, tolerance = 0.03)
  }
})

test_that("coxph recovers the true betas from the Weibull DGM", {
  set.seed(2)
  d <- dgm_weibull(20000, beta = c(b1 = 0.8, b2 = -0.5), cens_rate = 0.3)
  fit <- coxph(Surv(time, status) ~ x1 + x2, data = d$data)
  expect_equal(unname(coef(fit)), c(0.8, -0.5), tolerance = 0.06)
})

test_that("piecewise-exponential DGM censoring calibrates and betas recover", {
  set.seed(3)
  d <- dgm_piecewise(20000, cens_rate = 0.4)
  expect_equal(mean(d$data$status == 0), 0.4, tolerance = 0.03)
  fit <- coxph(Surv(time, status) ~ x1 + x2, data = d$data)
  expect_equal(unname(coef(fit)), c(0.8, -0.5), tolerance = 0.08)
})

test_that("oracle C lies in (0.5, 1) and IBS in (0, 0.25)", {
  hz <- c(1, 3, 5)
  orc <- oracle_metrics(function(n) dgm_weibull(n, cens_rate = 0.3),
                        horizons = hz, N = 3e4)
  expect_gt(orc$C, 0.5); expect_lt(orc$C, 1)
  expect_gt(orc$IBS, 0);  expect_lt(orc$IBS, 0.25)
  expect_equal(orc$cal_slope, 1)
})

test_that("evaluate_coverage computes coverage/width/bias correctly", {
  # Construct estimates with known coverage: truth = 0.7, half cover it.
  res <- data.frame(
    mode = "eval", metric = "C_uno", method = "bayes_bootstrap",
    est = 0.70, lower = c(0.65, 0.72), upper = c(0.75, 0.80),
    truth = 0.70, n = 300L, cens = 0.3, misspec = "none",
    rep = 1:2, stringsAsFactors = FALSE)
  expect_no_warning(ev <- evaluate_coverage(res))
  expect_equal(ev$coverage, 0.5)             # first covers, second does not
  expect_equal(ev$width, mean(c(0.10, 0.08)))
  expect_equal(ev$bias, 0)
  expect_equal(unique(ev$sep), "med")
})

test_that("bb_replicates returns a B x k matrix that varies with weights", {
  set.seed(4)
  d <- dgm_weibull(300, cens_rate = 0.3)
  hz <- as.numeric(quantile(d$data$time[d$data$status == 1], c(.25, .5, .75)))
  lp <- 0.8 * d$data$x1 - 0.5 * d$data$x2
  surv2d <- exp(-outer(exp(lp), c(0.1, 0.3, 0.5)))
  draw <- list(surv2d = surv2d, eta = lp, time = d$data$time,
               status = d$data$status, horizons = hz, n = 300L)
  bb <- bb_replicates(draw, function(dd, w) metric_uno_c(dd, w, tau = max(hz)),
                      B = 50L)
  expect_equal(dim(bb), c(50L, 1L))
  expect_gt(stats::sd(bb[, 1]), 0)           # weights genuinely perturb it
})
