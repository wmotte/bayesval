# metrics_calibration.R
# ---------------------------------------------------------------------------
# A2 -> B2. Per-draw calibration metrics. Calibration asks whether predicted
# risks pi_i(t) = 1 - S_i(t) NUMERICALLY match observed risk (level & shape) --
# distinct from discrimination. Each function has signature metric_fn(draw, w).
#
# The smooth "observed" calibration curve pi_hat^cal_i(t) is obtained via
# polspline::hare when available (the manuscript's primary estimator), else a
# grouped Kaplan-Meier fallback that needs only `survival`, so the pipeline runs
# out of the box.
# ---------------------------------------------------------------------------

#' Calibration-in-the-large at each horizon (per draw): mean predicted risk vs
#' observed (weighted Kaplan-Meier) risk.
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @return named vector CITL_t<h> = predicted_mean_risk - observed_risk.
#' @export
metric_citl <- function(draw, w = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  hz <- draw$horizons
  out <- numeric(length(hz))
  for (h in seq_along(hz)) {
    pred_risk <- stats::weighted.mean(1 - draw$surv2d[, h], w)
    obs_risk  <- 1 - .km_at(draw$time, draw$status, hz[h], w)
    out[h] <- pred_risk - obs_risk
  }
  names(out) <- paste0("CITL_t", seq_along(hz))
  out
}

#' Calibration slope (per draw): coefficient of a Cox model with eta as the
#' sole covariate. slope < 1 signals overfitting (van Houwelingen shrinkage).
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @return c(cal_slope = value)
#' @export
metric_cal_slope <- function(draw, w = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  df <- data.frame(time = draw$time, status = draw$status, eta = draw$eta)
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, status) ~ eta, data = df,
                    weights = w * draw$n),
    error = function(e) NULL)
  c(cal_slope = if (is.null(fit)) NA_real_ else unname(fit$coefficients[["eta"]]))
}

#' ICI / E50 / E90 at a horizon (per draw): mean / 50th / 90th percentile of
#' |pi_i(t) - pi_hat^cal_i(t)| (Austin, Harrell & van Klaveren 2020).
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @param horizon index into draw$horizons (default: last).
#' @param method  "auto" (hare if installed else group), "hare", or "group".
#' @return c(ICI = , E50 = , E90 = )
#' @export
metric_ici <- function(draw, w = NULL, horizon = NULL,
                       method = c("auto", "hare", "group")) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  method <- match.arg(method)
  h <- if (is.null(horizon)) draw$H else horizon
  t <- draw$horizons[h]
  pred_risk <- 1 - draw$surv2d[, h]
  obs_risk  <- .calibration_curve(pred_risk, draw$time, draw$status, t, w,
                                  method = method)
  ad <- abs(pred_risk - obs_risk)
  ok <- is.finite(ad)
  ww <- w[ok] / sum(w[ok]); ad <- ad[ok]
  c(ICI = sum(ww * ad),
    E50 = .wquantile(ad, ww, 0.50),
    E90 = .wquantile(ad, ww, 0.90))
}

#' D-calibration (per draw) — A2.6, Haider et al. (JMLR 2020). Distributional
#' calibration of the WHOLE predicted survival distribution, not one horizon: if
#' the model is calibrated, u_i = S_i(T_i) is Uniform(0,1). Bin `[0,1]` into
#' `bins` equal cells; an event contributes 1 to its cell, a censored subject
#' (known only that its true u < S_i(c_i)) spreads its unit mass uniformly over
#' `[0, S_i(c_i)]` across cells. The Pearson statistic vs the expected count n/bins
#' tests uniformity (chi-square, bins-1 df).
#'
#' Survival at each subject's own event/censoring time is interpolated from the
#' horizon-knot cumulative hazard, so D-calibration honours the full curve even
#' though the psurv array is stored at a few horizons.
#'
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @param bins number of equal-probability cells (default 10).
#' @return c(Dcal = chi-square statistic, Dcal_p = its p-value).
#' @export
metric_dcalibration <- function(draw, w = NULL, bins = 10L) {
  n <- draw$n
  if (is.null(w)) w <- rep(1 / n, n)
  cnt <- w * n                                   # weights on a "count" scale
  hz <- draw$horizons
  Lam <- .cumhaz_at_knots(draw$surv2d)
  ui <- .surv_at_own(Lam, c(0, hz), length(hz), draw$time)
  ui <- pmin(pmax(ui, 0), 1)

  edges <- seq(0, 1, length.out = bins + 1L)
  O <- numeric(bins)
  ev <- draw$status == 1
  # events: whole mass into the containing cell.
  if (any(ev)) {
    b <- pmin(pmax(findInterval(ui[ev], edges, rightmost.closed = TRUE), 1L),
              bins)
    O <- O + as.numeric(tapply(cnt[ev], factor(b, levels = seq_len(bins)),
                               sum, default = 0))
  }
  # censored: spread mass over [0, u_i] proportional to each cell's overlap.
  cen <- which(!ev)
  for (i in cen) {
    u <- ui[i]
    if (u <= 0) { O[1] <- O[1] + cnt[i]; next }
    ov <- pmax(0, pmin(edges[-1], u) - edges[-(bins + 1L)])
    O <- O + cnt[i] * ov / u
  }
  E <- sum(cnt) / bins
  stat <- sum((O - E)^2 / E)
  c(Dcal = unname(stat),
    Dcal_p = unname(stats::pchisq(stat, df = bins - 1L, lower.tail = FALSE)))
}

# ---------------------------------------------------------------------------
# Smooth "observed" calibration curve pi_hat^cal_i(t).
# ---------------------------------------------------------------------------

#' Per-subject calibrated observed risk at t as a function of predicted risk.
#' @keywords internal
.calibration_curve <- function(pred_risk, time, status, t, w,
                               method = "auto", groups = 10L) {
  use_hare <- method == "hare" ||
    (method == "auto" && requireNamespace("polspline", quietly = TRUE))
  if (use_hare) {
    out <- tryCatch(.calib_hare(pred_risk, time, status, t),
                    error = function(e) NULL)
    if (!is.null(out)) return(out)
    # fall through to grouped if hare failed numerically.
  }
  .calib_group(pred_risk, time, status, t, w, groups)
}

# hare: flexible hazard regression S(t | predictor). Predictor on cloglog scale.
#' @keywords internal
.calib_hare <- function(pred_risk, time, status, t) {
  pr <- pmin(pmax(pred_risk, 1e-6), 1 - 1e-6)
  cll <- log(-log(1 - pr))                 # complementary log-log of risk
  fit <- polspline::hare(data = time, delta = status, cov = as.matrix(cll))
  surv_t <- 1 - polspline::phare(q = t, cov = as.matrix(cll), fit = fit)
  1 - surv_t                                # observed risk at t
}

# Grouped Kaplan-Meier fallback (survival only): bin by predicted risk, take the
# weighted KM observed risk per bin, interpolate back to each subject.
#' @keywords internal
.calib_group <- function(pred_risk, time, status, t, w, groups = 10L) {
  n <- length(pred_risk)
  g <- min(groups, max(2L, floor(n / 20L)))
  br <- stats::quantile(pred_risk, probs = seq(0, 1, length.out = g + 1),
                        na.rm = TRUE)
  br[1] <- -Inf; br[length(br)] <- Inf
  bin <- cut(pred_risk, breaks = unique(br), include.lowest = TRUE,
             labels = FALSE)
  xg <- yg <- rep(NA_real_, max(bin, na.rm = TRUE))
  for (b in sort(unique(bin))) {
    ix <- which(bin == b)
    xg[b] <- stats::weighted.mean(pred_risk[ix], w[ix])
    yg[b] <- 1 - .km_at(time[ix], status[ix], t, w[ix])
  }
  ok <- is.finite(xg) & is.finite(yg)
  if (sum(ok) < 2) return(rep(stats::weighted.mean(1 -
      .km_at(time, status, t, w), 1), n))
  stats::approx(x = xg[ok], y = yg[ok], xout = pred_risk,
                method = "linear", rule = 2, ties = mean)$y
}

# --- helpers ---------------------------------------------------------------

#' Weighted Kaplan-Meier survival at a single time t.
#' @keywords internal
.km_at <- function(time, status, t, w = NULL) {
  if (is.null(w)) w <- rep(1, length(time))
  fit <- survival::survfit(survival::Surv(time, status) ~ 1,
                           weights = w * length(time))
  stats::approx(x = c(0, fit$time), y = c(1, fit$surv), xout = t,
                method = "constant", rule = 2, f = 0)$y
}

#' Weighted quantile (type-7-ish, linear on the weighted ECDF).
#' @keywords internal
.wquantile <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) - 0.5 * w
  cw <- cw / sum(w)
  stats::approx(x = cw, y = x, xout = p, rule = 2)$y
}
