# metrics_discrimination.R
# ---------------------------------------------------------------------------
# A1 -> B2. Per-draw discrimination metrics. Each is a function
#   metric_fn(draw, w) -> named numeric
# suitable for propagate_metric(). `draw` is a psurv_draw() list; `w` is a
# length-n weight vector (Dirichlet under the Bayesian bootstrap, uniform
# otherwise). Estimators are REUSED from `survival` (Harrell/Uno) or computed
# from the closed form (Gonen-Heller); the novelty is only the posterior loop.
# ---------------------------------------------------------------------------

#' Harrell's C (per draw). Wraps survival::concordance with case weights.
#'
#' survival::concordance(Surv ~ marker) uses the convention that a *larger*
#' marker predicts *longer* survival. eta is a risk score (larger = shorter
#' survival), so we feed -eta; the returned concordance is then Harrell's C for
#' the risk score (equivalently 1 - concordance(Surv ~ eta)).
#'
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @return c(C_harrell = value)
#' @export
metric_harrell_c <- function(draw, w = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  df <- data.frame(time = draw$time, status = draw$status, marker = -draw$eta)
  cc <- survival::concordance(survival::Surv(time, status) ~ marker,
                              data = df, weights = w * draw$n)
  c(C_harrell = unname(cc$concordance))
}

#' Uno's C (IPCW, truncated at tau) per draw. survival::concordance with
#' timewt = "n/G2" is the Uno estimator; ymax = tau truncates.
#'
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @param tau truncation time; defaults to the largest horizon.
#' @return c(C_uno = value)
#' @export
metric_uno_c <- function(draw, w = NULL, tau = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  if (is.null(tau)) tau <- max(draw$horizons)
  df <- data.frame(time = draw$time, status = draw$status, marker = -draw$eta)
  cc <- survival::concordance(survival::Surv(time, status) ~ marker,
                              data = df, weights = w * draw$n,
                              timewt = "n/G2", ymax = tau)
  c(C_uno = unname(cc$concordance))
}

#' Gonen-Heller K (per draw) — model-based concordance from eta only, so it is
#' independent of the study censoring distribution. Closed form (A1):
#'   `K = sum_{i!=j} w_i w_j I(eta_i > eta_j) / (1 + exp(eta_j - eta_i))`
#' normalised by the pair weight mass. O(n^2); intended for moderate n.
#'
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @param max_n guard: above this n the pair loop is subsampled (with a note in
#'              the returned attribute) to keep the posterior loop tractable.
#' @return c(K_gonen = value)
#' @export
metric_gonen_k <- function(draw, w = NULL, max_n = 3000L) {
  eta <- draw$eta; n <- length(eta)
  if (is.null(w)) w <- rep(1 / n, n)
  idx <- seq_len(n)
  if (n > max_n) idx <- sort(sample.int(n, max_n))
  e <- eta[idx]; ww <- w[idx]
  d <- outer(e, e, `-`)                 # d[i,j] = eta_i - eta_j
  # I(eta_i > eta_j) / (1 + exp(eta_j - eta_i)) = I(d>0) * plogis(d)
  contrib <- (d > 0) * stats::plogis(d)
  W <- outer(ww, ww)
  diag(W) <- 0
  num <- sum(W * contrib)
  den <- sum(W)                          # sum over ordered pairs i != j
  # Each unordered pair contributes from one side; multiply by 2 to match the
  # A1 definition that sums both I(eta_i>eta_j) and I(eta_j>eta_i) terms.
  c(K_gonen = unname(2 * num / den))
}

#' Antolini's time-dependent concordance (per draw) — A1.4. Uses the model's
#' time-specific survival: a usable pair (i event at Ytilde_i, j later) is
#' concordant when S_i(Ytilde_i) < S_j(Ytilde_i), i.e. the subject who failed
#' had the lower predicted survival AT THAT TIME. For a PH model this reduces to
#' Harrell's C; for flexible / non-PH models whose curves may cross it is the
#' natural single summary. Survival at the (arbitrary) event times is obtained
#' by interpolating each subject's cumulative hazard across the horizon knots.
#'
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @param max_n guard: subsample the pair set above this n (kept tractable in
#'              the posterior loop; noted, not silent — see the returned NA path).
#' @return c(C_antolini = value)
#' @export
metric_antolini_c <- function(draw, w = NULL, max_n = 2000L) {
  n <- draw$n
  if (is.null(w)) w <- rep(1 / n, n)
  idx <- seq_len(n)
  if (n > max_n) idx <- sort(sample.int(n, max_n))
  time <- draw$time[idx]; status <- draw$status[idx]; ww <- w[idx]
  hz <- draw$horizons
  Lam <- .cumhaz_at_knots(draw$surv2d[idx, , drop = FALSE])   # m x (H+1)
  knots <- c(0, hz); H <- length(hz)

  num <- 0; den <- 0
  cases <- which(status == 1)
  for (i in cases) {
    ti <- time[i]
    later <- time > ti
    if (!any(later)) next
    Sall <- .surv_at_scalar(Lam, knots, H, ti)   # S_k(ti) for all k in idx
    Si <- Sall[i]; Sj <- Sall[later]; wj <- ww[later]
    conc <- (Si < Sj) + 0.5 * (Si == Sj)          # ties count half
    num <- num + ww[i] * sum(wj * conc)
    den <- den + ww[i] * sum(wj)
  }
  c(C_antolini = if (den > 0) unname(num / den) else NA_real_)
}

#' Time-dependent (cumulative/dynamic) AUC at each horizon, IPCW, per draw.
#' Cases: T_i <= t & event; controls: T_j > t. Uses the predicted risk
#' 1 - S_i(t) as the marker and IPCW weights G-hat from the censoring KM.
#'
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @return named vector AUC_t<h> for each horizon plus iAUC (trapezoid over t).
#' @export
metric_td_auc <- function(draw, w = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  time <- draw$time; status <- draw$status; hz <- draw$horizons
  Ghat <- .km_censoring(time, status)          # step function P(C > t)
  out <- numeric(length(hz))
  for (h in seq_along(hz)) {
    t   <- hz[h]
    risk <- 1 - draw$surv2d[, h]                 # marker at horizon t
    case <- (time <= t) & (status == 1)
    ctrl <- (time >  t)
    # IPCW: weight cases by 1/G(T_i-), controls by 1/G(t).
    gt   <- .step_eval(Ghat, t)
    gci  <- .step_eval(Ghat, time)
    gci[gci <= 0] <- NA
    wc   <- w * case / gci
    wk   <- w * ctrl / ifelse(gt <= 0, NA, gt)
    out[h] <- .weighted_auc(risk, wc, wk)
  }
  names(out) <- paste0("AUC_t", seq_along(hz))
  iauc <- if (length(hz) > 1)
    .trapz(hz, out) / (max(hz) - min(hz)) else out[1]
  c(out, iAUC = unname(iauc))
}

# --- helpers ---------------------------------------------------------------

# Monotone cumulative hazard Lam_i(t) = -log S_i(t) at the knots [0, horizons]
# for every subject (rows). cummax enforces monotonicity against numerical
# non-monotone survival draws. Shared by Antolini C and D-calibration.
#' @keywords internal
.cumhaz_at_knots <- function(surv2d) {
  Lam <- cbind(0, -log(pmax(surv2d, 1e-12)))
  t(apply(Lam, 1L, cummax))
}

# S_k(t) for ALL subjects at one scalar time t, by linear interpolation of the
# cumulative hazard between knots (extrapolates past the last knot with the last
# piece's slope). `Lam` is from .cumhaz_at_knots(); `knots = c(0, horizons)`.
#' @keywords internal
.surv_at_scalar <- function(Lam, knots, H, t) {
  if (t <= 0) return(rep(1, nrow(Lam)))
  k <- min(max(findInterval(t, knots), 1L), H)     # piece [knots[k], knots[k+1]]
  frac <- (t - knots[k]) / (knots[k + 1L] - knots[k])
  exp(-(Lam[, k] + frac * (Lam[, k + 1L] - Lam[, k])))
}

# S_i(t_i): each subject's survival at their OWN time t_i (vectorised over
# subjects with a shared knot vector).
#' @keywords internal
.surv_at_own <- function(Lam, knots, H, times) {
  n <- nrow(Lam)
  k <- pmin(pmax(findInterval(times, knots), 1L), H)
  lo <- Lam[cbind(seq_len(n), k)]; hi <- Lam[cbind(seq_len(n), k + 1L)]
  frac <- (times - knots[k]) / (knots[k + 1L] - knots[k])
  exp(-(lo + frac * (hi - lo)))
}

# Kaplan-Meier of the CENSORING distribution G(t) = P(C > t): flip the event
# indicator (censoring becomes the "event").
#' @keywords internal
.km_censoring <- function(time, status) {
  fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  list(time = fit$time, surv = fit$surv)
}

#' Right-continuous step evaluation of a KM survival curve at query times.
#' @keywords internal
.step_eval <- function(km, xout) {
  if (!length(km$time)) return(rep(1, length(xout)))
  stats::approx(x = c(0, km$time), y = c(1, km$surv), xout = xout,
                method = "constant", rule = 2, f = 0)$y
}

#' Weighted AUC = P(marker_case > marker_control), Mann-Whitney form with
#' case/control weights (ties count 1/2). Vectorised over the smaller set.
#' @keywords internal
.weighted_auc <- function(marker, wc, wk) {
  ok_c <- which(!is.na(wc) & wc > 0)
  ok_k <- which(!is.na(wk) & wk > 0)
  if (!length(ok_c) || !length(ok_k)) return(NA_real_)
  mc <- marker[ok_c]; wcv <- wc[ok_c]
  mk <- marker[ok_k]; wkv <- wk[ok_k]
  num <- 0; den <- sum(wcv) * sum(wkv)
  for (a in seq_along(mc)) {
    gt <- mk <  mc[a]
    eq <- mk == mc[a]
    num <- num + wcv[a] * (sum(wkv[gt]) + 0.5 * sum(wkv[eq]))
  }
  num / den
}

#' Trapezoidal integral.
#' @keywords internal
.trapz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}
