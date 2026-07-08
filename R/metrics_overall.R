# metrics_overall.R
# ---------------------------------------------------------------------------
# A3 -> B2. Per-draw overall-accuracy (proper) scores. The Brier score is a
# proper score combining discrimination AND calibration, with IPCW weights
# Ŵ_i(t) handling censoring (Graf et al.). Signature metric_fn(draw, w).
# ---------------------------------------------------------------------------

#' IPCW Brier score at each horizon (per draw).
#'   BS(t) = sum_i w_i * Ŵ_i(t) * (I(Ytilde_i > t) - S_i(t))^2
#' Ŵ_i(t): 1/G(T_i-) for events before t, 1/G(t) for survivors past t, 0 for
#' subjects censored before t.
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @param Ghat optional precomputed censoring KM (from .km_censoring); supplying
#'   it avoids recomputing the fixed censoring distribution on every draw.
#' @return named vector Brier_t<h> plus IBS (integrated over horizons).
#' @export
metric_brier <- function(draw, w = NULL, Ghat = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  time <- draw$time; status <- draw$status; hz <- draw$horizons
  if (is.null(Ghat)) Ghat <- .km_censoring(time, status)
  bs <- numeric(length(hz))
  for (h in seq_along(hz)) {
    t   <- hz[h]
    S_t <- draw$surv2d[, h]
    bs[h] <- .brier_ipcw(time, status, S_t, t, Ghat, w)
  }
  names(bs) <- paste0("Brier_t", seq_along(hz))
  ibs <- if (length(hz) > 1)
    .trapz(hz, bs) / (max(hz) - min(hz)) else bs[1]
  c(bs, IBS = unname(ibs))
}

#' Index of Prediction Accuracy (scaled Brier) at each horizon (per draw):
#'   IPA(t) = 1 - BS_model(t) / BS_null(t), null = marginal Kaplan-Meier.
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @param Ghat optional precomputed censoring KM (from .km_censoring); supplying
#'   it avoids recomputing the fixed censoring distribution on every draw.
#' @return named vector IPA_t<h>.
#' @export
metric_ipa <- function(draw, w = NULL, Ghat = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  time <- draw$time; status <- draw$status; hz <- draw$horizons
  if (is.null(Ghat)) Ghat <- .km_censoring(time, status)
  out <- numeric(length(hz))
  for (h in seq_along(hz)) {
    t   <- hz[h]
    S_model <- draw$surv2d[, h]
    S_null  <- rep(.km_at(time, status, t, w), draw$n)  # marginal KM
    bm <- .brier_ipcw(time, status, S_model, t, Ghat, w)
    bn <- .brier_ipcw(time, status, S_null,  t, Ghat, w)
    out[h] <- if (bn > 0) 1 - bm / bn else NA_real_
  }
  names(out) <- paste0("IPA_t", seq_along(hz))
  out
}

#' Royston-Sauerbrei R^2_D (per draw) — explained variation from the ordered
#' linear predictor. D is the Cox coefficient on scaled normal rankit scores.
#'   R^2_D = (D^2/kappa^2) / (pi^2/6 + D^2/kappa^2),  kappa = sqrt(8/pi).
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @return c(R2_D = value)
#' @export
metric_r2d <- function(draw, w = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  eta <- draw$eta; n <- length(eta)
  kappa <- sqrt(8 / pi)
  z <- stats::qnorm((rank(eta) - 0.375) / (n + 0.25)) / kappa
  df <- data.frame(time = draw$time, status = draw$status, z = z)
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, status) ~ z, data = df,
                    weights = w * n),
    error = function(e) NULL)
  if (is.null(fit)) return(c(R2_D = NA_real_))
  D <- unname(fit$coefficients[["z"]])
  s2 <- pi^2 / 6
  c(R2_D = (D^2 / kappa^2) / (s2 + D^2 / kappa^2))
}

# --- helpers ---------------------------------------------------------------

#' Single-horizon IPCW Brier.
#' @keywords internal
.brier_ipcw <- function(time, status, S_t, t, Ghat, w) {
  ev_before <- (time <= t) & (status == 1)   # observed event by t
  surv_past <- (time >  t)                    # known event-free at t
  # censored before t -> contributes weight 0 (unknown status at t).
  g_ti <- .step_eval(Ghat, time)
  g_t  <- .step_eval(Ghat, t)
  ipcw <- numeric(length(time))
  ipcw[ev_before] <- 1 / pmax(g_ti[ev_before], 1e-8)
  ipcw[surv_past] <- 1 / max(g_t, 1e-8)
  resid2 <- numeric(length(time))
  resid2[ev_before] <- (0 - S_t[ev_before])^2   # true survival indicator = 0
  resid2[surv_past] <- (1 - S_t[surv_past])^2   # indicator = 1
  ww <- w * ipcw
  sum(ww * resid2) / sum(w)                       # normalise by base weights
}
