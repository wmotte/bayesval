# ppc_calibration.R
# ---------------------------------------------------------------------------
# B4 — calibration as a posterior predictive check (model criticism), not only
# as a summarised number.
#
# For each posterior draw s: simulate a replicate outcome y_rep^(s) from the
# posterior predictive (event times from S_i(.)^(s); censoring resampled
# independently from the empirical censoring KM). Recompute the ICI on the
# replicate and on the real data (both using the SAME predicted risks pi^(s),
# i.e. a realized discrepancy in Gelman's sense). The distribution of replicated
# ICIs vs the observed ICI gives a posterior predictive p-value; the ensemble of
# observed calibration curves gives credible bands.
# ---------------------------------------------------------------------------

#' Posterior predictive calibration check at a horizon.
#'
#' @param ps       a `psurv`.
#' @param horizon  index into ps$horizons (default: last).
#' @param method   calibration-curve method passed to .calibration_curve.
#' @param draws    optional draw subset/count (PPC is heavier; thin freely).
#' @param grid     predicted-risk grid on which observed curves are evaluated
#'                 for the credible band.
#' @param seed     RNG seed.
#' @return list with $ici_obs, $ici_rep, $ppp (posterior predictive p-value),
#'         $band (data.frame grid + observed-curve median/lower/upper), and
#'         $curve_obs (observed calibration curve on real data, per grid point).
#' @export
ppc_calibration <- function(ps, horizon = NULL,
                            method = c("auto", "hare", "group"),
                            draws = NULL, grid = seq(0.02, 0.98, by = 0.02),
                            seed = 1L) {
  stopifnot(inherits(ps, "psurv"))
  method <- match.arg(method)
  h <- if (is.null(horizon)) ps$H else horizon
  t <- ps$horizons[h]
  s_idx <- .resolve_draws(ps$S, draws)
  set.seed(seed)

  Gcens <- .km_censoring(ps$time, ps$status)   # to resample censoring times
   w <- rep(1 / ps$n, ps$n)

  ici_obs <- numeric(length(s_idx))
  ici_rep <- numeric(length(s_idx))
  curves_obs <- matrix(NA_real_, length(s_idx), length(grid))
  curves_rep <- matrix(NA_real_, length(s_idx), length(grid))

  for (k in seq_along(s_idx)) {
    s  <- s_idx[k]
    dr <- psurv_draw(ps, s)
    pred_risk <- 1 - dr$surv2d[, h]

    # observed curve + ICI on the real data, at this draw's predictions.
    oc_obs <- .calibration_curve(pred_risk, ps$time, ps$status, t, w, method)
    ici_obs[k] <- stats::weighted.mean(abs(pred_risk - oc_obs), w)
    curves_obs[k, ] <- .curve_on_grid(pred_risk, oc_obs, grid)

    # replicate outcome from the posterior predictive, then recompute.
    rep_yd <- .simulate_replicate(dr$surv2d, ps$horizons, Gcens)
    oc_rep <- .calibration_curve(pred_risk, rep_yd$time, rep_yd$status, t, w,
                                 method)
    ici_rep[k] <- stats::weighted.mean(abs(pred_risk - oc_rep), w)
    curves_rep[k, ] <- .curve_on_grid(pred_risk, oc_rep, grid)
  }

  ppp <- mean(ici_rep >= ici_obs, na.rm = TRUE)
  band <- tibble::tibble(
    pred   = grid,
    median = apply(curves_obs, 2L, stats::median, na.rm = TRUE),
    lower  = apply(curves_obs, 2L, stats::quantile, 0.025, na.rm = TRUE,
                   names = FALSE),
    upper  = apply(curves_obs, 2L, stats::quantile, 0.975, na.rm = TRUE,
                   names = FALSE)
  )
  list(horizon = t, method = method,
       ici_obs = ici_obs, ici_rep = ici_rep, ppp = ppp,
       band = band, curves_rep = curves_rep)
}

# --- helpers ---------------------------------------------------------------

#' Simulate a replicate (time, status) from a per-subject discrete survival
#' curve (piecewise-exponential interpolation) with independently resampled
#' censoring from the empirical censoring KM.
#' @keywords internal
.simulate_replicate <- function(surv2d, hz, Gcens) {
  n <- nrow(surv2d)
  Tt <- .sim_event_times(surv2d, hz)
  Cc <- .sim_from_km(Gcens, n)
  time <- pmin(Tt, Cc)
  status <- as.integer(Tt <= Cc)
  list(time = time, status = status)
}

#' Inverse-CDF event times from a discrete survival curve via piecewise-constant
#' hazard between horizon knots (extended beyond the last knot).
#' @keywords internal
.sim_event_times <- function(surv2d, hz) {
  n <- nrow(surv2d)
  H <- length(hz)
  tgrid <- c(0, hz)
  E <- stats::rexp(n)                         # target cumulative hazard
  out <- numeric(n)
  for (i in seq_len(n)) {
    Lam <- c(0, -log(pmax(surv2d[i, ], 1e-12)))   # cumhaz at 0, hz...
    Lam <- cummax(Lam)                             # enforce monotone
    if (E[i] <= Lam[H + 1]) {
      out[i] <- stats::approx(x = Lam, y = tgrid, xout = E[i], rule = 2)$y
    } else {
      # extrapolate with the last interval's hazard slope.
      lam_last <- (Lam[H + 1] - Lam[H]) / (tgrid[H + 1] - tgrid[H])
      lam_last <- max(lam_last, 1e-8)
      out[i] <- tgrid[H + 1] + (E[i] - Lam[H + 1]) / lam_last
    }
  }
  out
}

#' Sample n times from a KM survival curve by inverse transform.
#' @keywords internal
.sim_from_km <- function(km, n) {
  if (!length(km$time)) return(rep(Inf, n))
  u <- stats::runif(n)
  # survival curve is decreasing in time; find t with S(t) = u.
  tt <- c(0, km$time); ss <- c(1, km$surv)
  # Survival can be flat (tied levels); collapse to the earliest time reaching
  # each level (ties = min) so the inverse transform is well-defined & quiet.
  # for u below the last observed S, censoring beyond follow-up -> +Inf.
  out <- stats::approx(x = ss, y = tt, xout = u, method = "linear",
                       rule = 2, ties = min)$y
  out[u < min(ss)] <- Inf
  out
}

#' Evaluate a (predicted -> observed) point cloud on a fixed grid by monotone
#' interpolation (for building calibration bands).
#' @keywords internal
.curve_on_grid <- function(pred, obs, grid) {
  ok <- is.finite(pred) & is.finite(obs)
  if (sum(ok) < 2L) return(rep(if (any(ok)) obs[which(ok)[1L]] else NA_real_,
                                length(grid)))
  o <- order(pred[ok])
  stats::approx(x = pred[ok][o], y = obs[ok][o], xout = grid, rule = 2,
                ties = mean)$y
}
