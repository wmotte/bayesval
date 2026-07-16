# decision_curve.R
# ---------------------------------------------------------------------------
# B6---the decision-analytic layer. Net benefit at threshold p_t (Vickers &
# Elkin), survival version: among subjects flagged positive (pi_i(t) >= p_t) the
# event rate is estimated by 1 - KM within that group at t.
#   NB(p_t) = f * r - f * (1 - r) * p_t / (1 - p_t)
# with f = P(flagged positive), r = event rate among the flagged. Propagated
# over posterior draws -> net-benefit bands and P(model beats treat-all/none).
# ---------------------------------------------------------------------------

# Net benefit lives on the scale of an event probability, so a difference below
# this is a floating-point tie and not a decision-analytic gain. Any comparison
# of a model net benefit against a default strategy must go through it.
.nb_tol <- 1e-8

#' Net benefit across thresholds for a single draw.
#' @param draw a `psurv_draw()` list for a single posterior draw.
#' @param thresholds numeric vector of risk thresholds p_t.
#' @param horizon index into `draw$horizons` (default: last).
#' @param w    length-n case-weight vector (sums to 1); `NULL` uses uniform weights.
#' @return list(model = `numeric[T]`, all = `numeric[T]`) net benefit per threshold.
#' @export
net_benefit_draw <- function(draw, thresholds, horizon = NULL, w = NULL) {
  if (is.null(w)) w <- rep(1 / draw$n, draw$n)
  h <- if (is.null(horizon)) draw$H else horizon
  t <- draw$horizons[h]
  risk <- 1 - draw$surv2d[, h]
  r_all <- 1 - .km_at(draw$time, draw$status, t, w)   # overall event rate

  nb_model <- numeric(length(thresholds))
  nb_all   <- numeric(length(thresholds))
  for (j in seq_along(thresholds)) {
    pt  <- thresholds[j]
    odds <- pt / (1 - pt)
    pos <- risk >= pt
    f   <- sum(w[pos])                       # weighted P(flagged positive)
    nb_all[j] <- r_all - (1 - r_all) * odds
    if (f <= 0) {
      nb_model[j] <- 0
    } else if (all(pos)) {
      # Flagging everyone IS treat-all. Recomputing it as f * r - f * (1-r) * odds
      # differs from nb_all in the last bit, because f is sum(w) = 1 only up to
      # floating-point summation error. That noise is ~1e-16 with an arbitrary
      # sign, and a downstream `NB_model > NB_all` test reads half of it as a
      # real benefit. Assign the identical value so the difference is exactly 0.
      nb_model[j] <- nb_all[j]
    } else {
      r <- 1 - .km_at(draw$time[pos], draw$status[pos], t, w[pos])
      nb_model[j] <- f * r - f * (1 - r) * odds
    }
  }
  list(model = nb_model, all = nb_all)
}

#' Posterior net-benefit curves over draws (B6).
#'
#' @param ps         a `psurv`.
#' @param thresholds risk thresholds p_t.
#' @param horizon    index into ps$horizons (default last).
#' @param weights    "none" or "bootstrap" (adds evaluation-sample uncertainty).
#' @param draws      optional draw subset/count.
#' @param seed       RNG seed.
#' @return list with $thresholds, $model (band df), $all (band df), and
#'         $p_useful = P(NB_model > max(NB_all, 0)) per threshold.
#' @export
posterior_net_benefit <- function(ps, thresholds = seq(0.05, 0.6, by = 0.05),
                                  horizon = NULL,
                                  weights = c("none", "bootstrap"),
                                  draws = NULL, seed = 1L) {
  stopifnot(inherits(ps, "psurv"))
  wmode <- match.arg(weights)
  s_idx <- .resolve_draws(ps$S, draws)
  set.seed(seed)

  M <- matrix(NA_real_, length(s_idx), length(thresholds))  # model NB
  A <- matrix(NA_real_, length(s_idx), length(thresholds))  # treat-all NB
  for (k in seq_along(s_idx)) {
    dr <- psurv_draw(ps, s_idx[k])
    w  <- if (wmode == "bootstrap") .dirichlet1(ps$n) else rep(1 / ps$n, ps$n)
    nb <- net_benefit_draw(dr, thresholds, horizon = horizon, w = w)
    M[k, ] <- nb$model
    A[k, ] <- nb$all
  }

  band <- function(X) tibble::tibble(
    threshold = thresholds,
    median = apply(X, 2L, stats::median, na.rm = TRUE),
    lower  = apply(X, 2L, stats::quantile, 0.025, na.rm = TRUE, names = FALSE),
    upper  = apply(X, 2L, stats::quantile, 0.975, na.rm = TRUE, names = FALSE))

  best_default <- pmax(A, 0)                            # treat-all or treat-none
  # Strictly greater, but only beyond .nb_tol: at low thresholds the model flags
  # every subject and ties the better default exactly, and a bare `>` would score
  # the rounding residue of that tie as usefulness in about half of those draws.
  p_useful <- colMeans(M > best_default + .nb_tol, na.rm = TRUE)

  list(thresholds = thresholds, weights = wmode,
       model = band(M), all = band(A),
       p_useful = stats::setNames(p_useful, thresholds),
       draws_model = M, draws_all = A)
}
