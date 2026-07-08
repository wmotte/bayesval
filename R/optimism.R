# optimism.R
# ---------------------------------------------------------------------------
# B5 — optimism / overfitting correction, the Bayesian way.
#
#  * loo_survival(): PSIS-LOO elpd from a pointwise survival log-likelihood
#    matrix (e.g. rstanarm::log_lik(stan_surv_fit)) as an out-of-sample score.
#  * bb_optimism(): a Bayesian-bootstrap analogue of Harrell's rms::validate
#    optimism. Each Dirichlet(1) resample plays the role of a bootstrap sample:
#    fit on the resample, measure "apparent" performance on the resample and
#    "test" performance on the original data; optimism = apparent - test.
# ---------------------------------------------------------------------------

#' PSIS-LOO elpd from a pointwise log-likelihood matrix `[S x n]`.
#' @param log_lik matrix draws x subjects, e.g. rstanarm::log_lik(fit).
#' @return the loo object (contains elpd_loo, p_loo, pareto-k diagnostics).
#' @export
loo_survival <- function(log_lik) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("loo_survival() needs the 'loo' package.")
  stopifnot(is.matrix(log_lik))
  r_eff <- loo::relative_eff(exp(log_lik),
                             chain_id = rep(1L, nrow(log_lik)))
  loo::loo(log_lik, r_eff = r_eff)
}

#' Bayesian-bootstrap optimism correction for a scalar validation metric.
#'
#' @param fit_psurv function(data, weights) -> a `psurv` fitted on `data`
#'                  (weights are Dirichlet case weights; pass to the fitter or
#'                  ignore if it resamples internally).
#' @param data      the evaluation data.frame.
#' @param metric_fn a per-draw metric function metric_fn(draw, w) -> scalar
#'                  (single-value; if it returns a vector the first element is
#'                  used).
#' @param B         number of Bayesian-bootstrap replicates.
#' @param summarise function to collapse a metric's posterior draws to a point
#'                  (default median).
#' @param seed      RNG seed.
#' @return list with $apparent (on full data), $optimism (mean and per-rep),
#'         and $corrected = apparent - mean(optimism).
#' @export
bb_optimism <- function(fit_psurv, data, metric_fn, B = 50L,
                        summarise = stats::median, seed = 1L) {
  set.seed(seed)
  n <- nrow(data)

  point <- function(ps, w = NULL) {
    pm <- propagate_metric(ps, function(dr, ww)
      metric_fn(dr, if (is.null(w)) ww else w)[1], weights = "none")
    summarise(unclass(pm)[, 1L], na.rm = TRUE)
  }

  # Apparent performance of the model fitted on the full (unweighted) data.
  ps_full <- fit_psurv(data, rep(1 / n, n))
  apparent_full <- point(ps_full)

  optimism <- numeric(B)
  for (b in seq_len(B)) {
    g   <- .dirichlet1(n)                    # Bayesian-bootstrap weights
    psb <- fit_psurv(data, g)               # fit on the resample
    # apparent: performance measured on the resample (its own weights).
    app_b  <- point_weighted(psb, metric_fn, g, summarise)
    # test: same fitted model evaluated on the original data (uniform weights).
    test_b <- point_weighted(psb, metric_fn, rep(1 / n, n), summarise)
    optimism[b] <- app_b - test_b
  }

  corrected <- apparent_full - mean(optimism, na.rm = TRUE)
  list(apparent = apparent_full,
       optimism = list(mean = mean(optimism, na.rm = TRUE), reps = optimism),
       corrected = corrected)
}

#' Posterior point summary of a metric under fixed case weights.
#' @keywords internal
point_weighted <- function(ps, metric_fn, w, summarise = stats::median) {
  pm <- propagate_metric(ps, function(dr, ww) metric_fn(dr, w)[1],
                         weights = "none")
  summarise(unclass(pm)[, 1L], na.rm = TRUE)
}
