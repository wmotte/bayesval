# bayes_bootstrap.R
# ---------------------------------------------------------------------------
# B3 — the second variance source. The parameter posterior alone under-
# propagates uncertainty about the EVALUATION sample. Rubin's (1981) Bayesian
# bootstrap attaches a Dirichlet(1,...,1) weight vector g^(s) over the n
# evaluation subjects at each draw; the weighted metric M(S^(s), data; g^(s))
# then reflects BOTH model-parameter and finite-sample uncertainty.
#
# The mechanism lives in propagate_metric(weights = "bootstrap"); this file adds
# named entry points and a diagnostic for the double-counting risk flagged in
# the plan (parameter posterior + bootstrap may over-cover).
# ---------------------------------------------------------------------------

#' Generate a matrix of Bayesian-bootstrap weights `[S x n]` (rows sum to 1).
#' @param S,n dimensions.
#' @param seed RNG seed.
#' @export
bb_weights <- function(S, n, seed = 1L) {
  set.seed(seed)
  t(vapply(seq_len(S), function(i) .dirichlet1(n), numeric(n)))
}

#' Propagate a metric with the Bayesian bootstrap (B2 + B3 in one call).
#' Thin alias for propagate_metric(..., weights = "bootstrap").
#' @param ps        a `psurv` object.
#' @param metric_fn function(draw, w) -> numeric scalar OR named numeric vector.
#' @param draws     optional integer subset / count of draws to use (thinning).
#' @param seed      RNG seed for the bootstrap weights (reproducibility).
#' @param progress  print a dot every 10% of draws.
#' @export
propagate_bootstrap <- function(ps, metric_fn, draws = NULL, seed = 1L,
                                progress = FALSE) {
  propagate_metric(ps, metric_fn, weights = "bootstrap", draws = draws,
                   seed = seed, progress = progress)
}

#' Bayesian-bootstrap replicates of a metric on a SINGLE fixed prediction set.
#'
#' Isolates evaluation-set (B3) uncertainty when the model is treated as fixed
#' (e.g. the true model in the coverage study): the metric is recomputed under B
#' fresh Dirichlet(1) weightings of the n subjects.
#'
#' @param draw   a psurv_draw()-style list (fixed predictions + data).
#' @param metric_fn function(draw, w) -> scalar or named vector.
#' @param B      number of bootstrap replicates.
#' @param seed   RNG seed.
#' @return matrix `[B x K]` of metric values.
#' @export
bb_replicates <- function(draw, metric_fn, B = 500L, seed = 1L) {
  set.seed(seed)
  k <- length(metric_fn(draw, rep(1 / draw$n, draw$n)))
  # vapply returns k x B (or a length-B vector when k == 1); reshape to B x k.
  v <- vapply(seq_len(B), function(i)
    as.numeric(metric_fn(draw, .dirichlet1(draw$n))), numeric(k))
  matrix(as.vector(v), nrow = B, ncol = k, byrow = TRUE)
}

#' Decompose posterior variance into parameter vs evaluation-sample sources.
#'
#' Runs the same metric three ways and reports interval widths so the coverage
#' study (Part C) can decide whether to down-weight one source:
#'   * parameter-only  (weights = "none")
#'   * bootstrap-only  (fix the draw index, vary only Dirichlet weights)
#'   * combined        (weights = "bootstrap")
#'
#' @param ps        a `psurv` object.
#' @param metric_fn function(draw, w) -> numeric scalar OR named numeric vector.
#' @param draws     optional integer subset / count of draws to use (thinning).
#' @param seed      RNG seed for the bootstrap weights (reproducibility).
#' @return data.frame of source, sd, and 95% interval width per metric column.
#' @export
variance_decomposition <- function(ps, metric_fn, draws = NULL, seed = 1L) {
  param_only <- propagate_metric(ps, metric_fn, weights = "none",
                                 draws = draws, seed = seed)
  combined   <- propagate_metric(ps, metric_fn, weights = "bootstrap",
                                 draws = draws, seed = seed)

  # Bootstrap-only: hold the posterior at its mean draw, vary weights.
  s_star <- which.min(rowSums((ps$eta - matrix(colMeans(ps$eta), ps$S, ps$n,
                                               byrow = TRUE))^2))
  set.seed(seed + 1L)
  B <- if (is.null(draws)) ps$S else min(ps$S, max(draws))
  dr <- psurv_draw(ps, s_star)
  boot_only <- t(vapply(seq_len(B), function(i)
    as.numeric(metric_fn(dr, .dirichlet1(ps$n))),
    numeric(length(metric_fn(dr, rep(1 / ps$n, ps$n))))))

  width <- function(M) apply(M, 2L,
    function(x) diff(stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)))
  nm <- colnames(unclass(param_only))
  sources <- list(parameter  = unclass(param_only),
                  evaluation = boot_only,
                  combined   = unclass(combined))
  purrr::imap(sources, function(M, src)
    tibble::tibble(source = src, metric = nm,
                   sd = apply(M, 2L, stats::sd, na.rm = TRUE),
                   width95 = width(M))) |>
    dplyr::bind_rows()
}
