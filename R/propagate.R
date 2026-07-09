# propagate.R
# ---------------------------------------------------------------------------
# B2 — the single estimation principle of the framework.
#
# Any performance estimator M (Harrell C, Uno C, ICI, IBS, net benefit, ...) is
# a functional M(S^(s), data ; g^(s)) of the predicted-survival object and the
# data. For each posterior draw s we evaluate M with an existing estimator; the
# collection {M^(s)} is the posterior distribution induced for that finite-sample
# estimator. No new point estimators are invented.
#
# B3 — two variance sources. Passing weights = "bootstrap" attaches a fresh
# Dirichlet(1,...,1) weight vector g^(s) over the n evaluation subjects at each
# draw (Rubin's Bayesian bootstrap). The resulting posterior propagates BOTH
# model-parameter uncertainty (the draws) and finite evaluation-sample
# uncertainty (the weights) coherently.
# ---------------------------------------------------------------------------

#' Propagate a per-draw metric across the posterior (B2), optionally with the
#' Bayesian bootstrap (B3).
#'
#' @param ps        a `psurv` object.
#' @param metric_fn function(draw, w) -> numeric scalar OR named numeric vector,
#'                  where `draw` is the list from psurv_draw() and `w` is a
#'                  length-n weight vector (sums to 1). Metric functions may
#'                  ignore `w` if they cannot use case weights.
#' @param weights   "none" (parameter uncertainty only), "bootstrap" (add a
#'                  Dirichlet(1) reweight per draw), or a user matrix `[S x n]`.
#' @param draws     optional integer subset / count of draws to use (thinning).
#' @param seed      RNG seed for the bootstrap weights (reproducibility).
#' @param progress  print a dot every 10% of draws.
#' @return a `posterior_metric` object: matrix `[S_used x K]` of draws plus the
#'         metric names; summarise with summary(). The summary reports
#'         non-finite draw counts and flags metrics with more than 5% non-finite
#'         draws by default.
#' @export
propagate_metric <- function(ps, metric_fn,
                             weights = c("none", "bootstrap"),
                             draws = NULL, seed = 1L, progress = FALSE) {
  stopifnot(inherits(ps, "psurv"))
  metric_label <- paste(deparse(substitute(metric_fn)), collapse = " ")
  user_w <- NULL
  if (is.matrix(weights)) { user_w <- weights; wmode <- "user" }
  else { wmode <- match.arg(weights) }

  s_idx <- .resolve_draws(ps$S, draws)
  set.seed(seed)

  out <- vector("list", length(s_idx))
  tick <- max(1L, length(s_idx) %/% 10L)
  for (k in seq_along(s_idx)) {
    s  <- s_idx[k]
    dr <- psurv_draw(ps, s)
    w  <- switch(wmode,
                 none      = rep(1 / ps$n, ps$n),
                 bootstrap = .dirichlet1(ps$n),
                 user      = { v <- user_w[s, ]; v / sum(v) })
    val <- metric_fn(dr, w)
    out[[k]] <- val
    if (progress && k %% tick == 0L) cat(".")
  }
  if (progress) cat("\n")

  nm <- names(out[[1L]])
  mat <- do.call(rbind, lapply(out, function(v) as.numeric(v)))
  if (!is.null(nm)) colnames(mat) <- nm
  structure(mat, class = c("posterior_metric", "matrix"),
            weights = wmode, label = ps$label, seed = seed,
            metric = metric_label, package_version = .bayesval_version())
}

#' Installed package version, for attaching to propagated results (falls back
#' to "dev" when the package is loaded via `source()` and not installed).
#' @keywords internal
.bayesval_version <- function() {
  v <- tryCatch(as.character(utils::packageVersion("bayesval")),
                error = function(e) NA_character_)
  if (is.na(v)) "dev" else v
}

#' Resolve a draw subset: NULL = all; an integer = evenly thinned count; a
#' vector = explicit indices.
#' @keywords internal
.resolve_draws <- function(S, draws) {
  if (is.null(draws)) return(seq_len(S))
  if (length(draws) == 1L) {
    if (draws >= S) return(seq_len(S))
    return(unique(round(seq(1, S, length.out = draws))))
  }
  stopifnot(all(draws >= 1 & draws <= S))
  draws
}

#' Convenience: propagate several named metric functions in one pass over the
#' draws (cheaper than one pass each — the survival object is only sliced once).
#'
#' @param ps          a `psurv`.
#' @param metric_fns  named list of function(draw, w) -> scalar/named-vector.
#' @inheritParams propagate_metric
#' @export
propagate_many <- function(ps, metric_fns,
                           weights = c("none", "bootstrap"),
                           draws = NULL, seed = 1L, progress = FALSE) {
  combined <- function(draw, w) {
    unlist(lapply(metric_fns, function(f) f(draw, w)), use.names = TRUE)
  }
  propagate_metric(ps, combined, weights = weights, draws = draws,
                   seed = seed, progress = progress)
}

# ---------------------------------------------------------------------------
# Summaries (median + credible interval) — B2 reporting object.
# ---------------------------------------------------------------------------

#' @export
summary.posterior_metric <- function(object, probs = c(0.025, 0.5, 0.975),
                                      nonfinite_threshold = 0.05, ...) {
  m <- unclass(object)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
  nm <- colnames(m) %||% paste0("V", seq_len(ncol(m)))
  med_i <- which.min(abs(probs - 0.5))
  # One tidy row per metric column (K is small — this is not a hot per-draw loop).
  res <- purrr::map_dfr(seq_len(ncol(m)), function(j) {
    col <- m[, j]
    finite <- is.finite(col)
    n_draws <- length(col)
    n_finite <- sum(finite)
    n_nonfinite <- n_draws - n_finite
    pct_nonfinite <- n_nonfinite / n_draws
    q <- if (n_finite > 0L) {
      stats::quantile(col[finite], probs = probs, names = FALSE)
    } else {
      rep(NA_real_, length(probs))
    }
    avg <- if (n_finite > 0L) mean(col[finite]) else NA_real_
    tibble::tibble(metric = nm[j], mean = avg,
                   median = q[med_i], lower = q[1L], upper = q[length(q)],
                   n_draws = n_draws, n_finite = n_finite,
                   n_nonfinite = n_nonfinite,
                   pct_nonfinite = pct_nonfinite,
                   flag_unreliable = pct_nonfinite > nonfinite_threshold)
  })
  attr(res, "probs") <- probs
  attr(res, "weights") <- attr(object, "weights")
  attr(res, "nonfinite_threshold") <- nonfinite_threshold
  res
}

#' @export
print.posterior_metric <- function(x, ...) {
  cat(sprintf("<posterior_metric: %d draws x %d metric(s), weights = %s>\n",
              nrow(unclass(x)), ncol(unclass(x)),
              attr(x, "weights")))
  print(summary(x), digits = 4)
  invisible(x)
}
