# cache.R
# ---------------------------------------------------------------------------
# Optional disk caching of the (expensive) posterior survival array. The `psurv`
# object is the costly artefact — an MCMC fit plus an [S x n x H] array — and is
# consumed unchanged by every metric (B2-B6). Caching it keyed on the inputs lets
# the worked example, the simulation, and interactive work skip re-fitting.
#
# Pure base R + `digest` (for a content hash; a serialize+md5 fallback is used if
# `digest` is absent). Caching is entirely opt-in; nothing else in the framework
# depends on it.
# ---------------------------------------------------------------------------

#' Content hash of an R object (digest if available, else serialize + md5sum).
#' @keywords internal
.hash_object <- function(x) {
  if (requireNamespace("digest", quietly = TRUE))
    return(digest::digest(x, algo = "md5"))
  tf <- tempfile(); on.exit(unlink(tf))
  con <- file(tf, "wb"); serialize(x, con); close(con)
  unname(tools::md5sum(tf))
}

#' Deterministic cache key for a psurv fit from its defining inputs.
#'
#' @param formula,data,horizons,backend,S the fit specification.
#' @param extra any additional backend arguments that change the result
#'        (priors, n_pieces, seed, ...).
#' @return a short hex string.
#' @export
psurv_cache_key <- function(formula, data, horizons, backend, S, extra = NULL) {
  .hash_object(list(formula = paste(deparse(formula), collapse = " "),
                    data = data, horizons = as.numeric(horizons),
                    backend = backend, S = S, extra = extra))
}

#' Fit a psurv with transparent disk caching.
#'
#' On a cache hit the stored psurv is returned immediately; on a miss the model
#' is fitted via fit_bayes_survival() and written to `<cache_dir>/psurv-<key>.rds`.
#' The key is derived from the inputs (or supply your own via `key`).
#'
#' @inheritParams fit_bayes_survival
#' @param cache_dir directory for cache files (created if needed).
#' @param key       explicit cache key (default: hash of the inputs).
#' @param refresh   TRUE forces a re-fit even on a cache hit.
#' @param verbose   message on hit/miss.
#' @param ...       passed to fit_bayes_survival() (and folded into the key).
#' @return a `psurv` object.
#' @export
fit_bayes_survival_cached <- function(formula, data, horizons,
                                      backend = "coxph", S = 500L, ...,
                                      cache_dir = ".psurv_cache", key = NULL,
                                      refresh = FALSE, verbose = TRUE) {
  if (is.null(key))
    key <- psurv_cache_key(formula, data, horizons, backend, S,
                           extra = list(...))
  path <- file.path(cache_dir, paste0("psurv-", key, ".rds"))
  if (!refresh && file.exists(path)) {
    if (verbose) message("psurv cache hit: ", path)
    ps <- readRDS(path)
    if (inherits(ps, "psurv")) return(ps)
    if (verbose) message("  (stale/invalid cache; re-fitting)")
  }
  ps <- fit_bayes_survival(formula, data, horizons, backend = backend, S = S,
                           ...)
  save_psurv(ps, path)
  if (verbose) message("psurv cached: ", path)
  ps
}

#' Save a psurv to disk (creates the directory).
#' @param ps a `psurv`.
#' @param path destination .rds path.
#' @export
save_psurv <- function(ps, path) {
  stopifnot(inherits(ps, "psurv"))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(ps, path)
  invisible(path)
}

#' Load a psurv from disk (validates the class).
#' @param path an .rds file written by save_psurv() / the cache.
#' @export
load_psurv <- function(path) {
  ps <- readRDS(path)
  stopifnot(inherits(ps, "psurv"))
  ps
}

#' Remove cached psurv files.
#' @param cache_dir the cache directory.
#' @export
clear_psurv_cache <- function(cache_dir = ".psurv_cache") {
  fs <- list.files(cache_dir, pattern = "^psurv-.*\\.rds$", full.names = TRUE)
  unlink(fs)
  invisible(length(fs))
}
