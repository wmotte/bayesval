# fit_stan_pwexp.R
# ---------------------------------------------------------------------------
# B1 (primary, genuine full-Bayes). A piecewise-exponential proportional-hazards
# survival model fitted by MCMC with rstan (stan/pwexp_ph.stan). Because the CRAN
# release of rstanarm does NOT ship stan_surv (it lives only on an unmerged
# GitHub branch), this hand-written model is the framework's real Bayesian PH
# analogue of Cox. It returns the same model-agnostic `psurv` object as every
# other backend, so B2-B6 are unchanged.
# ---------------------------------------------------------------------------

# Compiled-model cache (compile once per session; persisted to stan/pwexp_ph.rds).
.stan_pwexp_model <- local({
  sm <- NULL
  function(stan_file = NULL) {
    if (!is.null(sm)) return(sm)
    if (!requireNamespace("rstan", quietly = TRUE))
      stop("fit_stan_pwexp() needs the 'rstan' package.")
    if (is.null(stan_file)) {
      stan_file <- system.file("stan", "pwexp_ph.stan", package = "bayesval")
      if (!nzchar(stan_file)) {
        root <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
        stan_file <- file.path(root, "..", "stan", "pwexp_ph.stan")
        if (!file.exists(stan_file)) stan_file <- file.path("stan", "pwexp_ph.stan")
      }
    }
    rds <- sub("\\.stan$", ".rds", stan_file)
    if (file.exists(rds)) {
      obj <- tryCatch(readRDS(rds), error = function(e) NULL)
      if (!is.null(obj)) { sm <<- obj; return(sm) }
    }
    sm <<- rstan::stan_model(file = stan_file)
    tryCatch(saveRDS(sm, rds), error = function(e) NULL)
    sm
  }
})

#' Subject exposure per baseline-hazard piece, and the event-containing piece.
#' @keywords internal
.pwexp_expose <- function(time, status, cuts) {
  edges <- c(0, cuts, Inf)
  K <- length(cuts) + 1L
  lo <- edges[seq_len(K)]; hi <- edges[seq_len(K) + 1L]
  E <- t(vapply(time, function(t) pmax(0, pmin(t, hi) - lo), numeric(K)))
  kev <- pmin(pmax(findInterval(time, edges, rightmost.closed = FALSE), 1L), K)
  list(E = E, kev = as.integer(kev), edges = edges, K = K)
}

#' Cumulative baseline hazard H0(t_h) per draw: exp(log_h0) %*% W, with
#' `W[k,h]` = width of piece k within `[0, t_h]`.
#' @keywords internal
.pwexp_H0 <- function(log_h0_draws, edges, horizons) {
  K <- ncol(log_h0_draws)
  lo <- edges[seq_len(K)]; hi <- edges[seq_len(K) + 1L]
  W <- vapply(horizons, function(t) pmax(0, pmin(t, hi) - lo), numeric(K)) # K x H
  exp(log_h0_draws) %*% W                                                  # S x H
}

#' Fit a full-Bayes piecewise-exponential PH model and return a rich object.
#'
#' @param formula   Surv(time, status) ~ covariates.
#' @param data      training data.frame.
#' @param horizons  `numeric[H]` evaluation horizons.
#' @param n_pieces  baseline-hazard pieces (cuts at event-time quantiles).
#' @param beta_scale prior SD on log-HRs (Normal(0, beta_scale)); plan default 2.5.
#' @param chains,iter,warmup,cores,seed  rstan sampling controls.
#' @param S         optional thinning to S posterior draws for the psurv array.
#' @param refresh   rstan::sampling() progress-print interval (0 = silent).
#' @param adapt_delta target HMC acceptance rate (higher = smaller steps, fewer
#'        divergences, slower sampling); passed to rstan's `control` list.
#' @return list with $psurv (internal), $fit (stanfit), $beta_draws, $log_h0_draws,
#'   $cuts, $edges, $terms/$xlevels/$contrasts (for external prediction), and
#'   $diagnostics (max Rhat, min ESS, divergences).
#' @export
fit_stan_pwexp <- function(formula, data, horizons, n_pieces = 8L,
                           beta_scale = 2.5, chains = 4L, iter = 1000L,
                           warmup = floor(iter / 2), cores = chains,
                           seed = 1L, S = NULL, refresh = 0L,
                           adapt_delta = 0.95) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("fit_stan_pwexp() needs the 'rstan' package.")
  mf <- stats::model.frame(formula, data)
  y  <- stats::model.response(mf)
  time <- as.numeric(y[, 1L]); status <- as.integer(y[, ncol(y)])
  Terms <- stats::terms(mf)
  X <- stats::model.matrix(Terms, mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  # Standardise covariates: unstandardised predictors on very different scales
  # give a badly-conditioned posterior geometry and MCMC fails to mix. Center &
  # scale here; the SAME transform is applied to external data at prediction.
  ctr <- colMeans(X)
  scl <- apply(X, 2L, stats::sd); scl[scl == 0 | !is.finite(scl)] <- 1
  Xs <- sweep(sweep(X, 2L, ctr, `-`), 2L, scl, `/`)

  ev <- time[status == 1]; if (!length(ev)) ev <- time
  cuts <- unique(as.numeric(stats::quantile(
    ev, seq_len(n_pieces - 1L) / n_pieces, names = FALSE)))
  cuts <- cuts[cuts > 0 & is.finite(cuts)]
  ex <- .pwexp_expose(time, status, cuts)

  standata <- list(N = nrow(Xs), P = ncol(Xs), K = ex$K, X = Xs, E = ex$E,
                   d = status, kev = ex$kev, beta_scale = beta_scale)
  sm <- .stan_pwexp_model()
  fit <- rstan::sampling(sm, data = standata, chains = chains, iter = iter,
                         warmup = warmup, cores = cores, seed = seed,
                         refresh = refresh,
                         control = list(adapt_delta = adapt_delta))

  # Extract beta and log_h0 in a SINGLE permuted call: rstan::extract() draws a
  # fresh permutation on every call, so two separate permuted=TRUE calls return
  # MISALIGNED draws (beta row s and log_h0 row s would be different posterior
  # draws). One call shares the permutation, keeping (beta, log_h0) aligned per
  # draw — essential because the centering shift and the survival array pair them.
  .ext <- rstan::extract(fit, pars = c("beta", "log_h0"), permuted = TRUE)
  beta   <- as.matrix(.ext$beta)
  log_h0 <- as.matrix(.ext$log_h0)
  # Back-transform draws from the standardised sampling space to the ORIGINAL
  # covariate units, so eta_i = X_raw_i %*% beta and S(t) are reproduced directly
  # from (beta_draws, log_h0_draws) with the raw design (no stored scaling
  # needed downstream). beta_raw = beta_std / scl; the per-draw baseline shift
  # c = beta_raw %*% ctr absorbs the centering into the baseline log-hazards.
  beta <- sweep(beta, 2L, scl, `/`)
  shift <- as.numeric(beta %*% ctr)                 # length = #draws
  log_h0 <- log_h0 - shift                          # recycle by column
  colnames(beta) <- colnames(X)
  if (!is.null(S) && S < nrow(beta)) {
    idx <- unique(round(seq(1, nrow(beta), length.out = S)))
    beta <- beta[idx, , drop = FALSE]; log_h0 <- log_h0[idx, , drop = FALSE]
  }

  diag <- .stan_diagnostics(fit)
  obj <- list(fit = fit, beta_draws = beta, log_h0_draws = log_h0,
              cuts = cuts, edges = ex$edges, terms = stats::delete.response(Terms),
              xlevels = stats::.getXlevels(Terms, mf),
              contrasts = attr(X, "contrasts"), colnames = colnames(X),
              horizons = horizons, beta_scale = beta_scale, diagnostics = diag)
  class(obj) <- "stan_pwexp"
  # draws are in raw covariate units, so the psurv uses the raw design directly.
  obj$psurv <- .stan_pwexp_psurv(obj, X, time, status, horizons,
                                 label = "stan_pwexp")
  attr(obj$psurv, "fit") <- fit                    # for mcmc_diagnostics(ps)
  obj
}

#' Build a psurv from stored draws on a given design matrix X and response.
#' @keywords internal
.stan_pwexp_psurv <- function(obj, X, time, status, horizons, label) {
  eta <- obj$beta_draws %*% t(X)                       # S x n
  H0  <- .pwexp_H0(obj$log_h0_draws, obj$edges, horizons)  # S x H
  S <- nrow(eta); n <- ncol(eta)
  surv <- array(NA_real_, dim = c(S, n, length(horizons)))
  for (s in seq_len(S)) surv[s, , ] <- exp(-outer(exp(eta[s, ]), H0[s, ]))
  new_psurv(surv, eta, time, status, horizons, label = label)
}

#' External-validation psurv from a fitted stan_pwexp model on new data.
#'
#' Uses the training posterior draws of beta and the baseline log-hazards; only
#' the covariate design changes (train on cohort A, validate on cohort B, B3).
#' @param obj      a fitted `stan_pwexp` object from fit_stan_pwexp().
#' @param newdata  the external validation data.frame.
#' @param time,status names of the time / status columns in `newdata`.
#' @param horizons `numeric[H]` evaluation horizons (default: the training horizons).
#' @export
stan_pwexp_predict <- function(obj, newdata, time, status, horizons = NULL) {
  stopifnot(inherits(obj, "stan_pwexp"))
  if (is.null(horizons)) horizons <- obj$horizons
  mf <- stats::model.frame(obj$terms, data = newdata, xlev = obj$xlevels)
  X  <- stats::model.matrix(obj$terms, mf, contrasts.arg = obj$contrasts)
  X  <- X[, obj$colnames, drop = FALSE]
  .stan_pwexp_psurv(obj, X, newdata[[time]], newdata[[status]], horizons,
                    label = "stan_pwexp(external)")
}

#' MCMC diagnostics gate: max Rhat, min bulk-ESS, divergent transitions.
#' @keywords internal
.stan_diagnostics <- function(fit) {
  s <- rstan::summary(fit)$summary
  pars <- grep("^(beta|log_h0|sigma_rw)", rownames(s))
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divs <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  list(max_rhat = max(s[pars, "Rhat"], na.rm = TRUE),
       min_ess  = min(s[pars, "n_eff"], na.rm = TRUE),
       divergences = divs)
}

#' @export
print.stan_pwexp <- function(x, ...) {
  d <- x$diagnostics
  cat(sprintf("<stan_pwexp: %d draws, %d pieces, %d covariates>\n",
              nrow(x$beta_draws), length(x$cuts) + 1L, ncol(x$beta_draws)))
  cat(sprintf("  MCMC: max Rhat = %.3f, min ESS = %.0f, divergences = %d\n",
              d$max_rhat, d$min_ess, d$divergences))
  invisible(x)
}
