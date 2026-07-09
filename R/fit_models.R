# fit_models.R
# ---------------------------------------------------------------------------
# B1 -- the model layer. Fit a Bayesian survival model and return a `psurv`
# (posterior survival array). One dispatcher covers the primary Stan
# path and the `survival`-only Laplace/Cox fallback used for fast, dependency-
# light testing of the rest of the framework.
# ---------------------------------------------------------------------------

#' Fit a Bayesian survival model and return a `psurv`.
#'
#' @param formula survival formula, Surv(time, status) ~ covariates.
#' @param data    training/evaluation data.frame.
#' @param horizons `numeric[H]` evaluation horizons.
#' @param backend "stan_pwexp" (primary rstan piecewise-exponential PH),
#'                "stan_surv" (rstanarm M-spline PH if available), "pwexp"
#'                (piecewise-exponential Poisson-trick PH; survival-only),
#'                "bart" (nonparametric BART survival), or "coxph"
#'                (Cox-Laplace approximation; fallback / testing).
#' @param time,status column names (needed for the stan_surv prediction step).
#' @param S       number of posterior draws (coxph backend) / thinning target.
#' @param ...     passed to the backend fitter (priors, chains, cuts, etc.).
#' @export
fit_bayes_survival <- function(formula, data, horizons,
                               backend = c("stan_surv", "stan_pwexp", "coxph",
                                           "pwexp", "bart"),
                               time = NULL, status = NULL, S = 500L, ...) {
  backend <- match.arg(backend)
  if (backend == "coxph") {
    return(psurv_from_coxph(formula, data, horizons, S = S, ...))
  }
  if (backend == "pwexp") {
    return(psurv_from_pwexp(formula, data, horizons, S = S, ...))
  }
  if (backend == "stan_pwexp") {
    # Full-likelihood Bayesian PH via hand-written rstan model (fit_stan_pwexp.R);
    # this is the working primary Bayesian path since CRAN rstanarm ships no
    # stan_surv. Return its psurv with the full fit attached for loo/diagnostics.
    obj <- fit_stan_pwexp(formula, data, horizons, S = S, ...)
    ps <- obj$psurv
    attr(ps, "fit") <- obj$fit
    attr(ps, "stan_pwexp") <- obj
    return(ps)
  }
  if (backend == "bart") {
    if (is.null(time) || is.null(status))
      stop("Provide `time` and `status` column names for the bart backend.")
    return(psurv_from_bart(formula, data, horizons, time = time,
                           status = status, S = S, ...))
  }
  if (!requireNamespace("rstanarm", quietly = TRUE) ||
      !exists("stan_surv", where = asNamespace("rstanarm"), inherits = FALSE))
    stop("backend = 'stan_surv' needs rstanarm::stan_surv, which CRAN rstanarm ",
         "does not currently ship (it lives on an unmerged GitHub branch). Use ",
         "backend = 'stan_pwexp' (the primary full-likelihood path) or 'coxph' instead.")
  if (is.null(time) || is.null(status))
    stop("Provide `time` and `status` column names for the stan_surv backend.")
  stan_surv_fn <- get("stan_surv", envir = asNamespace("rstanarm"))
  fit <- stan_surv_fn(formula = formula, data = data, basehaz = "ms", ...)
  ps <- psurv_from_stan_surv(fit, newdata = data, time = time,
                             status = status, horizons = horizons, draws = S)
  attr(ps, "fit") <- fit                      # keep for loo/log_lik
  ps
}

#' Build a fit_psurv closure for bb_optimism(): function(data, weights) -> psurv.
#'
#' The Laplace/Cox backend accepts Dirichlet case weights directly, so the
#' Bayesian-bootstrap optimism loop reweights the fit coherently.
#'
#' @param formula,horizons,S as in fit_bayes_survival.
#' @param backend "coxph" (weighted refit) is the supported optimism backend.
#' @export
make_psurv_fitter <- function(formula, horizons, S = 300L,
                              backend = c("coxph")) {
  backend <- match.arg(backend)
  function(data, weights = NULL) {
    if (is.null(weights)) weights <- rep(1 / nrow(data), nrow(data))
    # coxph does non-standard evaluation on `weights`; pass the numeric vector
    # as a literal via do.call so it is not looked up in the data environment.
    wv <- weights * nrow(data)
    fit <- do.call(survival::coxph,
                   list(formula = formula, data = data, weights = wv,
                        x = TRUE, y = TRUE))
    .psurv_from_coxph_fit(fit, horizons = horizons, S = S)
  }
}

#' External-validation psurv: evaluate a coxph's Laplace posterior on NEW data.
#'
#' Fits nothing new---takes the training coxph, draws coefficients from its
#' Laplace posterior N(beta_hat, V) and uses its Breslow baseline, then predicts
#' S_i(t)^(s) and eta_i^(s) on `newdata`. This is how B3's external-set
#' validation design is realised (train on cohort A, validate on cohort B).
#'
#' @param fit      a survival::coxph fit (from the training cohort).
#' @param newdata  external validation data.frame (with the model covariates).
#' @param time,status column names of the response in newdata.
#' @param horizons `numeric[H]` evaluation horizons.
#' @param S,seed   posterior draws and RNG seed.
#' @export
psurv_external <- function(fit, newdata, time, status, horizons, S = 500L,
                           seed = 1L) {
  set.seed(seed)
  Terms <- stats::delete.response(stats::terms(fit))
  mf <- stats::model.frame(Terms, data = newdata, xlev = fit$xlevels)
  X  <- stats::model.matrix(Terms, mf, contrasts.arg = fit$contrasts)
  cn <- names(fit$coefficients)
  if (!all(cn %in% colnames(X)))
    stop("newdata design does not match the fitted model's coefficients.")
  X <- X[, cn, drop = FALSE]

  beta_draws <- .rmvn(S, fit$coefficients, as.matrix(stats::vcov(fit)))
  bh  <- survival::basehaz(fit, centered = FALSE)
  H0h <- .basehaz_at_horizons(bh, horizons)

  eta  <- beta_draws %*% t(X)
  surv <- array(NA_real_, dim = c(S, nrow(X), length(horizons)))
  for (s in seq_len(S)) surv[s, , ] <- exp(-outer(exp(eta[s, ]), H0h))
  new_psurv(surv, eta, newdata[[time]], newdata[[status]], horizons,
            x = newdata, label = "coxph-laplace(external)")
}

#' Internal: build a psurv from an already-fitted (possibly weighted) coxph.
#' Shares the Laplace-posterior logic with psurv_from_coxph().
#' @keywords internal
.psurv_from_coxph_fit <- function(fit, horizons, S = 300L) {
  X <- fit$x; n <- nrow(X)
  beta_hat <- fit$coefficients
  V <- as.matrix(stats::vcov(fit))
  beta_draws <- .rmvn(S, beta_hat, V)
  bh <- survival::basehaz(fit, centered = FALSE)
  H0h <- .basehaz_at_horizons(bh, horizons)
  y <- fit$y; time <- y[, 1L]; status <- y[, ncol(y)]
  eta <- beta_draws %*% t(X)
  surv <- array(NA_real_, dim = c(S, n, length(horizons)))
  for (s in seq_len(S)) surv[s, , ] <- exp(-outer(exp(eta[s, ]), H0h))
  new_psurv(surv, eta, time, status, horizons, label = "coxph-laplace(refit)")
}
