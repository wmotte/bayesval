# diagnostics.R
# ---------------------------------------------------------------------------
# B1 quality gates (checklist 1): the MCMC-convergence gate that must pass
# before any metric is propagated from an MCMC fit, and the prior-
# sensitivity harness. Both operate on the primary Bayesian path; the
# survival-only coxph/pwexp Laplace backends have no MCMC, so the gate reports
# "not applicable" for them and passes.
# ---------------------------------------------------------------------------

#' MCMC convergence diagnostics for a fitted Bayesian survival model.
#'
#' Extracts split-Rhat, bulk/tail ESS and the divergent-transition count. Works
#' on a stanreg (rstanarm) fit or anything the `posterior` package can coerce to
#' draws. For the Laplace backends (label starts with "coxph" / "pwexp") there is
#' no MCMC and the function returns a row with applicable = FALSE.
#'
#' @param fit a stan_surv/stanreg fit, OR a `psurv` carrying attr("fit").
#' @return list(applicable, max_rhat, min_ess_bulk, min_ess_tail, n_divergent,
#'         table = per-parameter data.frame).
#' @export
mcmc_diagnostics <- function(fit) {
  if (inherits(fit, "psurv")) {
    lbl <- fit$label %||% "model"
    f <- attr(fit, "fit")
    if (is.null(f)) {
      applicable <- !grepl("^(coxph|pwexp|bart)", lbl)
      return(list(applicable = applicable, reason = if (!applicable) NULL else
        "no stan fit attached", max_rhat = NA_real_, min_ess_bulk = NA_real_,
        min_ess_tail = NA_real_, n_divergent = NA_integer_, table = NULL))
    }
    fit <- f
  }
  # Unwrap the user's stan_pwexp object (or any list carrying a $fit stanfit).
  if (!methods::is(fit, "stanfit") && is.list(fit) &&
      methods::is(fit$fit, "stanfit")) fit <- fit$fit

  if (!requireNamespace("posterior", quietly = TRUE))
    stop("mcmc_diagnostics() needs the 'posterior' package for a stan fit.")

  drw <- posterior::as_draws_df(fit)
  smry <- posterior::summarise_draws(drw, "rhat", "ess_bulk", "ess_tail") |>
    tibble::as_tibble() |>
    dplyr::filter(!grepl("^\\.", variable))
  ndiv <- .n_divergent(fit)
  list(
    applicable   = TRUE,
    max_rhat     = max(smry$rhat, na.rm = TRUE),
    min_ess_bulk = min(smry$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(smry$ess_tail, na.rm = TRUE),
    n_divergent  = ndiv,
    table        = smry)
}

#' Convergence gate: TRUE only if Rhat < max_rhat, ESS >= min_ess, no divergences.
#'
#' Metrics computed on a fit that fails this gate must not be reported (plan
#' "MCMC diagnostics gate"). For a non-MCMC backend the gate passes vacuously
#' with a note.
#'
#' @param fit as in mcmc_diagnostics().
#' @param max_rhat maximum acceptable split-Rhat (default 1.01).
#' @param min_ess minimum acceptable bulk & tail ESS (default 400).
#' @param error   if TRUE, stop() on failure instead of returning FALSE.
#' @return logical (invisibly), with attribute "diagnostics".
#' @export
check_mcmc_gate <- function(fit, max_rhat = 1.01, min_ess = 400,
                            error = FALSE) {
  d <- mcmc_diagnostics(fit)
  if (!isTRUE(d$applicable)) {
    message("MCMC gate: not applicable (Laplace/analytic backend) - passing.")
    return(structure(TRUE, diagnostics = d))
  }
  ok <- is.finite(d$max_rhat) && d$max_rhat < max_rhat &&
        d$min_ess_bulk >= min_ess && d$min_ess_tail >= min_ess &&
        (is.na(d$n_divergent) || d$n_divergent == 0)
  msg <- sprintf(
    "MCMC gate %s: max Rhat = %.4f (<%.3f), min ESS = %.0f (>=%.0f), divergences = %s",
    if (ok) "PASS" else "FAIL", d$max_rhat, max_rhat,
    min(d$min_ess_bulk, d$min_ess_tail), min_ess,
    ifelse(is.na(d$n_divergent), "NA", d$n_divergent))
  if (!ok && error) stop(msg)
  message(msg)
  structure(ok, diagnostics = d)
}

#' Prior-sensitivity analysis (checklist 1): refit the primary Bayesian PH
#' model (fit_stan_pwexp) under several log-HR prior SDs and compare the
#' posterior of each performance estimator. A robust conclusion should be stable
#' across weakly-informative priors; a large shift flags prior-dependence.
#'
#' Uses `beta_scale` (the Normal(0, beta_scale) SD on the log-hazard-ratios) as
#' the prior knob---the working analogue of stan_surv's coefficient prior, since
#' CRAN rstanarm ships no stan_surv.
#'
#' @param formula,data,horizons as in fit_stan_pwexp.
#' @param metric_fns named list of per-draw metric functions to compare.
#' @param scales named numeric of prior SDs; default contrasts the plan default
#'        (2.5) with a tight (1.0) and a loose (5.0) log-HR prior.
#' @param ... extra args passed to fit_stan_pwexp (chains, iter, n_pieces, ...).
#' @param flag threshold on the posterior-median range across priors that
#'        triggers the prior-sensitive message (default 0.02).
#' @return data.frame: prior x metric x (median, lower, upper), plus a printed
#'         note if any metric's median moves by more than `flag` across priors.
#' @export
prior_sensitivity <- function(formula, data, horizons, metric_fns,
                              scales = c(tight = 1.0, default = 2.5,
                                         loose = 5.0), flag = 0.02, ...) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("prior_sensitivity() needs 'rstan' (fit_stan_pwexp, the primary path).")
  # Each prior refits the (expensive) Stan model, so the loop stays sequential;
  # only the posterior-summary assembly is tidied.
  res <- purrr::map(names(scales), function(pn) {
    obj <- fit_stan_pwexp(formula = formula, data = data, horizons = horizons,
                          beta_scale = scales[[pn]], ...)
    summary(propagate_many(obj$psurv, metric_fns)) |>
      dplyr::mutate(prior = pn, beta_scale = scales[[pn]]) |>
      dplyr::select(prior, beta_scale, metric, median, lower, upper)
  }) |>
    dplyr::bind_rows()
  # flag metrics whose posterior median is sensitive to the prior.
  sens <- res |>
    dplyr::group_by(metric) |>
    dplyr::summarise(range = diff(range(median)), .groups = "drop") |>
    dplyr::filter(range > flag) |>
    dplyr::pull(metric)
  if (length(sens))
    message("Prior-sensitive metric(s) (median range > ", flag, "): ",
            paste(sens, collapse = ", "))
  attr(res, "sensitive") <- sens
  res
}

#' Alias for check_mcmc_gate() (the convergence gate).
#' @rdname check_mcmc_gate
#' @export
mcmc_gate <- function(fit, max_rhat = 1.01, min_ess = 400, error = FALSE)
  check_mcmc_gate(fit, max_rhat = max_rhat, min_ess = min_ess, error = error)

# --- helpers ---------------------------------------------------------------

#' @keywords internal
.n_divergent <- function(fit) {
  if (!requireNamespace("rstan", quietly = TRUE)) return(NA_integer_)
  # a raw stanfit, a stanreg ($stanfit), or a stan_pwexp/list carrying $fit.
  sf <- if (methods::is(fit, "stanfit")) fit
        else if (!is.null(fit$stanfit)) fit$stanfit
        else if (!is.null(fit$fit) && methods::is(fit$fit, "stanfit")) fit$fit
        else NULL
  if (is.null(sf)) return(NA_integer_)
  sp <- tryCatch(rstan::get_sampler_params(sf, inc_warmup = FALSE),
                 error = function(e) NULL)
  if (is.null(sp)) return(NA_integer_)
  as.integer(sum(vapply(sp, function(m) sum(m[, "divergent__"]), numeric(1))))
}

#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
