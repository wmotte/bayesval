# config.R
# ---------------------------------------------------------------------------
# Central configuration for the framework (checklist 0): the global seed, the
# default evaluation horizons and truncation tau, the number of posterior draws
# S, and thinning. Everything reproducible flows from here so the worked
# example, the simulation, and the tests share one set of defaults.
#
# Access the singleton with bc_config(); override per call by passing named
# fields, e.g. bc_config(S = 1000L). Nothing here fits or predicts---it only
# holds settings.
# ---------------------------------------------------------------------------

#' The framework's default settings.
#'
#' @field seed      global RNG seed for reproducible draws/bootstraps.
#' @field horizons  default evaluation horizons t (NULL => choose from the data,
#'                  see bc_horizons()).
#' @field horizons_q if horizons is NULL, the event-time quantiles used to pick
#'                  data-driven horizons.
#' @field tau       truncation time for time-dependent measures (NULL => the
#'                  largest horizon).
#' @field S         number of posterior draws (coxph-Laplace backend / thinning
#'                  target for stan_surv).
#' @field thin      keep every `thin`-th draw when summarising heavy metrics.
#' @field bins_dcal number of cells for D-calibration.
#' @field thresholds default decision-curve threshold probabilities.
#' @field ci_probs  credible-interval probabilities (median + 95% CrI).
.BC_DEFAULTS <- list(
  seed       = 1L,
  horizons   = NULL,
  horizons_q = c(0.25, 0.50, 0.75),
  tau        = NULL,
  S          = 500L,
  thin       = 1L,
  bins_dcal  = 10L,
  thresholds = seq(0.05, 0.60, by = 0.05),
  ci_probs   = c(0.025, 0.50, 0.975)
)

#' Get the framework configuration, with optional per-call overrides.
#'
#' @param ... named fields to override (validated against the known settings).
#' @return a named list of settings.
#' @export
bc_config <- function(...) {
  ov <- list(...)
  cfg <- .BC_DEFAULTS
  if (length(ov)) {
    unknown <- setdiff(names(ov), names(cfg))
    if (length(unknown))
      stop("Unknown config field(s): ", paste(unknown, collapse = ", "))
    cfg[names(ov)] <- ov
  }
  cfg
}

#' Resolve evaluation horizons from data (or return the configured horizons).
#'
#' @param time,status observed follow-up and event indicator.
#' @param cfg a bc_config() list.
#' @return numeric horizons (event-time quantiles when cfg$horizons is NULL).
#' @export
bc_horizons <- function(time, status, cfg = bc_config()) {
  if (!is.null(cfg$horizons)) return(cfg$horizons)
  ev <- time[status == 1]
  if (!length(ev)) ev <- time
  as.numeric(stats::quantile(ev, cfg$horizons_q, na.rm = TRUE))
}

#' Resolve the truncation tau (largest horizon when cfg$tau is NULL).
#' @param horizons evaluation horizons.
#' @param cfg a bc_config() list.
#' @export
bc_tau <- function(horizons, cfg = bc_config()) {
  if (!is.null(cfg$tau)) cfg$tau else max(horizons)
}
