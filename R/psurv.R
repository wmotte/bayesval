# psurv.R
# ---------------------------------------------------------------------------
# The central data structure of the Bayesian validation framework (Part B).
#
# A `psurv` object is the tidy posterior survival array described in B1/B2:
# for every posterior draw s = 1..S, every subject i = 1..n and every horizon
# t_h it stores the predicted survival S_i(t_h)^(s) and the linear predictor
# eta_i^(s), together with the observed evaluation data (time, status, X).
#
# Every performance estimator in Part B is a *functional* of this object plus
# the data (B2). Keeping a single, model-agnostic representation is what makes the
# framework work identically for stan_surv, a piecewise-exponential model, BART
# survival, or the Laplace/Cox fallback used for fast testing.
# ---------------------------------------------------------------------------

#' Construct a posterior survival object.
#'
#' @param surv    array `[S x n x H]` of predicted survival S_i(t_h)^(s) in (0,1].
#' @param eta     matrix `[S x n]` of linear predictors eta_i^(s).
#' @param time    `numeric[n]` observed follow-up time Ytilde_i = min(T_i, C_i).
#' @param status  `integer[n]` event indicator Delta_i (1 = event, 0 = censored).
#' @param horizons `numeric[H]` the evaluation horizons t_h (same order as dim 3).
#' @param x       optional data.frame`[n x p]` of covariates (needed by some
#'                calibration/PPC steps that refit on eta).
#' @param label   character tag identifying the model that produced the draws.
#' @return an object of class "psurv".
#' @export
new_psurv <- function(surv, eta, time, status, horizons, x = NULL,
                      label = "model") {
  stopifnot(is.array(surv), length(dim(surv)) == 3L)
  S <- dim(surv)[1L]; n <- dim(surv)[2L]; H <- dim(surv)[3L]
  stopifnot(is.matrix(eta), nrow(eta) == S, ncol(eta) == n)
  stopifnot(length(time) == n, length(status) == n, length(horizons) == H)
  if (any(surv < 0 | surv > 1, na.rm = TRUE))
    stop("surv must lie in [0, 1].")
  structure(
    list(surv = surv, eta = eta,
         time = as.numeric(time), status = as.integer(status),
         horizons = as.numeric(horizons), x = x,
         S = S, n = n, H = H, label = label),
    class = "psurv"
  )
}

#' @export
print.psurv <- function(x, ...) {
  cat(sprintf("<psurv: %s>\n", x$label))
  cat(sprintf("  draws S = %d, subjects n = %d, horizons H = %d\n",
              x$S, x$n, x$H))
  cat(sprintf("  horizons: %s\n",
              paste(signif(x$horizons, 4), collapse = ", ")))
  cat(sprintf("  events: %d (%.1f%% censored)\n",
              sum(x$status), 100 * mean(x$status == 0)))
  invisible(x)
}

#' Pull a single draw as a light list, for per-draw metric functions.
#' @keywords internal
psurv_draw <- function(ps, s) {
  list(
    surv   = ps$surv[s, , , drop = FALSE][1, , , drop = FALSE], # [n x H]
    surv2d = matrix(ps$surv[s, , ], nrow = ps$n, ncol = ps$H),
    eta    = ps$eta[s, ],
    time   = ps$time, status = ps$status,
    horizons = ps$horizons, x = ps$x, n = ps$n, H = ps$H
  )
}

# ---------------------------------------------------------------------------
# Constructors (B1) — turn a fitted model into a psurv object.
# ---------------------------------------------------------------------------

#' Build a psurv from an rstanarm::stan_surv fit (the primary Bayesian PH model).
#'
#' Uses posterior_survfit() to obtain S_i(t)^(s) and the posterior linear
#' predictor. Requires rstanarm; kept thin so the model layer is swappable.
#'
#' @param fit      a stan_surv fit.
#' @param newdata  evaluation data (defaults to the fit's model frame).
#' @param time,status names of the time / status columns in newdata.
#' @param horizons `numeric[H]` evaluation horizons.
#' @param draws    optional number of posterior draws to thin to.
#' @export
psurv_from_stan_surv <- function(fit, newdata, time, status, horizons,
                                 draws = NULL) {
  if (!requireNamespace("rstanarm", quietly = TRUE))
    stop("psurv_from_stan_surv() needs the 'rstanarm' package.")
  n <- nrow(newdata)
  # Posterior linear predictor: draws x n.
  eta <- rstanarm::posterior_linpred(fit, newdata = newdata, draws = draws)
  S <- nrow(eta)
  H <- length(horizons)
  surv <- array(NA_real_, dim = c(S, n, H))
  for (h in seq_len(H)) {
    # standardise = FALSE -> subject-specific survival at horizon h.
    ps <- rstanarm::posterior_survfit(
      fit, newdata = newdata, times = horizons[h],
      standardise = FALSE, draws = draws, return_matrix = TRUE
    )
    # ps is a list (one matrix per subject) or matrix draws x n; normalise.
    surv[, , h] <- .as_draws_by_subject(ps, S = S, n = n)
  }
  new_psurv(surv, eta, newdata[[time]], newdata[[status]], horizons,
            x = newdata, label = "stan_surv")
}

# rstanarm's return shapes vary by version; coerce to draws x n.
.as_draws_by_subject <- function(ps, S, n) {
  if (is.list(ps) && !is.data.frame(ps)) {
    out <- vapply(ps, function(m) as.numeric(m), numeric(S))
    return(matrix(out, nrow = S, ncol = n))
  }
  m <- as.matrix(ps)
  if (nrow(m) == S && ncol(m) == n) return(m)
  if (nrow(m) == n && ncol(m) == S) return(t(m))
  stop("Unexpected posterior_survfit shape.")
}

#' Fallback constructor: posterior draws from a frequentist Cox fit.
#'
#' This is the framework's *testing harness* and a legitimate large-sample
#' (Laplace / Bayesian-central-limit) approximation to the posterior of a Cox
#' model under flat priors: coefficients are drawn from N(beta_hat, vcov_hat)
#' and combined with the Breslow baseline cumulative hazard. It lets the entire
#' Part B pipeline (B2-B6) run and be unit-tested with only the `survival`
#' package installed — no Stan toolchain required. The primary manuscript path
#' remains psurv_from_stan_surv(); this documents the model-agnostic claim.
#'
#' @param formula a survival::coxph formula, e.g. Surv(time, status) ~ x1 + x2.
#' @param data    the evaluation data.frame.
#' @param horizons `numeric[H]` evaluation horizons.
#' @param S       number of posterior draws.
#' @param baseline "fixed" (Breslow at the MLE) or "bootstrap" (Bayesian
#'                bootstrap of the baseline hazard for extra coherence).
#' @param seed    RNG seed for reproducibility.
#' @export
psurv_from_coxph <- function(formula, data, horizons, S = 500L,
                             baseline = c("fixed", "bootstrap"), seed = 1L) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("psurv_from_coxph() needs the 'survival' package.")
  baseline <- match.arg(baseline)
  set.seed(seed)

  fit <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)
  X   <- fit$x                              # n x p model matrix (centered below)
  n   <- nrow(X); p <- ncol(X)
  beta_hat <- fit$coefficients
  V   <- as.matrix(stats::vcov(fit))

  # Draw coefficients from the Laplace posterior N(beta_hat, V).
  beta_draws <- .rmvn(S, beta_hat, V)       # S x p

  # Breslow baseline cumulative hazard H0(t) at the MLE, evaluated at horizons.
  bh  <- survival::basehaz(fit, centered = FALSE)  # cols: hazard, time
  H0h <- stats::approx(x = bh$time, y = bh$hazard, xout = horizons,
                       method = "constant", rule = 2, f = 0)$y
  H0h[is.na(H0h)] <- 0

  # Observed response.
  y      <- fit$y
  time   <- y[, 1L]
  status <- y[, ncol(y)]

  eta  <- beta_draws %*% t(X)               # S x n linear predictors
  surv <- array(NA_real_, dim = c(S, n, length(horizons)))
  for (s in seq_len(S)) {
    H0s <- H0h
    if (baseline == "bootstrap") {
      # Bayesian bootstrap on the baseline: Dirichlet reweight of increments.
      w <- .dirichlet1(length(H0h))
      H0s <- cumsum(diff(c(0, H0h)) * w * length(H0h))
    }
    lin <- exp(eta[s, ])                     # n
    # S_i(t_h) = exp(-H0(t_h) * exp(eta_i))
    surv[s, , ] <- exp(-outer(lin, H0s))
  }

  xdf <- as.data.frame(X)
  xdf[[".time"]] <- time; xdf[[".status"]] <- status
  new_psurv(surv, eta, time, status, horizons, x = xdf,
            label = sprintf("coxph-laplace(%s)", baseline))
}

#' Piecewise-exponential (Poisson-trick) Bayesian PH model -> psurv.
#'
#' A genuine Bayesian proportional-hazards analogue with a piecewise-constant
#' baseline hazard, fit WITHOUT Stan: split the follow-up at cut points into
#' person-time pieces, fit a Poisson GLM with piece-specific intercepts (the log
#' baseline hazards) and covariate log-HRs, and draw the posterior from the
#' GLM's Laplace approximation N(coef, vcov). This is the plan's robustness
#' comparator (B1) and shows the framework is model-agnostic — it produces the
#' identical `psurv` object as stan_surv or the Cox-Laplace fallback, so B2-B6
#' are unchanged. Needs only the `survival` package.
#'
#' @param formula Surv(time, status) ~ covariates.
#' @param data    the data.frame.
#' @param horizons `numeric[H]` evaluation horizons.
#' @param cuts    interior cut points; NULL => event-time quantiles giving
#'                `n_pieces` pieces.
#' @param n_pieces number of baseline-hazard pieces when `cuts` is NULL.
#' @param S,seed  posterior draws and RNG seed.
#' @export
psurv_from_pwexp <- function(formula, data, horizons, cuts = NULL,
                             n_pieces = 8L, S = 500L, seed = 1L) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("psurv_from_pwexp() needs the 'survival' package.")
  set.seed(seed)

  mf <- stats::model.frame(formula, data)
  y  <- stats::model.response(mf)
  time <- as.numeric(y[, 1L]); status <- as.integer(y[, ncol(y)])
  X <- stats::model.matrix(formula, data)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  cn <- colnames(X); n <- nrow(X)

  if (is.null(cuts)) {
    ev <- time[status == 1]
    if (!length(ev)) ev <- time
    cuts <- unique(as.numeric(stats::quantile(
      ev, seq_len(n_pieces - 1L) / n_pieces, names = FALSE)))
    cuts <- cuts[cuts > 0 & is.finite(cuts)]
  }
  edges <- c(0, cuts, Inf)
  K <- length(cuts) + 1L

  # Expand to person-time pieces (survSplit), then Poisson-trick GLM.
  dd <- data.frame(ytime = time, ystatus = status, X, check.names = FALSE)
  # survSplit requires the formula LHS to be a bare `Surv(...)` call (it rejects
  # `survival::Surv`), so build the formula in an env where Surv resolves.
  fenv <- new.env(parent = environment())
  fenv$Surv <- survival::Surv
  fsplit <- stats::as.formula("Surv(ytime, ystatus) ~ .", env = fenv)
  sp <- survival::survSplit(fsplit, data = dd,
                            cut = cuts, episode = ".piece", start = ".t0")
  sp$.expo <- sp$ytime - sp$.t0
  sp$.d <- sp$ystatus
  sp <- sp[sp$.expo > 0, , drop = FALSE]
  sp$.piece <- factor(sp$.piece, levels = seq_len(K))
  keep <- levels(sp$.piece)[table(sp$.piece) > 0]          # non-empty pieces
  sp$.piece <- factor(sp$.piece, levels = keep)

  form <- stats::as.formula(paste(".d ~ 0 + .piece",
    if (length(cn)) paste("+", paste(sprintf("`%s`", cn), collapse = " + "))
    else ""))
  fit <- stats::glm(form, family = stats::poisson(), data = sp,
                    offset = log(sp$.expo))
  co <- stats::coef(fit); V <- as.matrix(stats::vcov(fit))
  ok <- is.finite(co) & !is.na(diag(V))
  co <- co[ok]; V <- V[ok, ok, drop = FALSE]

  draws <- .rmvn(S, co, V)
  colnames(draws) <- names(co)
  a_names <- grep("^\\.piece", names(co), value = TRUE)
  b_names <- setdiff(names(co), a_names)
  # map present pieces back to their integer index
  piece_idx <- as.integer(sub("^\\.piece", "", a_names))

  alpha <- draws[, a_names, drop = FALSE]                  # S x K_present
  beta  <- if (length(b_names)) draws[, b_names, drop = FALSE]
           else matrix(0, S, 0)
  Xb <- if (length(b_names)) X[, sub("^`|`$", "", b_names), drop = FALSE] else X

  # width of each PRESENT piece within [0, t_h]
  widths <- vapply(horizons, function(t) {
    lo <- edges[piece_idx]; hi <- edges[piece_idx + 1L]
    pmax(0, pmin(t, hi) - lo)
  }, numeric(length(piece_idx)))                           # K_present x H
  if (is.null(dim(widths))) widths <- matrix(widths, nrow = length(piece_idx))

  eta  <- if (ncol(beta)) beta %*% t(Xb) else matrix(0, S, n)  # S x n
  surv <- array(NA_real_, dim = c(S, n, length(horizons)))
  for (s in seq_len(S)) {
    H0 <- as.numeric(exp(alpha[s, ]) %*% widths)           # length H
    surv[s, , ] <- exp(-outer(exp(eta[s, ]), H0))
  }
  xdf <- as.data.frame(X); xdf[[".time"]] <- time; xdf[[".status"]] <- status
  new_psurv(surv, eta, time, status, horizons, x = xdf,
            label = sprintf("pwexp(%d pieces)", K))
}

#' BART survival (nonparametric Bayesian) model -> psurv.
#'
#' The plan's nonparametric comparator (B1): demonstrates the framework is
#' model-agnostic by feeding a Bayesian machine-learning survival model through
#' the SAME `psurv` interface. Requires the `BART` package (Sparapani et al.);
#' guarded like the stan_surv backend so the framework still loads without it.
#'
#' @param formula Surv(time, status) ~ covariates.
#' @param data    the data.frame.
#' @param horizons `numeric[H]` evaluation horizons.
#' @param time,status response column names.
#' @param S       posterior draws to thin to.
#' @param seed    RNG seed.
#' @export
psurv_from_bart <- function(formula, data, horizons, time, status,
                            S = 500L, seed = 1L) {
  if (!requireNamespace("BART", quietly = TRUE))
    stop("psurv_from_bart() needs the 'BART' package. Use backend = 'pwexp' ",
         "or 'coxph' to run with only the 'survival' package.")
  X <- stats::model.matrix(formula, data)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  n <- nrow(X)
  ti <- data[[time]]; di <- data[[status]]

  fit <- BART::surv.bart(x.train = X, times = ti, delta = di,
                         ndpost = S, seed = seed)
  # Predict the survival matrix at the model's time grid for each subject.
  pre <- BART::surv.pre.bart(times = ti, delta = di, x.train = X, x.test = X)
  pred <- stats::predict(fit, newdata = pre$tx.test, mc.cores = 1L)
  grid <- fit$times                                   # BART's event-time grid
  Sdraws <- pred$surv.test                            # (S) x (n * |grid|)
  Sd <- nrow(Sdraws)
  surv <- array(NA_real_, dim = c(Sd, n, length(horizons)))
  for (i in seq_len(n)) {
    cols <- ((i - 1L) * length(grid) + 1L):(i * length(grid))
    Si <- Sdraws[, cols, drop = FALSE]                # Sd x |grid|
    for (h in seq_along(horizons)) {
      j <- max(1L, findInterval(horizons[h], grid))
      surv[, i, h] <- Si[, j]
    }
  }
  # BART has no single linear predictor; use -log cumulative hazard proxy at the
  # largest horizon as the risk score for rank-based metrics.
  eta <- -log(pmax(surv[, , length(horizons)], 1e-12))
  new_psurv(surv, eta, ti, di, horizons, x = as.data.frame(X),
            label = "bart-survival")
}

# --- small numerical helpers (no hard dependency on MASS) -------------------

#' Draw n rows from a multivariate normal via Cholesky (falls back to eigen).
#' @keywords internal
.rmvn <- function(n, mu, Sigma) {
  p <- length(mu)
  L <- tryCatch(chol(Sigma), error = function(e) NULL)
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  if (!is.null(L)) {
    out <- Z %*% L
  } else {
    e <- eigen(Sigma, symmetric = TRUE)
    A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), p) %*% t(e$vectors)
    out <- Z %*% A
  }
  sweep(out, 2L, mu, `+`)
}

#' One Dirichlet(1,...,1) weight vector of length k (Bayesian bootstrap).
#' @keywords internal
.dirichlet1 <- function(k) {
  g <- stats::rexp(k)
  g / sum(g)
}
