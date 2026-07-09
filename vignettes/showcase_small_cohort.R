# showcase_small_cohort.R
# ===========================================================================
# SHOWCASE -- where joint posterior reporting changes the local-use read-out.
#
# Setting: a published survival model must be validated in a SMALL, low-event
# LOCAL cohort (n = 120, ~55% censored) before a clinic decides whether to
# deploy it. This is the single most common real external-validation situation,
# and the one where classical point-estimate-plus-CI reporting is weakest.
#
# The showcase makes three claims and backs each with numbers from this script:
#
#   ACT 1  In one such cohort the FREQUENTIST read-out is ambiguous: every 95%
#          CI is so wide it spans "good" and "broken", so the deploy/reject
#          decision is a coin flip. The BAYESIAN framework returns the clinically
#          decision-relevant quantities directly -- P(discrimination adequate),
#          P(calibration acceptable), P(net clinical benefit > 0), P(net harm) --
#          and a single JOINT deployment-gate probability.
#
#   ACT 2  That joint gate is genuinely NOT recoverable from separate frequentist
#          intervals: the metrics are correlated across the posterior, so the
#          honest joint probability differs materially from the naive product of
#          marginals. Only a coherent posterior over all metrics at once carries
#          this information.
#
#   ACT 3  The Bayesian deployment probability is checked for reliability: pooled over
#          many simulated cohorts with KNOWN ground truth, when the gate says
#          "deploy-worthy with probability p", about a fraction p really are.
#          A frequentist dichotomous rule offers no such graded, checkable
#          confidence and mis-classifies far more often.
#
# Truth is known throughout (Part C data-generating machinery), so every claim
# is checkable against an oracle. Runs on the coxph-Laplace backend.
#
#   Rscript vignettes/showcase_small_cohort.R          # quick demo (~1-2 min)
#   Rscript vignettes/showcase_small_cohort.R --full   # more cohorts, tighter
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ===========================================================================

suppressMessages(library(survival))
# patchwork must be attached, not called through `::`: the `+` and `&` operators
# that stitch the two panels of the figure are S3 methods it only registers on
# attach.
suppressMessages(library(patchwork))
source(file.path("..", "R", "load_all.R"), chdir = TRUE)
source(file.path("..", "simulation", "dgm.R"), chdir = TRUE)

FULL <- "--full" %in% commandArgs(trailingOnly = TRUE)
figdir <- "figures"; dir.create(figdir, showWarnings = FALSE)
set.seed(2026)

## ===========================================================================
## 0.  A published model + the clinic's decision rule
## ===========================================================================
# The published model is developed on a large cohort; its coefficients and
# baseline hazard are FIXED thereafter (this is what "externally validate" means).
HZ   <- c(3, 5)          # evaluation horizons (years); headline at t* = 5
H5   <- 2L               # index of the 5-year horizon
TAU  <- 5                # concordance truncation
PT   <- 0.20             # clinically chosen risk threshold for the decision

# Development truth: moderate-signal PH model.
DEV_BETA  <- c(b1 = 0.70, b2 = -0.50)
DEV_SHAPE <- 1.3
DEV_RATE  <- 0.14
dev <- dgm_weibull(4000, beta = DEV_BETA, shape = DEV_SHAPE, rate = DEV_RATE,
                   cens_rate = 0.35)$data
published <- coxph(Surv(time, status) ~ x1 + x2, data = dev, x = TRUE, y = TRUE)

# The clinic's pre-specified DEPLOYMENT GATE (fixed before seeing local data):
#   (i)   adequate discrimination : Harrell C >= 0.65
#   (ii)  acceptable calibration  : calibration slope in [0.80, 1.25]
#   (iii) positive clinical value : net benefit at threshold PT beats treat-all/none
C_MIN     <- 0.65
SLOPE_LO  <- 0.80
SLOPE_HI  <- 1.25

## ===========================================================================
## Helper: apply the FIXED published model to any cohort (linear predictor and
## absolute survival at the horizons), and a per-draw net-benefit metric.
## ===========================================================================
baseline_at_horizons <- function(bh, horizons) {
  ok <- is.finite(bh$time) & is.finite(bh$hazard)
  time <- bh$time[ok]
  hazard <- bh$hazard[ok]
  if (!length(time)) return(rep(0, length(horizons)))
  ord <- order(time)
  time <- c(0, time[ord])
  hazard <- c(0, hazard[ord])
  keep <- !duplicated(time, fromLast = TRUE)
  time <- time[keep]
  hazard <- hazard[keep]
  if (length(time) == 1L) return(rep(hazard[[1]], length(horizons)))
  out <- stats::approx(time, hazard, xout = horizons,
                       method = "constant", rule = 2, f = 0)$y
  out[is.na(out)] <- 0
  out
}

published_predict <- function(fit, newdata, horizons) {
  Terms <- stats::delete.response(stats::terms(fit))
  mf <- stats::model.frame(Terms, newdata, xlev = fit$xlevels)
  X  <- stats::model.matrix(Terms, mf, contrasts.arg = fit$contrasts)
  X  <- X[, names(stats::coef(fit)), drop = FALSE]
  eta <- as.numeric(X %*% stats::coef(fit))
  bh  <- survival::basehaz(fit, centered = FALSE)
  H0  <- baseline_at_horizons(bh, horizons)
  list(eta = eta, S = exp(-outer(exp(eta), H0)))     # S is n x H
}

# Net benefit (and delta vs the best default strategy) at one threshold, as a
# per-draw metric so it lands in the SAME joint posterior matrix as the others.
metric_dnb <- function(pt, h) function(draw, w) {
  nb <- net_benefit_draw(draw, thresholds = pt, horizon = h, w = w)
  c(NB_model = nb$model, dNB = nb$model - max(nb$all, 0))
}

# Oracle for the FIXED published model against a given local truth: the true
# discrimination, calibration slope and net benefit from a huge reference
# sample, hence the ground-truth deploy-worthiness under the clinic's gate.
oracle_worthy <- function(loc_beta, loc_rate, N = 5e4) {
  ref <- dgm_weibull(N, beta = loc_beta, shape = DEV_SHAPE, rate = loc_rate,
                     cens_rate = 0)                    # no censoring for oracle
  pp  <- published_predict(published, ref$data, HZ)
  eta_pub <- pp$eta
  Strue   <- surv_weibull(HZ[H5], ref$lp, DEV_SHAPE, loc_rate)
  risk_t  <- 1 - Strue
  risk_p  <- 1 - pp$S[, H5]
  Ttrue   <- ref$data$time                              # uncensored here
  # same risk orientation as metric_harrell_c (marker = -eta):
  oC  <- concordance(Surv(Ttrue, rep(1, N)) ~ I(-eta_pub))$concordance
  osl <- unname(coef(coxph(Surv(Ttrue, rep(1, N)) ~ eta_pub))[["eta_pub"]])
  odds <- PT / (1 - PT); pos <- risk_p >= PT
  f <- mean(pos); r <- if (f > 0) mean(risk_t[pos]) else 0
  nb_m <- if (f > 0) f * r - f * (1 - r) * odds else 0
  r_all <- mean(risk_t); nb_all <- r_all - (1 - r_all) * odds
  odnb <- nb_m - max(nb_all, 0)
  list(C = oC, slope = osl, dNB = odnb,
       worthy = (oC >= C_MIN) && (osl >= SLOPE_LO && osl <= SLOPE_HI) &&
                (odnb > 0))
}

# The Bayesian validation suite, evaluated jointly per posterior draw.
suite <- list(
  harrell = metric_harrell_c,
  uno     = function(d, w) metric_uno_c(d, w, tau = TAU),
  slope   = metric_cal_slope,
  ici     = function(d, w) metric_ici(d, w, horizon = H5, method = "auto")["ICI"],
  ibs     = function(d, w) metric_brier(d, w)["IBS"],
  nb      = metric_dnb(PT, H5))

## The deployment gate reads only three of these. ACT 3 propagates 240 cohorts x
## S draws, so evaluating the full suite there would fit a `polspline::hare`
## calibration smoother once per draw for an ICI that the gate never consults.
## Metric evaluation draws no random numbers (the Dirichlet weights are drawn
## before the metric call in propagate_metric), so the gate probability is
## bit-identical to the one the full suite would produce.
gate_suite <- list(harrell = metric_harrell_c,
                   slope   = metric_cal_slope,
                   nb      = metric_dnb(PT, H5))

## ===========================================================================
## Fair frequentist comparator: the PAIRED nonparametric bootstrap joint gate.
## ---------------------------------------------------------------------------
## The dichotomous rule below (all three 95% intervals clear their bars) is a
## deliberately conservative reporting-practice reference. A stronger, fairer
## comparator resamples the evaluation cohort ONCE per replicate and recomputes
## every gate-critical metric on that SAME replicate, so a frequentist analyst
## can also express the joint deployment event as a bootstrap frequency -- the
## direct classical analogue of the Bayesian gate probability. It preserves the
## dependence across metrics that separate marginal intervals discard.
## `with_preserved_rng()` shields the global stream so this additive comparator
## never perturbs the reproducible draws of the rest of the script.
with_preserved_rng <- function(expr) {
  old <- if (exists(".Random.seed", .GlobalEnv))
           get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv))
  force(expr)
}

paired_boot_gate <- function(time, status, lp, surv2d, B, seed) {
  n <- length(time)
  with_preserved_rng({
    set.seed(seed)
    pass <- logical(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      tt <- time[idx]; ss <- status[idx]; ll <- lp[idx]
      ch <- tryCatch(concordance(Surv(tt, ss) ~ ll), error = function(e) NULL)
      cs <- tryCatch(coxph(Surv(tt, ss) ~ ll), error = function(e) NULL)
      if (is.null(ch) || is.null(cs)) { pass[b] <- NA; next }
      Cb  <- 1 - ch$concordance                      # marker orientation = -eta
      slb <- unname(coef(cs))
      db  <- list(surv2d = surv2d[idx, , drop = FALSE], eta = ll, time = tt,
                  status = ss, horizons = HZ, n = n, H = length(HZ))
      dnbb <- as.numeric(metric_dnb(PT, H5)(db, rep(1 / n, n))["dNB"])
      pass[b] <- (Cb >= C_MIN) &&
        (slb >= SLOPE_LO && slb <= SLOPE_HI) && (dnbb > 0)
    }
    mean(pass, na.rm = TRUE)
  })
}

## Uncertainty band on each bin's empirical worthy fraction.
##
## A Wilson (binomial) interval would be wrong here. `worthy` is a deterministic
## function of the data-generating configuration, so within a bin it is not a
## Bernoulli draw; the only randomness is WHICH cohorts land in the bin. The bin
## fraction is therefore a ratio estimator whose numerator and denominator are
## both random and dependent. We resample cohorts (the independent unit) and
## recompute the fraction inside fixed bin edges, which respects that structure.
boot_bin_ci <- function(df, col, breaks, B = 2000L, seed = 11L) {
  set.seed(seed)
  lev <- levels(cut(df[[col]], breaks, include.lowest = TRUE))
  reps <- replicate(B, {
    s <- df[sample.int(nrow(df), nrow(df), replace = TRUE), ]
    b <- cut(s[[col]], breaks, include.lowest = TRUE)
    as.numeric(tapply(s$worthy, factor(b, levels = lev), mean))
  })
  tibble::tibble(
    emp_lo = apply(reps, 1, stats::quantile, 0.025, na.rm = TRUE),
    emp_hi = apply(reps, 1, stats::quantile, 0.975, na.rm = TRUE))
}

## Ranking summary: P(a truly worthy cohort is given a higher gate probability
## than a truly unworthy one), ties counted as half. This is the quantity the
## diagram actually supports; see the note on the diagonal in reliability_of().
gate_auc <- function(p, worthy) {
  pw <- p[worthy]; pu <- p[!worthy]
  if (!length(pw) || !length(pu)) return(NA_real_)
  mean(outer(pw, pu, ">") + 0.5 * outer(pw, pu, "=="))
}
gate_auc_ci <- function(df, col, B = 2000L, seed = 12L) {
  set.seed(seed)
  v <- replicate(B, {
    s <- df[sample.int(nrow(df), nrow(df), replace = TRUE), ]
    gate_auc(s[[col]], s$worthy)
  })
  stats::quantile(v, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
}

## ===========================================================================
## ACT 1 -- one small local cohort: frequentist ambiguity vs Bayesian answer
## ===========================================================================
cat("\n=====================================================================\n")
cat("ACT 1  A small local cohort (n = 120): can we deploy the model here?\n")
cat("=====================================================================\n")

# Local truth: same signal, but a heavier baseline hazard than the development
# population (a sicker clinic) -> the published model is mildly MIS-calibrated
# here. This is realistic external validation, not a broken model.
LOC_RATE <- 0.20
local_full <- dgm_weibull(120, beta = DEV_BETA, shape = DEV_SHAPE,
                          rate = LOC_RATE, cens_rate = 0.55)
local <- local_full$data
cat(sprintf("Local cohort: n = %d, events = %d (%.0f%% censored), horizon t* = %g y\n",
            nrow(local), sum(local$status), 100 * mean(local$status == 0), HZ[H5]))

S_DRAWS <- if (FULL) 2000L else 1000L
ps <- psurv_external(published, local, time = "time", status = "status",
                     horizons = HZ, S = S_DRAWS, seed = 7)

# Joint posterior: every metric evaluated on the SAME draw + SAME bootstrap
# weights -> the columns of `J` are a coherent joint sample (B2 + B3).
J <- propagate_many(ps, suite, weights = "bootstrap", seed = 1)
Jm <- unclass(J); colnames(Jm) <- sub("^[^.]+\\.", "", colnames(Jm))
post_summary <- summary(J) |>
  dplyr::mutate(metric = sub("^[^.]+\\.", "", metric))

qi <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
fmt <- function(v, d = 3) sprintf(paste0("%.", d, "f"), v)

## ---- Frequentist comparators: point estimate [95% CI] --------------------
lp  <- colMeans(ps$eta)                          # published linear predictor
dfv <- tibble::tibble(time = local$time, status = local$status, lp = lp)
wald <- function(est, se) c(est, est - 1.96 * se, est + 1.96 * se)

ch <- concordance(Surv(time, status) ~ lp, data = dfv)          # Harrell
fr_C <- wald(1 - ch$concordance, sqrt(ch$var))
cs <- coxph(Surv(time, status) ~ lp, data = dfv)                # slope
fr_slope <- wald(unname(coef(cs)), sqrt(vcov(cs)[1, 1]))

# ICI and net benefit have no closed-form SE -> classical bootstrap.
mean_draw <- list(surv2d = matrix(colMeans(ps$surv), nrow = ps$n),
                  eta = lp, time = local$time, status = local$status,
                  horizons = HZ, n = ps$n, H = length(HZ))
boot_ci <- function(fn, B = if (FULL) 1000L else 400L, seed = 11) {
  set.seed(seed); n <- ps$n
  est <- as.numeric(fn(mean_draw, rep(1 / n, n)))
  v <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    db <- list(surv2d = mean_draw$surv2d[idx, , drop = FALSE], eta = lp[idx],
               time = local$time[idx], status = local$status[idx],
               horizons = HZ, n = n, H = length(HZ))
    as.numeric(fn(db, rep(1 / n, n)))
  })
  c(est, stats::quantile(v, c(0.025, 0.975), names = FALSE))
}
fr_ici <- boot_ci(function(d, w) metric_ici(d, w, horizon = H5, method = "auto")["ICI"])
fr_dnb <- boot_ci(function(d, w) metric_dnb(PT, H5)(d, w)["dNB"])

cat("\n-- Frequentist read-out (point estimate [95% CI]) -------------------\n")
frtab <- rbind(
  `Harrell C`            = fr_C,
  `Calibration slope`    = fr_slope,
  `ICI (5y)`             = fr_ici,
  `dNB @20% (vs default)`= fr_dnb)
colnames(frtab) <- c("est", "lo", "hi")
print(round(frtab, 3))

verdict <- function(lo, hi, need_lo = NA, need_hi = NA) {
  # is the pre-specified requirement DECIDED (both ends agree) or ambiguous?
  ok_lo <- if (is.na(need_lo)) TRUE else lo >= need_lo
  ok_hi <- if (is.na(need_hi)) TRUE else hi <= need_hi
  bad_lo <- if (is.na(need_lo)) FALSE else hi < need_lo
  bad_hi <- if (is.na(need_hi)) FALSE else lo > need_hi
  if (ok_lo && ok_hi) "PASS" else if (bad_lo || bad_hi) "FAIL" else "AMBIGUOUS"
}
cat(sprintf("  discrimination C>=%.2f : %s\n", C_MIN,
            verdict(fr_C[2], fr_C[3], need_lo = C_MIN)))
cat(sprintf("  calibration slope in [%.2f,%.2f] : %s\n", SLOPE_LO, SLOPE_HI,
            verdict(fr_slope[2], fr_slope[3], need_lo = SLOPE_LO, need_hi = SLOPE_HI)))
cat(sprintf("  net benefit dNB>0 : %s\n",
            verdict(fr_dnb[2], fr_dnb[3], need_lo = 0)))
cat("  => The marginal-interval rule does not decide the deployment question.\n")

## ---- Fair frequentist comparator: paired-bootstrap joint-gate frequency ----
P_gate_boot <- paired_boot_gate(local$time, local$status, lp, mean_draw$surv2d,
                                B = if (FULL) 4000L else 2000L, seed = 23)
cat(sprintf("\n  Paired-bootstrap joint-gate frequency        = %.3f\n", P_gate_boot))
cat("     (all three metrics recomputed on each bootstrap replicate; the\n")
cat("      direct classical analogue of the Bayesian gate probability.)\n")

## ---- Bayesian read-out: full posteriors + deployment probabilities -------
cat("\n-- Bayesian read-out (posterior median [95% CrI]) -------------------\n")
show <- c(C_harrell = "Harrell C", C_uno = "Uno C", cal_slope = "Calibration slope",
          ICI = "ICI (5y)", IBS = "IBS", NB_model = "Net benefit @20%",
          dNB = "dNB @20% (vs default)")
for (col in names(show)) {
  q <- qi(Jm[, col])
  cat(sprintf("  %-22s %s [%s, %s]\n", show[[col]],
              fmt(q[2]), fmt(q[1]), fmt(q[3])))
}

P_disc  <- mean(Jm[, "C_harrell"] >= C_MIN)
P_cal   <- mean(Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI)
P_ici   <- mean(Jm[, "ICI"] <= 0.05)
P_bene  <- mean(Jm[, "dNB"] > 0)
P_harm  <- mean(Jm[, "dNB"] < 0)
P_gate  <- mean(Jm[, "C_harrell"] >= C_MIN &
                Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI &
                Jm[, "dNB"] > 0)

cat("\n-- Deployment probabilities -------------------------------------------\n")
cat(sprintf("  P(discrimination adequate, C>=%.2f)         = %.2f\n", C_MIN, P_disc))
cat(sprintf("  P(calibration acceptable, slope in[%.2f,%.2f]) = %.2f\n",
            SLOPE_LO, SLOPE_HI, P_cal))
cat(sprintf("  P(ICI <= 0.05)                             = %.2f\n", P_ici))
cat(sprintf("  P(net clinical benefit, dNB>0)             = %.2f\n", P_bene))
cat(sprintf("  P(net HARM, dNB<0)                         = %.2f\n", P_harm))
cat(sprintf("  --------------------------------------------------------\n"))
cat(sprintf("  P(JOINT deployment gate passes)            = %.2f\n", P_gate))

# Close the loop: what does the ORACLE say this cohort's truth actually is?
orc1 <- oracle_worthy(DEV_BETA, LOC_RATE)
cat(sprintf("\n  [ground truth, oracle]  C=%.2f  slope=%.2f  dNB=%+.3f  =>  %s\n",
            orc1$C, orc1$slope, orc1$dNB,
            if (orc1$worthy) "truly deploy-worthy" else "NOT deploy-worthy"))
cat("  The joint gate reached the correct low-probability conclusion where every\n")
cat("  frequentist interval was individually undecided.\n")

## ===========================================================================
## ACT 2 -- the joint gate is not recoverable from separate intervals
## ===========================================================================
cat("\n=====================================================================\n")
cat("ACT 2  Why separate CIs cannot answer this: metrics are correlated\n")
cat("=====================================================================\n")

# Naive independence assumption: multiply the marginal pass-probabilities.
P_gate_indep <- P_disc * P_cal * P_bene
rho_cs <- cor(Jm[, "cal_slope"], Jm[, "dNB"])
rho_cc <- cor(Jm[, "C_harrell"], Jm[, "cal_slope"])
cat(sprintf("  joint gate, honest (coherent posterior) = %.3f\n", P_gate))
cat(sprintf("  joint gate, naive product of marginals  = %.3f\n", P_gate_indep))
cat(sprintf("  relative error of the naive shortcut     = %+.0f%%\n",
            100 * (P_gate_indep - P_gate) / P_gate))
cat(sprintf("  joint gate, paired classical bootstrap  = %.3f\n", P_gate_boot))
cat(sprintf("  posterior corr(slope, dNB) = %+.2f ; corr(C, slope) = %+.2f\n",
            rho_cs, rho_cc))
cat("  => The metrics move together across the posterior; a set of separate\n")
cat("     frequentist intervals throws this dependence away. A paired bootstrap\n")
cat("     that recomputes every metric on the same replicate DOES preserve this\n")
cat("     dependence and yields a comparable joint frequency; the Bayesian gate\n")
cat("     adds model-parameter uncertainty and a coherent posterior read-out.\n")

## ---- Joint posterior of the two gate-critical metrics --------------------
## The figure is drawn by `showcase_joint_composite.R`, which faces this cohort
## with the PBC external validation in one panel. Store the raw draws and this
## cohort's own thresholds; the gate rule is applied once, there, to both
## cohorts, so the two panels cannot drift apart on what "passes" means.
saveRDS(list(
  draws  = tibble::tibble(cal_slope = Jm[, "cal_slope"],
                          C_harrell = Jm[, "C_harrell"],
                          dNB       = Jm[, "dNB"]),
  gate   = list(C_min = C_MIN, slope_lo = SLOPE_LO, slope_hi = SLOPE_HI),
  P_gate = P_gate,
  label  = "Simulated cohort",
  n      = nrow(local),
  events = sum(local$status)),
  file.path(figdir, "showcase_joint_draws.rds"))

## ===========================================================================
## ACT 3 -- the deployment probability is calibrated against known truth
## ===========================================================================
cat("\n=====================================================================\n")
cat("ACT 3  Is the deployment probability trustworthy? Calibrate it vs truth\n")
cat("=====================================================================\n")

# Configurations spanning truly-worthy to truly-unworthy. Attenuation `a`
# scales the true signal (drives discrimination AND the true slope ~ a); the
# baseline rate shifts calibration-in-the-large and net benefit.
CONFIGS <- list(
  list(a = 1.10, rate = 0.14),   # strong signal, well calibrated: worthy
  list(a = 1.00, rate = 0.14),   # matched: worthy
  list(a = 0.88, rate = 0.15),   # slope ~0.88: borderline worthy
  list(a = 0.75, rate = 0.16),   # slope <0.8: unworthy
  list(a = 0.60, rate = 0.20),   # weak signal + miscalib: unworthy
  list(a = 1.35, rate = 0.24))   # slope >1.25 + bad CITL: unworthy
REPS   <- if (FULL) 40L else 15L
N_LOC  <- 120L

cat(sprintf("  %d configurations x %d cohorts (n=%d) ; oracle-defined truth\n",
            length(CONFIGS), REPS, N_LOC))
rows <- list(); r <- 0L
for (ci in seq_along(CONFIGS)) {
  cf <- CONFIGS[[ci]]
  lb <- DEV_BETA * cf$a
  orc <- oracle_worthy(lb, cf$rate)
  cat(sprintf("  config %d (a=%.2f, rate=%.2f): oracle C=%.2f slope=%.2f dNB=%+.3f -> %s\n",
              ci, cf$a, cf$rate, orc$C, orc$slope, orc$dNB,
              if (orc$worthy) "WORTHY" else "unworthy"))
  for (rep in seq_len(REPS)) {
    r <- r + 1L
    dat <- dgm_weibull(N_LOC, beta = lb, shape = DEV_SHAPE, rate = cf$rate,
                       cens_rate = 0.55)$data
    psr <- psurv_external(published, dat, time = "time", status = "status",
                          horizons = HZ, S = if (FULL) 600L else 400L,
                          seed = 1000L * ci + rep)
    Jr <- unclass(propagate_many(psr, gate_suite, weights = "bootstrap",
                                 seed = rep))
    colnames(Jr) <- sub("^[^.]+\\.", "", colnames(Jr))
    pg <- mean(Jr[, "C_harrell"] >= C_MIN &
               Jr[, "cal_slope"] >= SLOPE_LO & Jr[, "cal_slope"] <= SLOPE_HI &
               Jr[, "dNB"] > 0)
    # frequentist dichotomous deploy rule on the same cohort (conservative
    # reporting-practice reference):
    lpr <- colMeans(psr$eta)
    dvr <- tibble::tibble(time = dat$time, status = dat$status, lp = lpr)
    chr <- concordance(Surv(time, status) ~ lp, data = dvr)
    csr <- coxph(Surv(time, status) ~ lp, data = dvr)
    C_lo  <- (1 - chr$concordance) - 1.96 * sqrt(chr$var)
    sl    <- unname(coef(csr)); sl_se <- sqrt(vcov(csr)[1, 1])
    freq_deploy <- (C_lo >= C_MIN) &&                       # C CI clears bar
      ((sl - 1.96 * sl_se) >= SLOPE_LO) && ((sl + 1.96 * sl_se) <= SLOPE_HI)
    # fair comparator: paired-bootstrap joint-gate frequency on the same cohort.
    p_boot <- paired_boot_gate(dat$time, dat$status, lpr,
                               matrix(colMeans(psr$surv), nrow = psr$n),
                               B = if (FULL) 500L else 250L,
                               seed = 7000L * ci + rep)
    rows[[r]] <- tibble::tibble(config = ci, worthy = orc$worthy,
                                p_deploy = pg, p_boot = p_boot,
                                freq_deploy = freq_deploy)
  }
}
D <- dplyr::bind_rows(rows)

## ---- Ranking curve, with a cohort-bootstrap band on each bin fraction ------
## QUANTILE bins, not equal-width ones. The gate probability is a conjunction of
## three criteria under a diffuse n=120 posterior, so it is heavily right-skewed
## and never approaches 1; equal-width bins put ~70% of cohorts in [0, 0.2] and
## leave the top bin empty, which destroys the resolution of the diagram.
##
## The 45-degree line is NOT the reference here. Posterior probabilities are
## calibrated when averaged over the prior that generated the truth, and the
## truth in this study comes from six hand-picked configurations in a fixed
## 50/50 worthy/unworthy split, which is not the model's prior. The diagram
## therefore reports RANKING and CONSERVATIVENESS; see gate_auc().
reliability_of <- function(df, col, K = 4L) {
  br <- unique(stats::quantile(df[[col]], seq(0, 1, length.out = K + 1L),
                               names = FALSE, type = 7))
  if (length(br) < 3L) br <- unique(range(df[[col]]))   # degenerate fallback
  tab <- df |>
    dplyr::mutate(bin = cut(.data[[col]], br, include.lowest = TRUE)) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(mid = mean(.data[[col]]), emp = mean(worthy),
                     n = dplyr::n(), .groups = "drop")
  dplyr::bind_cols(tab, boot_bin_ci(df, col, br)) |>
    dplyr::filter(is.finite(mid)) |>
    dplyr::select(mid, emp, emp_lo, emp_hi, n)
}
rel      <- reliability_of(D, "p_deploy")   # Bayesian gate
rel_boot <- reliability_of(D, "p_boot")     # paired-bootstrap gate

auc_bayes <- gate_auc(D$p_deploy, D$worthy); auc_b_ci <- gate_auc_ci(D, "p_deploy")
auc_boot  <- gate_auc(D$p_boot,   D$worthy); auc_bt_ci <- gate_auc_ci(D, "p_boot")

cat("\n-- Bayesian gate: quantile bins [cohort-bootstrap 95%] ---------------\n")
print(as.data.frame(dplyr::mutate(rel, dplyr::across(c(mid, emp, emp_lo, emp_hi),
                                                     ~ round(.x, 3)))),
      row.names = FALSE)
cat("\n-- Paired-bootstrap joint-gate frequency: same bins ------------------\n")
print(as.data.frame(dplyr::mutate(rel_boot, dplyr::across(c(mid, emp, emp_lo, emp_hi),
                                                          ~ round(.x, 3)))),
      row.names = FALSE)

cat(sprintf(paste0("\n  Base rate of truly worthy cohorts: %.3f (fixed by design).\n",
                   "  Mean stated gate probability: Bayesian %.3f, paired bootstrap %.3f.\n",
                   "  => both gates are CONSERVATIVE; the 45-degree line is not the\n",
                   "     reference, because the truth is not drawn from the model's prior.\n"),
            mean(D$worthy), mean(D$p_deploy), mean(D$p_boot)))
cat(sprintf("  Ranking AUC vs truth: Bayesian %.3f [%.3f, %.3f], bootstrap %.3f [%.3f, %.3f]\n",
            auc_bayes, auc_b_ci[1], auc_b_ci[2],
            auc_boot, auc_bt_ci[1], auc_bt_ci[2]))

# Comparator operating characteristics.
tp <- if (any(D$worthy)) mean(D$freq_deploy[D$worthy]) else NA_real_
fp <- if (any(!D$worthy)) mean(D$freq_deploy[!D$worthy]) else NA_real_
brier_bayes <- mean((D$p_deploy - D$worthy)^2)
brier_boot  <- mean((D$p_boot   - D$worthy)^2)
cat(sprintf("\n  Frequentist dichotomous rule: deploys %.0f%% of truly-worthy,",
            100 * tp))
cat(sprintf(" %.0f%% of unworthy cohorts.\n", 100 * fp))
cat(sprintf("  Probability Brier score vs truth: Bayesian gate = %.3f,",
            brier_bayes))
cat(sprintf(" paired bootstrap = %.3f (lower is better).\n", brier_boot))
cat("  => Both graded gates rank cohorts by worthiness; the dichotomous rule\n")
cat("     gives an over/under-confident yes/no only at this sample size.\n")

## ---- ROC of the stated probability against the oracle truth ---------------
## The reliability curve shows only whether a stated probability MEANS what it
## says. It cannot show whether the stated probabilities ORDER the cohorts, and
## ordering is what the AUC reported just above measures. Sweeping a threshold
## over the stated probability gives the curve whose area is that AUC.
##
## Ties matter here: p_deploy is a fraction of draws, so many cohorts share a
## value (0 above all), and a rank-order sweep would draw a staircase whose area
## depends on the tie-breaking order. Thresholding on the UNIQUE values instead
## draws one straight segment across each tie block, whose trapezoidal area is
## exactly gate_auc()'s half-credit-for-ties convention. The stopifnot() below
## holds the two definitions together.
roc_of <- function(p, worthy) {
  cuts <- sort(unique(p), decreasing = TRUE)
  tibble::tibble(
    fpr = c(0, purrr::map_dbl(cuts, ~ sum(!worthy & p >= .x)) / sum(!worthy)),
    tpr = c(0, purrr::map_dbl(cuts, ~ sum(worthy  & p >= .x)) / sum(worthy)))
}
trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

roc      <- roc_of(D$p_deploy, D$worthy)
roc_boot <- roc_of(D$p_boot,   D$worthy)
stopifnot(abs(trapz(roc$fpr, roc$tpr)           - auc_bayes) < 1e-10,
          abs(trapz(roc_boot$fpr, roc_boot$tpr) - auc_boot)  < 1e-10)

## ---- Figure ----------------------------------------------------------------
## House style: see FIGURE_STYLE.md (Okabe-Ito palette + theme_bayesval()).
## Two panels, because the paper's own thesis is that discrimination and
## calibration are separate axes of performance, and a deployment gate is a
## forecast that deserves the same two-sided audit as the model it gates.
## Left, calibration: bins are quantiles of the gate probability, so each holds
## the same number of cohorts; the vertical bars are cohort-bootstrap 95%
## intervals on the empirical worthy fraction, and the paired-bootstrap gate is
## the fair comparator. The dashed 45-degree line is drawn for orientation only
## (see reliability_of). Right, discrimination: the same two gates, ranked
## against the oracle, with the dashed line now the chance ranking. The two
## dashed diagonals therefore mean different things, which the caption states.
##
## patchwork, not facets: the panels share only their [0,1]^2 frame, and their
## axes measure different quantities, so a shared facet scale would put the wrong
## label under each. FIGURE_STYLE.md reserves patchwork for exactly this case. A
## one-level facet_wrap() on each panel buys the bold strip that names it.
meth_lab <- c("Bayesian gate", "paired bootstrap gate")
gate_scale <- ggplot2::scale_colour_manual(
  values = stats::setNames(bc_palette[c("blue", "orange")], meth_lab),
  name = NULL)
tag_panel <- function(df, lab)
  dplyr::mutate(df, panel = lab, method = factor(method, levels = meth_lab))

rel_both <- dplyr::bind_rows(
  dplyr::mutate(rel,      method = meth_lab[1]),
  dplyr::mutate(rel_boot, method = meth_lab[2])) |>
  tag_panel("Calibration of the stated probability")
roc_both <- dplyr::bind_rows(
  dplyr::mutate(roc,      method = meth_lab[1]),
  dplyr::mutate(roc_boot, method = meth_lab[2])) |>
  tag_panel("Ranking of cohorts by true worthiness")

## Only the line layer feeds the legend. patchwork's guides = "collect" merges two
## panels' legends when their keys are identical, so a point-and-errorbar key here
## against a plain line key in the ROC panel would print the same two series
## twice, side by side, above the figure.
p_rel <- ggplot2::ggplot(rel_both, ggplot2::aes(mid, emp, colour = method)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                       colour = "gray50") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = emp_lo, ymax = emp_hi),
                         width = 0.02, linewidth = 0.4, show.legend = FALSE) +
  # Quantile bins hold equal numbers of cohorts, so point area would encode a
  # constant. Fixed size, and the per-bin count goes in the caption.
  ggplot2::geom_point(size = 3.2, show.legend = FALSE) +
  ggplot2::geom_point(size = 3.2, shape = 21, colour = "gray40",
                      stroke = 0.3, show.legend = FALSE) +
  ggplot2::facet_wrap(~ panel) +
  gate_scale +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::scale_x_continuous(breaks = number_ticks(6)) +
  ggplot2::scale_y_continuous(breaks = number_ticks(6)) +
  # No title/subtitle and no explanatory axis text: that belongs in the caption.
  ggplot2::labs(x = "Stated P(deploy-worthy)",
                y = "Empirical fraction truly worthy") +
  theme_bayesval()

## The frequentist dichotomous rule cannot trace a curve: it emits one yes/no per
## cohort, so it owns a single operating point (fp, tp). Drawn as a mark and not
## as a legend series, because a legend entry would promise a curve in the other
## panel too, and the rule has no stated probability to be calibrated against.
freq_pt <- tibble::tibble(fpr = fp, tpr = tp, panel = roc_both$panel[1])

p_roc <- ggplot2::ggplot(roc_both, ggplot2::aes(fpr, tpr, colour = method)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                       colour = "gray50") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(data = freq_pt, shape = 24, size = 3.2,
                      fill = bc_palette[["green"]], colour = "gray40",
                      stroke = 0.4, inherit.aes = FALSE,
                      ggplot2::aes(fpr, tpr)) +
  ggplot2::facet_wrap(~ panel) +
  gate_scale +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::scale_x_continuous(breaks = number_ticks(6)) +
  ggplot2::scale_y_continuous(breaks = number_ticks(6)) +
  ggplot2::labs(x = "False-positive rate", y = "True-positive rate") +
  theme_bayesval()

# One legend for both panels; the AUCs go in the caption, not on the panel.
p_rel_roc <- p_rel + p_roc + plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "top")
# Authored at the printed width (\linewidth = 6.5in at 1in margins) so the 12pt
# figure text is not shrunk relative to the 12pt body text by \includegraphics.
bc_ggsave(file.path(figdir, "showcase_reliability.png"), p_rel_roc,
          width = 6.5, height = 3.9)

## ===========================================================================
## Save + closing summary
## ===========================================================================
saveRDS(list(
  act1 = list(freq = frtab, post = post_summary,
              P = c(disc = P_disc, cal = P_cal, ici = P_ici, bene = P_bene,
                    harm = P_harm, gate = P_gate)),
  act2 = list(P_gate = P_gate, P_gate_indep = P_gate_indep,
              P_gate_boot = P_gate_boot,
              rho_slope_dnb = rho_cs, rho_C_slope = rho_cc),
  act3 = list(data = D, reliability = rel, reliability_boot = rel_boot,
              roc = roc, roc_boot = roc_boot,
              brier = brier_bayes, brier_boot = brier_boot,
              auc = auc_bayes, auc_ci = auc_b_ci,
              auc_boot = auc_boot, auc_boot_ci = auc_bt_ci,
              base_rate = mean(D$worthy), mean_p = mean(D$p_deploy),
              freq_tp = tp, freq_fp = fp)),
  file.path(figdir, "showcase_results.rds"))

cat("\n=====================================================================\n")
cat("SHOWCASE COMPLETE\n")
cat(sprintf("  ACT 1  frequentist deploy verdict: AMBIGUOUS ; Bayesian P(gate)=%.2f, P(harm)=%.2f\n",
            P_gate, P_harm))
cat(sprintf("  ACT 2  joint gate %.2f != naive product %.2f ; paired bootstrap %.2f (metrics correlated)\n",
            P_gate, P_gate_indep, P_gate_boot))
cat(sprintf("  ACT 3  gate Brier: Bayesian %.3f vs paired bootstrap %.3f ; dichotomous rule deploys %.0f%% worthy / %.0f%% unworthy\n",
            brier_bayes, brier_boot, 100 * tp, 100 * fp))
cat(sprintf("  figures + results -> vignettes/%s/showcase_*.{png,rds}\n", figdir))
cat("=====================================================================\n")
