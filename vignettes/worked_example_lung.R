# worked_example_lung.R
# ---------------------------------------------------------------------------
# Part D — extra worked example on the NCCTG lung cohort (survival::lung).
#
# A small, event-heavy cohort with missing covariates: exactly the setting the
# manuscript's small-cohort deployment gate is built for, and it links to the
# multiple-imputation section (complete-case here; the MI variants are optional).
# Train a Cox PH model on the largest recruiting institutes and EXTERNALLY
# validate it on the remaining institutes (a centre split), producing the full
# Bayesian validation suite (posterior median [95% CrI]) alongside the
# frequentist estimates, the posterior calibration / decision-curve objects, and
# a joint deployment-gate probability. Runs on the coxph-Laplace backend (no
# Stan toolchain required); reproducible with the repo seed.
#
#   Rscript vignettes/worked_example_lung.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ---------------------------------------------------------------------------

suppressMessages({library(survival); library(dplyr); library(tidyr)})
source(file.path("..", "R", "load_all.R"), chdir = TRUE)
YEAR <- 365.25
figdir <- "figures"; dir.create(figdir, showWarnings = FALSE)

## ---- Data preparation: NCCTG lung -----------------------------------------
# status is coded 1 = censored, 2 = dead; recode to the framework's 0/1. Drop
# meal.cal (too many NA) and complete-case on the retained covariates.
lung_cc <- lung |>
  transmute(inst,
            time   = time / YEAR,                 # days -> years
            status = status - 1L,                 # 1/2 -> 0/1
            age,
            sex    = factor(sex, levels = c(1, 2), labels = c("male", "female")),
            ph.ecog, ph.karno, wt.loss) |>
  drop_na()

## ---- External-validation split: largest institutes vs the rest -------------
# Centre (geographic) split: the three largest recruiting institutes develop the
# model, the remaining institutes validate it. This keeps the validation cohort
# genuinely independent of development and leaves >40 events (no random-split
# fallback needed); see Runbook B.
inst_by_size <- as.integer(names(sort(table(lung_cc$inst), decreasing = TRUE)))
dev_inst <- inst_by_size[1:3]
train <- lung_cc |> filter(inst %in% dev_inst)  |> select(-inst)
test  <- lung_cc |> filter(!inst %in% dev_inst) |> select(-inst)

horizons <- c(1, 2)                               # years; headline at 1 year
hH  <- 1L                                         # headline horizon index
tau <- 2                                          # Uno-C truncation (< last event)
cat(sprintf(paste0("NCCTG lung external validation (centre split)\n",
                   "  Train (institutes %s): n=%d, %d events\n",
                   "  Test  (other institutes): n=%d, %d events, %.0f%% censored\n"),
            paste(sort(dev_inst), collapse = ", "),
            nrow(train), sum(train$status),
            nrow(test), sum(test$status), 100 * mean(test$status == 0)))

form <- Surv(time, status) ~ age + sex + ph.ecog + ph.karno + wt.loss
fit  <- coxph(form, data = train, x = TRUE, y = TRUE)

## ---- External-validation posterior survival (B1) --------------------------
S <- 500L
ps <- psurv_external(fit, test, time = "time", status = "status",
                     horizons = horizons, S = S, seed = 1)
print(ps)

## ---- Deployment gate + net-benefit metric ---------------------------------
PT       <- 0.20                                  # clinical risk threshold
C_MIN    <- 0.65
SLOPE_LO <- 0.80
SLOPE_HI <- 1.25
NB_TOL   <- 1e-8
metric_dnb <- function(pt, h) function(draw, w) {
  nb <- net_benefit_draw(draw, thresholds = pt, horizon = h, w = w)
  c(NB_model = nb$model, dNB = nb$model - max(nb$all, 0))
}

## ---- Bayesian validation suite (B2 + B3) at the headline horizon ----------
uno   <- function(d, w) metric_uno_c(d, w, tau = tau)
citl  <- function(d, w) metric_citl(d, w)[paste0("CITL_t", hH)]
# Grouped Kaplan-Meier calibration (survival-only). The cohort's short survival
# packs 1-year risks into the upper range, where the adaptive-hazard smoother is
# unstable; the grouped estimator is the documented fallback and keeps the run
# reproducible (no C-level solver).
ici   <- function(d, w) metric_ici(d, w, horizon = hH, method = "group")
brier <- function(d, w) metric_brier(d, w)[c(paste0("Brier_t", hH), "IBS")]
ipa   <- function(d, w) metric_ipa(d, w)[paste0("IPA_t", hH)]

suite <- list(harrell = metric_harrell_c, uno = uno, gonen = metric_gonen_k,
              citl = citl, slope = metric_cal_slope, ici = ici,
              brier = brier, ipa = ipa, r2d = metric_r2d,
              nb = metric_dnb(PT, hH))

cat("\nPropagating posterior (parameter + Bayesian bootstrap)...\n")
post_bb    <- propagate_many(ps, suite, weights = "bootstrap", seed = 1)
post_param <- propagate_many(ps, suite, weights = "none", seed = 1)

# Strip the suite-list prefix and relabel horizon indices to years.
relabel <- function(m) {
  m <- sub("^[^.]+\\.", "", m)
  m <- sub("_t1\\b", "@1y", m); m <- sub("_t2\\b", "@2y", m)
  m
}
tab_bb <- posterior_table(post_bb) |> dplyr::mutate(metric = relabel(metric))
cat("\n== Bayesian posterior median [95% CrI] (param + bootstrap) ==\n")
print(as.data.frame(dplyr::select(tab_bb, metric, summary)), row.names = FALSE)

## ---- Frequentist comparators (point estimate [95% CI]) --------------------
lp  <- colMeans(ps$eta)
dfv <- tibble::tibble(time = test$time, status = test$status)
freq_wald <- function(est, se) c(estimate = est, lower = est - 1.96 * se,
                                 upper = est + 1.96 * se)
ch <- concordance(Surv(time, status) ~ lp, data = dfv)
cu <- concordance(Surv(time, status) ~ lp, data = dfv, timewt = "n/G2", ymax = tau)
fr_harrell <- freq_wald(1 - ch$concordance, sqrt(ch$var))
fr_uno     <- freq_wald(1 - cu$concordance, sqrt(cu$var))
cs <- coxph(Surv(time, status) ~ lp, data = dfv)
fr_slope <- freq_wald(unname(coef(cs)), sqrt(vcov(cs)[1, 1]))

mean_draw <- list(surv2d = matrix(colMeans(ps$surv), nrow = ps$n),
                  eta = lp, time = test$time, status = test$status,
                  horizons = horizons, n = ps$n, H = length(horizons))
boot_ci <- function(metric_fn, B = 500L, seed = 7) {
  set.seed(seed); n <- ps$n
  v <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    db <- list(surv2d = mean_draw$surv2d[idx, , drop = FALSE], eta = lp[idx],
               time = test$time[idx], status = test$status[idx],
               horizons = horizons, n = n, H = length(horizons))
    as.numeric(metric_fn(db, rep(1 / n, n))[1])
  })
  c(estimate = as.numeric(metric_fn(mean_draw, rep(1 / n, n))[1]),
    stats::quantile(v, c(0.025, 0.975), names = FALSE))
}
fr_ibs <- boot_ci(function(d, w) metric_brier(d, w)["IBS"])
fr_ici <- boot_ci(function(d, w) metric_ici(d, w, horizon = hH, method = "group")["ICI"])
names(fr_ibs) <- names(fr_uno); names(fr_ici) <- names(fr_uno)

freq <- tibble::tibble(
  metric   = c("C_harrell", "C_uno", "cal_slope", "IBS", "ICI"),
  estimate = c(fr_harrell[1], fr_uno[1], fr_slope[1], fr_ibs[1], fr_ici[1]),
  lower    = c(fr_harrell[2], fr_uno[2], fr_slope[2], fr_ibs[2], fr_ici[2]),
  upper    = c(fr_harrell[3], fr_uno[3], fr_slope[3], fr_ibs[3], fr_ici[3]))

headline <- posterior_vs_frequentist(tab_bb, freq)
cat("\n== Headline: Bayesian [95% CrI] vs frequentist [95% CI] ==\n")
print(as.data.frame(headline), row.names = FALSE)

## ---- Deployment-gate probabilities (B6) -----------------------------------
Jm <- unclass(post_bb); colnames(Jm) <- sub("^[^.]+\\.", "", colnames(Jm))
P_disc <- mean(Jm[, "C_harrell"] >= C_MIN)
P_cal  <- mean(Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI)
P_bene <- mean(Jm[, "dNB"] >  NB_TOL)
P_harm <- mean(Jm[, "dNB"] < -NB_TOL)
P_gate <- mean(Jm[, "C_harrell"] >= C_MIN &
               Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI &
               Jm[, "dNB"] > NB_TOL)
P_gate_indep <- P_disc * P_cal * P_bene
cat("\n== Deployment gate (C>=0.65, slope in [0.80,1.25], dNB>0 @20%) ==\n")
cat(sprintf("  P(discrimination adequate)  = %.2f\n", P_disc))
cat(sprintf("  P(calibration acceptable)   = %.2f\n", P_cal))
cat(sprintf("  P(net clinical benefit)     = %.2f\n", P_bene))
cat(sprintf("  P(net harm)                 = %.2f\n", P_harm))
cat(sprintf("  P(JOINT gate passes)        = %.2f  (naive product %.2f)\n",
            P_gate, P_gate_indep))

## ---- PPC calibration (B4) and posterior decision curves (B6) --------------
# This cohort's short survival puts most 1-year risks in the upper range, where a
# smoothed observed calibration curve is barely estimable. Use the grouped
# Kaplan-Meier curve (survival-only, no C-level solver) and guard the plot so a
# non-estimable curve never aborts the run; the suite ICI above keeps the primary
# adaptive-hazard smoother.
ppc <- tryCatch(
  ppc_calibration(ps, horizon = hH, method = "group", draws = 200L, seed = 3),
  error = function(e) {message("PPC calibration skipped: ", conditionMessage(e)); NULL})
pnb <- posterior_net_benefit(ps, thresholds = seq(0.1, 0.7, 0.05),
                             horizon = hH, weights = "bootstrap",
                             draws = 300L, seed = 4)
if (!is.null(ppc))
  bc_ggsave(file.path(figdir, "lung_calibration.png"),
            plot_calibration_band(ppc), width = 4.4, height = 4.2)
bc_ggsave(file.path(figdir, "lung_decision_curve.png"),
          plot_net_benefit(pnb), width = 4.8, height = 4.2)

saveRDS(list(headline = headline, tab_bb = tab_bb, freq = freq,
             ppc = ppc, pnb = pnb,
             gate = c(disc = P_disc, cal = P_cal, bene = P_bene, harm = P_harm,
                      gate = P_gate, gate_indep = P_gate_indep),
             n = list(train = nrow(train), test = nrow(test),
                      ev_train = sum(train$status), ev_test = sum(test$status),
                      cens_test = mean(test$status == 0), dev_inst = sort(dev_inst))),
        file.path(figdir, "lung_results.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(figdir, "lung_session_info.txt"))
cat(sprintf("\nResults + figures written to vignettes/%s/ (lung_*)\n", figdir))
