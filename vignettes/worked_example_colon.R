# worked_example_colon.R
# ---------------------------------------------------------------------------
# Part D — extra worked example on the colon adjuvant-trial cohort
# (survival::colon, death endpoint).
#
# A larger cohort with rich covariates and a natural case-mix / transport shift:
# train on the observation and levamisole arms, then EXTERNALLY validate on the
# levamisole+5-FU arm. Because that arm receives the effective therapy, its
# baseline risk shifts, so calibration-in-the-large and the calibration slope
# reveal the transport cost that discrimination alone hides. Produces the full
# Bayesian validation suite (posterior median [95% CrI]) alongside the
# frequentist estimates, the posterior calibration / decision-curve objects, and
# a joint deployment-gate probability. Runs on the coxph-Laplace backend (no
# Stan toolchain required); reproducible with the repo seed.
#
#   Rscript vignettes/worked_example_colon.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ---------------------------------------------------------------------------

suppressMessages({library(survival); library(dplyr); library(tidyr)})
source(file.path("..", "R", "load_all.R"), chdir = TRUE)
YEAR <- 365.25
figdir <- "figures"; dir.create(figdir, showWarnings = FALSE)

## ---- Data preparation: colon adjuvant trial, death endpoint ---------------
# colon has two rows per patient (recurrence + death); keep the death endpoint.
# status is already 0/1. Complete-case on the retained covariates (a handful of
# NA in nodes).
colon_cc <- colon |>
  filter(etype == 2) |>
  transmute(rx,
            time   = time / YEAR,                 # days -> years
            status,
            age,
            sex      = factor(sex, levels = c(0, 1), labels = c("female", "male")),
            obstruct, adhere, nodes,
            extent, surg, node4) |>
  drop_na()

## ---- External-validation split: treatment arm (case-mix / transport shift) -
# Develop on the observation + levamisole arms; validate on the levamisole+5-FU
# arm, whose effective therapy lowers baseline risk. See Runbook B.
train <- colon_cc |> filter(rx %in% c("Obs", "Lev")) |> select(-rx)
test  <- colon_cc |> filter(rx == "Lev+5FU")        |> select(-rx)

horizons <- c(3, 5)                               # years; headline at 5 years
hH  <- 2L                                         # headline horizon index
tau <- 5                                          # Uno-C truncation (< last event)
cat(sprintf(paste0("Colon adjuvant-trial external validation (arm split)\n",
                   "  Train (Obs + Lev): n=%d, %d events\n",
                   "  Test  (Lev+5FU):   n=%d, %d events, %.0f%% censored\n"),
            nrow(train), sum(train$status),
            nrow(test), sum(test$status), 100 * mean(test$status == 0)))

form <- Surv(time, status) ~ age + sex + obstruct + adhere + nodes +
  extent + surg + node4
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
ici   <- function(d, w) metric_ici(d, w, horizon = hH, method = "auto")
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
  m <- sub("_t1\\b", "@3y", m); m <- sub("_t2\\b", "@5y", m)
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
fr_ici <- boot_ci(function(d, w) metric_ici(d, w, horizon = hH, method = "auto")["ICI"])
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

## ---- Joint-posterior scatter of the two gate-critical axes (B2 + B3) -------
# The joint gate sits well above the product of the three marginal
# pass-probabilities because the metrics are dependent across the posterior:
# Harrell C and the calibration slope move together, while the slope and net
# benefit are close to independent. The C-slope plane therefore shows directly
# why separate marginal intervals cannot recover the joint event. Every point
# carries the full three-criterion gate (a passing draw also clears dNB > 0).
# Style matches the house joint-gate figure (showcase_joint.png).
rho_Cslope <- cor(Jm[, "C_harrell"], Jm[, "cal_slope"])
cat(sprintf("  corr(C, slope) = %.2f ; corr(slope, dNB) = %.2f\n",
            rho_Cslope, cor(Jm[, "cal_slope"], Jm[, "dNB"])))
# Legend uses the manuscript's incremental-net-benefit notation (Delta NB).
gate_lab <- c("gate passes (also ΔNB > 0)", "gate fails")
joint_df <- tibble::tibble(
  cal_slope = Jm[, "cal_slope"], C_harrell = Jm[, "C_harrell"],
  cohort = "Colon adjuvant trial (Lev+5FU arm)",
  gate = factor(dplyr::if_else(
    Jm[, "C_harrell"] >= C_MIN & Jm[, "cal_slope"] >= SLOPE_LO &
      Jm[, "cal_slope"] <= SLOPE_HI & Jm[, "dNB"] > NB_TOL,
    gate_lab[1], gate_lab[2]), levels = gate_lab)) |>
  # passing draws are the minority; draw them last so they are not painted over.
  dplyr::arrange(gate == gate_lab[1])
p_joint <- ggplot2::ggplot(joint_df,
                           ggplot2::aes(cal_slope, C_harrell, fill = gate)) +
  ggplot2::geom_vline(xintercept = c(SLOPE_LO, SLOPE_HI), linetype = 2,
                      colour = "gray40") +
  ggplot2::geom_hline(yintercept = C_MIN, linetype = 2, colour = "gray40") +
  ggplot2::geom_point(shape = 21, colour = "gray40", stroke = 0.3,
                      size = 1.4, alpha = 0.7) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::scale_fill_manual(
    values = stats::setNames(bc_palette[c("blue", "orange")], gate_lab)) +
  ggplot2::scale_x_continuous(breaks = number_ticks(5)) +
  ggplot2::scale_y_continuous(breaks = number_ticks(5)) +
  ggplot2::labs(x = "Calibration slope", y = "Harrell C", fill = NULL) +
  theme_bayesval()
bc_ggsave(file.path("..", "manuscript", "tex", "figures", "colon_joint.png"),
          p_joint, width = 4.9, height = 4.3)
saveRDS(list(
  draws  = tibble::tibble(cal_slope = Jm[, "cal_slope"],
                          C_harrell = Jm[, "C_harrell"], dNB = Jm[, "dNB"]),
  gate   = list(C_min = C_MIN, slope_lo = SLOPE_LO, slope_hi = SLOPE_HI),
  P_gate = P_gate, P_gate_indep = P_gate_indep, rho_Cslope = rho_Cslope,
  label  = "Colon adjuvant trial", n = nrow(test), events = sum(test$status)),
  file.path(figdir, "colon_joint_draws.rds"))

## ---- PPC calibration (B4) and posterior decision curves (B6) --------------
ppc <- tryCatch(
  ppc_calibration(ps, horizon = hH, method = "auto", draws = 200L, seed = 3),
  error = function(e) {message("PPC calibration skipped: ", conditionMessage(e)); NULL})
pnb <- posterior_net_benefit(ps, thresholds = seq(0.1, 0.7, 0.05),
                             horizon = hH, weights = "bootstrap",
                             draws = 300L, seed = 4)
if (!is.null(ppc))
  bc_ggsave(file.path(figdir, "colon_calibration.png"),
            plot_calibration_band(ppc), width = 4.4, height = 4.2)
bc_ggsave(file.path(figdir, "colon_decision_curve.png"),
          plot_net_benefit(pnb), width = 4.8, height = 4.2)

saveRDS(list(headline = headline, tab_bb = tab_bb, freq = freq,
             ppc = ppc, pnb = pnb,
             gate = c(disc = P_disc, cal = P_cal, bene = P_bene, harm = P_harm,
                      gate = P_gate, gate_indep = P_gate_indep),
             n = list(train = nrow(train), test = nrow(test),
                      ev_train = sum(train$status), ev_test = sum(test$status),
                      cens_test = mean(test$status == 0))),
        file.path(figdir, "colon_results.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(figdir, "colon_session_info.txt"))
cat(sprintf("\nResults + figures written to vignettes/%s/ (colon_*)\n", figdir))
