# worked_example.R
# ---------------------------------------------------------------------------
# Part D — worked example on public breast-cancer data. Train a survival model
# on the Rotterdam cohort and EXTERNALLY validate it on the German Breast Cancer
# Study Group (GBSG) cohort at a 5-year horizon, producing the full Bayesian
# validation suite (posterior median [95% CrI]) alongside the frequentist
# estimates, a posterior-predictive calibration plot, and posterior decision
# curves. Runs on the coxph-Laplace backend (no Stan toolchain required).
#
#   Rscript vignettes/worked_example.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ---------------------------------------------------------------------------

suppressMessages(library(survival))
source(file.path("..", "R", "load_all.R"), chdir = TRUE)
YEAR <- 365.25
figdir <- "figures"; dir.create(figdir, showWarnings = FALSE)

## ---- Data preparation: harmonise Rotterdam (train) & GBSG (test) ----------
prep <- function(df, time, status) {
  size_cat <- if (is.factor(df$size)) df$size else
    cut(df$size, breaks = c(-Inf, 20, 50, Inf), labels = c("<=20", "20-50", ">50"))
  tibble::tibble(
    time   = df[[time]] / YEAR,                 # days -> years
    status = df[[status]],
    age    = df$age,
    nodes  = df$nodes,
    grade3 = as.integer(df$grade >= 3),          # grade 3 vs lower
    size   = factor(size_cat, levels = c("<=20", "20-50", ">50")),
    logpgr = log(df$pgr + 1),
    hormon = df$hormon,
    meno   = df$meno)
}

rot <- rotterdam |>                             # recurrence-free survival
  dplyr::mutate(rfs     = as.integer(recur == 1 | death == 1),
                rfstime = ifelse(recur == 1, rtime, dtime))
train <- prep(rot, "rfstime", "rfs")  |> tidyr::drop_na()
test  <- prep(gbsg, "rfstime", "status") |> tidyr::drop_na()

horizons <- c(2, 5)                              # years; headline at 5
h5 <- which(horizons == 5)
tau <- 5
cat(sprintf("Train (Rotterdam): n=%d, %d events\nTest (GBSG): n=%d, %d events, %.0f%% censored\n",
            nrow(train), sum(train$status),
            nrow(test), sum(test$status), 100 * mean(test$status == 0)))

form <- Surv(time, status) ~ age + nodes + grade3 + size + logpgr + hormon + meno
fit  <- coxph(form, data = train, x = TRUE, y = TRUE)

## ---- External-validation posterior survival on GBSG (B1) ------------------
S <- 500L
ps <- psurv_external(fit, test, time = "time", status = "status",
                     horizons = horizons, S = S, seed = 2026)
print(ps)

## ---- Bayesian validation suite (B2 + B3) at the 5-year horizon ------------
uno   <- function(d, w) metric_uno_c(d, w, tau = tau)
har   <- metric_harrell_c
gonen <- metric_gonen_k
citl  <- function(d, w) metric_citl(d, w)[paste0("CITL_t", h5)]
slope <- metric_cal_slope
ici   <- function(d, w) metric_ici(d, w, horizon = h5, method = "auto")
brier <- function(d, w) metric_brier(d, w)[c(paste0("Brier_t", h5), "IBS")]
ipa   <- function(d, w) metric_ipa(d, w)[paste0("IPA_t", h5)]
r2d   <- metric_r2d

suite <- list(harrell = har, uno = uno, gonen = gonen, citl = citl,
              slope = slope, ici = ici, brier = brier, ipa = ipa, r2d = r2d)

cat("\nPropagating posterior (parameter + Bayesian bootstrap)...\n")
post_bb    <- propagate_many(ps, suite, weights = "bootstrap", seed = 1)
post_param <- propagate_many(ps, suite, weights = "none", seed = 1)

tab_bb <- posterior_table(post_bb)
# Strip the suite-list prefix ("harrell.C_harrell" -> "C_harrell") and relabel
# horizon indices to years (t1 -> 2y, t2 -> 5y) for readable output.
relabel <- function(m) {
  m <- sub("^[^.]+\\.", "", m)
  m <- sub("_t1\\b", "@2y", m); m <- sub("_t2\\b", "@5y", m)
  m
}
tab_bb <- dplyr::mutate(tab_bb, metric = relabel(metric))
cat("\n== Bayesian posterior median [95% CrI] (param + bootstrap) ==\n")
print(as.data.frame(dplyr::select(tab_bb, metric, summary)), row.names = FALSE)

## ---- Frequentist comparators (point estimate [95% CI]) --------------------
lp   <- colMeans(ps$eta)                          # posterior-mean linear predictor
S5   <- colMeans(ps$surv[, , h5])                 # posterior-mean 5y survival
dfv  <- tibble::tibble(time = test$time, status = test$status)

freq_wald <- function(est, se) c(estimate = est, lower = est - 1.96 * se,
                                 upper = est + 1.96 * se)
# Harrell & Uno C (concordance gives var of the reported statistic).
ch <- concordance(Surv(time, status) ~ lp, data = dfv)
cu <- concordance(Surv(time, status) ~ lp, data = dfv, timewt = "n/G2", ymax = tau)
fr_harrell <- freq_wald(1 - ch$concordance, sqrt(ch$var))
fr_uno     <- freq_wald(1 - cu$concordance, sqrt(cu$var))
# Calibration slope: Cox on lp.
cs <- coxph(Surv(time, status) ~ lp, data = dfv)
fr_slope <- freq_wald(unname(coef(cs)), sqrt(vcov(cs)[1, 1]))
# IBS & ICI: point at posterior mean + classical bootstrap CI.
mean_draw <- list(surv2d = matrix(colMeans(ps$surv), nrow = ps$n),
                  eta = lp, time = test$time, status = test$status,
                  horizons = horizons, n = ps$n)
boot_ci <- function(metric_fn, B = 500L, seed = 7) {
  set.seed(seed); n <- ps$n
  v <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    db <- list(surv2d = mean_draw$surv2d[idx, , drop = FALSE], eta = lp[idx],
               time = test$time[idx], status = test$status[idx],
               horizons = horizons, n = n)
    as.numeric(metric_fn(db, rep(1 / n, n))[1])
  })
  c(estimate = as.numeric(metric_fn(mean_draw, rep(1 / n, n))[1]),
    stats::quantile(v, c(0.025, 0.975), names = FALSE))
}
fr_ibs <- boot_ci(function(d, w) metric_brier(d, w)["IBS"])
fr_ici <- boot_ci(function(d, w) metric_ici(d, w, horizon = h5, method = "auto"))
names(fr_ibs) <- names(fr_uno); names(fr_ici) <- names(fr_uno)

freq <- tibble::tibble(
  metric = c("C_harrell", "C_uno", "cal_slope", "IBS", "ICI"),
  estimate = c(fr_harrell[1], fr_uno[1], fr_slope[1], fr_ibs[1], fr_ici[1]),
  lower    = c(fr_harrell[2], fr_uno[2], fr_slope[2], fr_ibs[2], fr_ici[2]),
  upper    = c(fr_harrell[3], fr_uno[3], fr_slope[3], fr_ibs[3], fr_ici[3]))
# freq metric names (C_harrell, C_uno, cal_slope, IBS, ICI) already match the
# relabelled Bayesian names, so posterior_vs_frequentist() merges cleanly.

## ---- Headline side-by-side table ------------------------------------------
headline <- posterior_vs_frequentist(tab_bb, freq)
cat("\n== Headline: Bayesian [95% CrI] vs frequentist [95% CI] ==\n")
print(as.data.frame(headline), row.names = FALSE)

## ---- PPC calibration (B4) and posterior decision curves (B6) --------------
## Computed here, drawn by `worked_example_composite.R`, which faces both against
## the full-Bayes backend in one faceted figure. Both objects go into the results
## file below; a standalone PNG from this script would duplicate half of that
## figure. `plot_calibration_band()` / `plot_net_benefit()` draw the single-backend
## version, and worked_example.Rmd renders them inline.
ppc <- ppc_calibration(ps, horizon = h5, method = "auto", draws = 200L, seed = 3)
pnb <- posterior_net_benefit(ps, thresholds = seq(0.1, 0.7, 0.05),
                             horizon = h5, weights = "bootstrap",
                             draws = 300L, seed = 4)

## ---- Optimism (B5, internal on Rotterdam) ---------------------------------
fitter <- make_psurv_fitter(form, horizons, S = 150L)
opt <- bb_optimism(fitter, train, metric_harrell_c, B = 20L, seed = 5)
cat(sprintf("\n== Internal optimism (Rotterdam, Harrell C) ==\napparent %.3f, optimism %.3f, corrected %.3f\n",
            opt$apparent, opt$optimism$mean, opt$corrected))

saveRDS(list(headline = headline, tab_bb = tab_bb, freq = freq, ppc = ppc,
             pnb = pnb, optimism = opt),
        file.path(figdir, "worked_example_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(figdir, "session_info.txt"))
cat(sprintf("\nFigures + results written to vignettes/%s/\n", figdir))
