# temporal_validation.R
# ---------------------------------------------------------------------------
# Part D companion -- temporal versus independent external validation.
#
# The headline worked example trains on all of Rotterdam and validates on GBSG,
# where the calibration slope falls below one. That single comparison cannot say
# whether the miscalibration is a property of the MODEL or of the TRANSPORT to a
# different population. This script separates the two by training ONE model and
# validating it on TWO cohorts:
#
#   temporal      Rotterdam 1988-1993  -- same institution, later accrual era
#   independent   GBSG                 -- different population entirely
#
# The training set is the earlier Rotterdam era (1978-1987), so both validation
# cohorts are out-of-sample for the same fitted model. Any difference in
# calibration is then attributable to the cohort, not to a different fit.
#
# The temporal cohort also carries roughly three times the events of GBSG, so it
# is the better-powered of the two validations.
#
#   Rscript vignettes/temporal_validation.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ---------------------------------------------------------------------------

suppressMessages(library(survival))
source(file.path("..", "R", "load_all.R"), chdir = TRUE)
YEAR <- 365.25

## ---- Data preparation: harmonise Rotterdam & GBSG -------------------------
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

# Split on accrual year. 1987 is the cut that keeps ~100 events per model
# degree of freedom in training while maximising validation events.
SPLIT_YEAR <- 1987
rot_early <- dplyr::filter(rot, year <= SPLIT_YEAR)
rot_late  <- dplyr::filter(rot, year >  SPLIT_YEAR)

train <- prep(rot_early, "rfstime", "rfs") |> tidyr::drop_na()
temporal <- prep(rot_late, "rfstime", "rfs") |> tidyr::drop_na()
independent <- prep(gbsg, "rfstime", "status") |> tidyr::drop_na()

horizons <- c(2, 5)                              # years; headline at 5
h5 <- which(horizons == 5)
tau <- 5

cat(sprintf("Train (Rotterdam <=%d):  n=%d, %d events\n", SPLIT_YEAR,
            nrow(train), sum(train$status)))
cat(sprintf("Temporal (Rotterdam >%d): n=%d, %d events, %.0f%% censored\n", SPLIT_YEAR,
            nrow(temporal), sum(temporal$status), 100 * mean(temporal$status == 0)))
cat(sprintf("Independent (GBSG):       n=%d, %d events, %.0f%% censored\n",
            nrow(independent), sum(independent$status),
            100 * mean(independent$status == 0)))

form <- Surv(time, status) ~ age + nodes + grade3 + size + logpgr + hormon + meno
fit  <- coxph(form, data = train, x = TRUE, y = TRUE)

## ---- Bayesian validation suite (B2 + B3) at the 5-year horizon ------------
S <- 500L
uno   <- function(d, w) metric_uno_c(d, w, tau = tau)
slope <- metric_cal_slope
ici   <- function(d, w) metric_ici(d, w, horizon = h5, method = "auto")
citl  <- function(d, w) metric_citl(d, w)[paste0("CITL_t", h5)]

suite <- list(harrell = metric_harrell_c, uno = uno,
              citl = citl, slope = slope, ici = ici)

relabel <- function(m) {
  m <- sub("^[^.]+\\.", "", m)
  m <- sub("_t1\\b", "@2y", m)
  sub("_t2\\b", "@5y", m)
}

# Same fitted model, propagated through each validation cohort in turn. Keep the
# raw draws: the slope's posterior mass below one is the quantity of interest,
# and a credible interval that straddles one does not by itself settle it.
validate <- function(newdata, label) {
  ps <- psurv_external(fit, newdata, time = "time", status = "status",
                       horizons = horizons, S = S, seed = 2026)
  post <- propagate_many(ps, suite, weights = "bootstrap", seed = 1)
  tab <- posterior_table(post) |>
    dplyr::mutate(metric = relabel(metric), cohort = label) |>
    dplyr::select(cohort, metric, summary)
  draws <- unclass(post)                          # [S x metrics], suite-prefixed
  list(tab = tab, slope = as.numeric(draws[, grep("^slope\\.", colnames(draws))]))
}

cat("\nPropagating posterior (parameter + Bayesian bootstrap)...\n")
v_temp <- validate(temporal,    "temporal (Rotterdam >1987)")
v_indep <- validate(independent, "independent (GBSG)")

cat("\n== Posterior median [95% CrI], one model, two validation cohorts ==\n")
print(as.data.frame(dplyr::bind_rows(v_temp$tab, v_indep$tab)), row.names = FALSE)

## ---- Posterior probability of calibration drift ---------------------------
# P(slope < 1) is "risk predictions are too extreme"; P(slope < 0.9) adds a
# margin so that trivially small drift does not count as drift.
drift <- function(sl) sprintf("P(slope<1) = %.3f   P(slope<0.9) = %.3f",
                              mean(sl < 1), mean(sl < 0.9))
cat("\n== Posterior probability of calibration drift ==\n")
cat("temporal    (Rotterdam >1987):", drift(v_temp$slope), "\n")
cat("independent (GBSG):           ", drift(v_indep$slope), "\n")

## ---- Frequentist calibration slope on each cohort, for reference ----------
freq_slope <- function(newdata) {
  ps <- psurv_external(fit, newdata, time = "time", status = "status",
                       horizons = horizons, S = S, seed = 2026)
  lp <- colMeans(ps$eta)
  cs <- coxph(Surv(time, status) ~ lp,
              data = tibble::tibble(time = newdata$time, status = newdata$status))
  b  <- unname(coef(cs)); se <- sqrt(vcov(cs)[1, 1])
  sprintf("%.3f [%.3f, %.3f]", b, b - 1.96 * se, b + 1.96 * se)
}

cat("\n== Frequentist calibration slope [95% CI] ==\n")
cat("temporal    (Rotterdam >1987):", freq_slope(temporal), "\n")
cat("independent (GBSG):           ", freq_slope(independent), "\n")
